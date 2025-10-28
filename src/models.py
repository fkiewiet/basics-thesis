# src/models.py
from __future__ import annotations
import math
import torch
import torch.nn as nn
import torch.nn.functional as F


# ----------------------------- Spectral Convolution -----------------------------
class SpectralConv2d(nn.Module):
    """
    2D spectral convolution on low Fourier modes (FNO-style).
    Expects real input of shape (B, C_in, H, W). Internally uses rFFT:
      X = rfft2(x) -> shape (B, C_in, H, W//2+1), complex
    Learns complex weights on a rectangular band of modes in both the
    'top' and 'bottom' vertical bands to respect rFFT symmetry.

    Args
    ----
    in_ch, out_ch : int
    modes_h, modes_w : int
        Number of low-frequency modes to keep along H and W.
        They will be clamped to available sizes at runtime.
    """

    def __init__(self, in_ch: int, out_ch: int, modes_h: int = 12, modes_w: int = 12):
        super().__init__()
        self.in_ch = in_ch
        self.out_ch = out_ch
        self.mh = int(modes_h)
        self.mw = int(modes_w)

        # Complex weights (real+imag parts stored separately)
        # We need two vertical bands: top ([:mh]) and bottom ([-mh:])
        # Weight shapes: (in_ch, out_ch, mh, mw, 2)
        scale = 1.0 / math.sqrt(in_ch * out_ch)
        self.weight_top = nn.Parameter(scale * torch.randn(in_ch, out_ch, self.mh, self.mw, 2))
        self.weight_bot = nn.Parameter(scale * torch.randn(in_ch, out_ch, self.mh, self.mw, 2))

    def compl(self, w: torch.Tensor) -> torch.Tensor:
        """Pack real+imag (..,2) into complex tensor (...)."""
        return torch.view_as_complex(w)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        B, C, H, W = x.shape
        assert C == self.in_ch, f"SpectralConv2d: expected {self.in_ch} in-ch, got {C}"

        X = torch.fft.rfft2(x, norm="ortho")  # (B, C, H, W//2+1), complex
        out_ft = torch.zeros(B, self.out_ch, H, W // 2 + 1, dtype=X.dtype, device=X.device)

        # Determine how many modes we can actually use
        mh = min(self.mh, H)
        mw = min(self.mw, W // 2 + 1)

        if mh > 0 and mw > 0:
            Wt = self.compl(self.weight_top[:, :, :mh, :mw, :])   # (Cin, Cout, mh, mw)
            Wb = self.compl(self.weight_bot[:, :, :mh, :mw, :])   # (Cin, Cout, mh, mw)

            # Top band
            x_top = X[:, :, :mh, :mw]                       # (B, Cin, mh, mw)
            out_ft[:, :, :mh, :mw] = torch.einsum(          # → (B, Cout, mh, mw)
                "bcij,coij->boij",                          # c,i,j = Cin,mh,mw
                x_top, Wt
            )

            # Bottom band
            if H - mh > 0:
                x_bot = X[:, :, -mh:, :mw]                  # (B, Cin, mh, mw)
                out_ft[:, :, -mh:, :mw] = torch.einsum(     # → (B, Cout, mh, mw)
                    "bcij,coij->boij",
                    x_bot, Wb
                )

        # Back to spatial
        y = torch.fft.irfft2(out_ft, s=(H, W), norm="ortho")  # (B, Cout, H, W), real
        return y


# ----------------------------- FNO Block & Model -----------------------------
class FNOBlock(nn.Module):
    """
    One residual FNO layer: spectral conv + 1x1 conv skip + activation.
    """
    def __init__(self, width: int, modes_h: int, modes_w: int, act=nn.GELU()):
        super().__init__()
        self.spectral = SpectralConv2d(width, width, modes_h, modes_w)
        self.pointwise = nn.Conv2d(width, width, kernel_size=1)
        self.act = act

        # lightweight init
        nn.init.kaiming_uniform_(self.pointwise.weight, a=math.sqrt(5))
        if self.pointwise.bias is not None:
            fan_in, _ = nn.init._calculate_fan_in_and_fan_out(self.pointwise.weight)
            bound = 1 / math.sqrt(fan_in)
            nn.init.uniform_(self.pointwise.bias, -bound, bound)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.act(self.spectral(x) + self.pointwise(x))


class SimpleFNO(nn.Module):
    """
    Minimal, stable FNO for (B, C_in, H, W) -> (B, C_out, H, W)

    Args
    ----
    in_ch, out_ch : channels
    width : hidden channel width used throughout
    modes : (modes_h, modes_w)
    layers : number of FNO residual blocks
    """

    def __init__(self, in_ch: int = 3, out_ch: int = 2, width: int = 48,
                 modes: tuple[int, int] = (12, 12), layers: int = 4, act: nn.Module = nn.GELU()):
        super().__init__()
        mh, mw = modes
        self.fc0 = nn.Conv2d(in_ch, width, kernel_size=1)
        self.blocks = nn.ModuleList([FNOBlock(width, mh, mw, act=act) for _ in range(layers)])
        self.fc1 = nn.Conv2d(width, max(32, width // 2), kernel_size=1)
        self.fc2 = nn.Conv2d(max(32, width // 2), out_ch, kernel_size=1)
        self.act = act

        # Kaiming init
        for m in [self.fc0, self.fc1, self.fc2]:
            nn.init.kaiming_uniform_(m.weight, a=math.sqrt(5))
            if m.bias is not None:
                fan_in, _ = nn.init._calculate_fan_in_and_fan_out(m.weight)
                bound = 1 / math.sqrt(fan_in)
                nn.init.uniform_(m.bias, -bound, bound)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        x = self.fc0(x)
        for blk in self.blocks:
            x = blk(x)
        x = self.act(self.fc1(x))
        x = self.fc2(x)
        return x


# ----------------------------- Local CNN (operator-like) -----------------------------
class LocalCNN(nn.Module):
    """
    Compact local CNN suitable for learning local operators/preconditioners.
    Keeps your original signature and behavior, but adds residual path and
    (optional) normalization via GroupNorm for stability on small batches.

    Args
    ----
    in_ch, out_ch : channels
    width : hidden width
    ksize : kernel size (odd)
    layers : number of conv layers (>= 2 recommended)
    use_gn : use lightweight GroupNorm on hidden layers
    """

    def __init__(self, in_ch: int = 3, out_ch: int = 2, width: int = 32,
                 ksize: int = 5, layers: int = 3, use_gn: bool = False, act: nn.Module = nn.GELU()):
        super().__init__()
        assert ksize % 2 == 1, "ksize must be odd for 'same' padding"
        pads = ksize // 2
        self.act = act
        self.use_residual = (in_ch == out_ch)

        blocks = []
        ch = in_ch
        for li in range(layers - 1):
            blocks.append(nn.Conv2d(ch, width, ksize, padding=pads))
            if use_gn:
                blocks.append(nn.GroupNorm(num_groups=min(8, width), num_channels=width))
            blocks.append(act)
            ch = width
        # last conv to out
        blocks.append(nn.Conv2d(ch, out_ch, ksize, padding=pads))
        self.net = nn.Sequential(*blocks)

        # init
        for m in self.modules():
            if isinstance(m, nn.Conv2d):
                nn.init.kaiming_uniform_(m.weight, a=math.sqrt(5))
                if m.bias is not None:
                    fan_in, _ = nn.init._calculate_fan_in_and_fan_out(m.weight)
                    bound = 1 / math.sqrt(fan_in)
                    nn.init.uniform_(m.bias, -bound, bound)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        y = self.net(x)
        if self.use_residual:
            # light residual “identity” if in_ch == out_ch
            y = y + x
        return y
