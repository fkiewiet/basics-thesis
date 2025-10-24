# src/models.py
from __future__ import annotations
import torch
import torch.nn as nn

class SpectralConv2d(nn.Module):
    def __init__(self, in_ch, out_ch, modes_h=12, modes_w=12):
        super().__init__()
        self.mh, self.mw = modes_h, modes_w
        self.weight = nn.Parameter(torch.randn(in_ch, out_ch, modes_h, modes_w, 2) * (1/(in_ch*out_ch)))

    def forward(self, x):
        B,C,H,W = x.shape
        x_ft = torch.fft.rfft2(x, norm="ortho")
        out_ft = torch.zeros(B, C, H, W//2+1, dtype=torch.cfloat, device=x.device)
        w = self.weight[...,0] + 1j*self.weight[...,1]
        mh, mw = self.mh, self.mw
        out_ft[:, :, :mh, :mw] = torch.einsum("bchw,cohw->bohw", x_ft[:, :, :mh, :mw], w)
        return torch.fft.irfft2(out_ft, s=(H,W), norm="ortho")

class SimpleFNO(nn.Module):
    def __init__(self, in_ch=3, out_ch=2, width=48, modes=(12,12), layers=4):
        super().__init__()
        self.fc0 = nn.Conv2d(in_ch, width, 1)
        self.sconvs = nn.ModuleList([SpectralConv2d(width, width, modes[0], modes[1]) for _ in range(layers)])
        self.ws     = nn.ModuleList([nn.Conv2d(width, width, 1) for _ in range(layers)])
        self.fc1 = nn.Conv2d(width, 64, 1)
        self.fc2 = nn.Conv2d(64, out_ch, 1)
        self.act = nn.GELU()
    def forward(self, x):
        x = self.fc0(x)
        for sc,w in zip(self.sconvs, self.ws):
            x = self.act(sc(x) + w(x))
        x = self.act(self.fc1(x))
        return self.fc2(x)

class LocalCNN(nn.Module):
    def __init__(self, in_ch=3, out_ch=2, width=32, ksize=5, layers=3):
        super().__init__()
        pads = ksize//2
        mods = []
        ch = in_ch
        for _ in range(layers-1):
            mods += [nn.Conv2d(ch, width, ksize, padding=pads), nn.GELU()]
            ch = width
        mods += [nn.Conv2d(ch, out_ch, ksize, padding=pads)]
        self.net = nn.Sequential(*mods)
    def forward(self, x): return self.net(x)
