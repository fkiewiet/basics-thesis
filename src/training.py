# src/training.py
from __future__ import annotations
from typing import Dict, List, Tuple, Optional
import math
import numpy as np
import torch
import torch.nn as nn
from torch.utils.data import Dataset, DataLoader, random_split

try:
    import matplotlib.pyplot as plt
except Exception:  # pragma: no cover
    plt = None  # pragma: no cover


# ----------------------------- Dataset wrapper -----------------------------

class XYDictDataset(Dataset):
    """
    Wrap {"x": (N,C,H,W), "y": (N,2,H,W)} into a torch Dataset.
    Keeps baseline behavior: tensors are created eagerly (CPU) as float32.
    """
    def __init__(self, data: Dict[str, np.ndarray]):
        self.x = torch.from_numpy(np.asarray(data["x"], dtype=np.float32))
        self.y = torch.from_numpy(np.asarray(data["y"], dtype=np.float32))

    def __len__(self):  # type: ignore[override]
        return self.x.shape[0]

    def __getitem__(self, i):  # type: ignore[override]
        return self.x[i], self.y[i]


# ----------------------------- Loader helper -------------------------------

def make_loaders(
    data: Dict[str, np.ndarray],
    batch_size: int = 8,
    val_frac: float = 0.2,
    seed: int = 0,
    num_workers: int = 0,
    pin_memory: bool = False,
) -> Tuple[DataLoader, DataLoader]:
    """
    Deterministically split into train/val and return DataLoaders.
    """
    ds = XYDictDataset(data)
    n_total = len(ds)
    n_val = max(1, int(math.floor(val_frac * n_total))) if 0 < val_frac < 1 else 1
    n_tr = max(1, n_total - n_val)
    g = torch.Generator().manual_seed(seed)
    tr, va = random_split(ds, [n_tr, n_val], generator=g)

    tl = DataLoader(tr, batch_size=batch_size, shuffle=True,
                    num_workers=num_workers, pin_memory=pin_memory)
    vl = DataLoader(va, batch_size=batch_size, shuffle=False,
                    num_workers=num_workers, pin_memory=pin_memory)
    return tl, vl


# ----------------------------- Train / Eval --------------------------------

def train_model(
    model: torch.nn.Module,
    data: Dict[str, np.ndarray],
    *,
    epochs: int = 20,
    batch_size: int = 8,
    lr: float = 1e-3,
    device: Optional[torch.device | str] = None,
    # extras (optional, all have safe defaults)
    val_frac: float = 0.2,
    seed: int = 0,
    weight_decay: float = 0.0,
    clip_grad_norm: Optional[float] = None,
    use_amp: bool = False,
    scheduler: Optional[str] = None,      # {"cosine","step",None}
    scheduler_kw: Optional[dict] = None,  # e.g. {"step_size":10,"gamma":0.5}
    early_stopping_patience: Optional[int] = None,
    restore_best: bool = True,
    num_workers: int = 0,
    pin_memory: bool = False,
    verbose: bool = True,
):
    """
    Extended trainer that preserves your baseline API and return shape.
    Returns (model, hist) where hist is a dict: {"train":[...], "val":[...]} (compatible).
    """
    device = device or torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = model.to(device)

    tl, vl = make_loaders(
        data, batch_size=batch_size, val_frac=val_frac, seed=seed,
        num_workers=num_workers, pin_memory=pin_memory
    )
    crit = nn.MSELoss()
    opt = torch.optim.Adam(model.parameters(), lr=lr, weight_decay=weight_decay)

    # optional schedulers
    scheduler_kw = scheduler_kw or {}
    if scheduler == "cosine":
        # CosineAnnealingLR over full epoch count
        sch = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=max(1, epochs))
    elif scheduler == "step":
        sch = torch.optim.lr_scheduler.StepLR(
            opt,
            step_size=int(scheduler_kw.get("step_size", 10)),
            gamma=float(scheduler_kw.get("gamma", 0.5)),
        )
    else:
        sch = None

    scaler = torch.cuda.amp.GradScaler(enabled=use_amp)

    hist: Dict[str, List[float]] = {"train": [], "val": []}
    best_val: float = float("inf")
    best_state: Optional[dict] = None
    patience_left = early_stopping_patience

    for ep in range(1, epochs + 1):
        # ---------------- train ----------------
        model.train()
        tot, n = 0.0, 0
        for xb, yb in tl:
            xb, yb = xb.to(device, non_blocking=pin_memory), yb.to(device, non_blocking=pin_memory)
            opt.zero_grad(set_to_none=True)
            with torch.cuda.amp.autocast(enabled=use_amp):
                pred = model(xb)
                loss = crit(pred, yb)
            scaler.scale(loss).backward()
            if clip_grad_norm is not None:
                scaler.unscale_(opt)
                nn.utils.clip_grad_norm_(model.parameters(), max_norm=clip_grad_norm)
            scaler.step(opt)
            scaler.update()
            tot += float(loss.detach()) * xb.size(0)
            n += xb.size(0)
        tr_loss = tot / max(1, n)

        # ---------------- val ------------------
        model.eval()
        tot, n = 0.0, 0
        with torch.no_grad():
            for xb, yb in vl:
                xb, yb = xb.to(device, non_blocking=pin_memory), yb.to(device, non_blocking=pin_memory)
                pred = model(xb)
                loss = crit(pred, yb)
                tot += float(loss) * xb.size(0)
                n += xb.size(0)
        va_loss = tot / max(1, n)

        hist["train"].append(tr_loss)
        hist["val"].append(va_loss)

        if verbose:
            lr_now = opt.param_groups[0]["lr"]
            print(f"Epoch {ep:03d}/{epochs} | train {tr_loss:.4e} | val {va_loss:.4e} | lr {lr_now:.2e}")

        # scheduler step (epoch-wise)
        if sch is not None:
            sch.step()

        # early stopping tracking
        if va_loss < best_val - 1e-12:
            best_val = va_loss
            if restore_best:
                best_state = {k: v.detach().cpu().clone() for k, v in model.state_dict().items()}
        else:
            if early_stopping_patience is not None:
                patience_left = early_stopping_patience if patience_left is None else patience_left - 1
                if patience_left <= 0:
                    if verbose:
                        print(f"Early stopping at epoch {ep}. Best val {best_val:.4e}.")
                    break

    # restore best weights if requested
    if restore_best and best_state is not None:
        model.load_state_dict(best_state)

    return model, hist


# ----------------------------- Extra evaluation -----------------------------

def evaluate_mse(model: torch.nn.Module, data: Dict[str, np.ndarray]) -> float:
    """
    Full-dataset MSE. (Convenience for quick checks.)
    """
    device = next(model.parameters()).device
    model.eval()
    x = torch.from_numpy(np.asarray(data["x"], dtype=np.float32)).to(device)
    y = torch.from_numpy(np.asarray(data["y"], dtype=np.float32)).to(device)
    with torch.no_grad():
        p = model(x)
        mse = torch.mean((p - y) ** 2).item()
    return float(mse)


def eval_relative_metrics(model: torch.nn.Module, data: Dict[str, np.ndarray]) -> Dict[str, float]:
    """
    Relative L2 stats and simple mag/phase RMSEs for interpretability.
    """
    device = next(model.parameters()).device
    model.eval()
    x = torch.from_numpy(np.asarray(data["x"], dtype=np.float32)).to(device)
    y = torch.from_numpy(np.asarray(data["y"], dtype=np.float32)).to(device)
    with torch.no_grad():
        p = model(x).cpu()
    y = y.cpu()

    num = torch.linalg.norm((p - y).flatten(1), dim=1)
    den = torch.linalg.norm(y.flatten(1), dim=1) + 1e-12
    rel = (num / den).numpy()

    p_np, y_np = p.numpy(), y.numpy()
    u_p = p_np[:, 0] + 1j * p_np[:, 1]
    u_y = y_np[:, 0] + 1j * y_np[:, 1]
    dmag_rmse = float(np.sqrt(np.mean((np.abs(u_p) - np.abs(u_y)) ** 2)))
    dphi = np.angle(np.exp(1j * (np.angle(u_p) - np.angle(u_y))))
    phase_rmse = float(np.sqrt(np.mean(dphi ** 2)))

    return {
        "mse": float(torch.mean((p - y) ** 2).item()),
        "rel_L2_mean": float(rel.mean()),
        "rel_L2_median": float(np.median(rel)),
        "rel_L2_p90": float(np.percentile(rel, 90)),
        "mag_rmse": dmag_rmse,
        "phase_rmse_rad": phase_rmse,
    }


# ----------------------------- Visualization -------------------------------

def _imshow(ax, Z, title="", cmap="viridis"):
    im = ax.imshow(Z, origin="lower", cmap=cmap)
    ax.set_title(title)
    ax.set_xticks([]); ax.set_yticks([])
    if plt is not None:
        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    return im


def show_example(model: torch.nn.Module, data: Dict[str, np.ndarray], index: int = 0, title: str = ""):
    if plt is None:
        return

    # --- robust device detection for models without parameters ---
    try:
        device = next(model.parameters()).device
    except StopIteration:
        device = torch.device("cpu")

    model.eval()

    x_np = np.asarray(data["x"][index], dtype=np.float32)   # (C,H,W)
    y_np = np.asarray(data["y"][index], dtype=np.float32)   # (2,H,W)

    pred_ok = True
    try:
        x = torch.from_numpy(x_np[None]).to(device)
        with torch.no_grad():
            p = model(x).cpu().numpy()[0]  # expect (2,H,W)
        if p.shape[0] != 2:
            pred_ok = False
    except Exception:
        pred_ok = False

    C, H, W = x_np.shape
    fig, ax = plt.subplots(2, 3, figsize=(11, 6))
    if title:
        fig.suptitle(title)

    # Inputs (first 3 channels)
    for j in range(min(3, C)):
        im = ax[0, j].imshow(x_np[j], origin="lower")
        ax[0, j].set_title(f"input ch{j}"); ax[0, j].axis("off")
        plt.colorbar(im, ax=ax[0, j], fraction=0.046, pad=0.04)

    # Truth
    im = ax[1,0].imshow(y_np[0], origin="lower", cmap="coolwarm"); ax[1,0].set_title("Re(u) truth"); ax[1,0].axis("off"); plt.colorbar(im, ax=ax[1,0], fraction=0.046, pad=0.04)
    im = ax[1,1].imshow(y_np[1], origin="lower", cmap="coolwarm"); ax[1,1].set_title("Im(u) truth"); ax[1,1].axis("off"); plt.colorbar(im, ax=ax[1,1], fraction=0.046, pad=0.04)

    # Pred magnitude (if available), else leave note
    if pred_ok:
        magp = np.abs(p[0] + 1j*p[1])
        im = ax[1,2].imshow(magp, origin="lower"); ax[1,2].set_title("|u| pred"); ax[1,2].axis("off"); plt.colorbar(im, ax=ax[1,2], fraction=0.046, pad=0.04)
    else:
        ax[1,2].axis("off")
        ax[1,2].text(0.5, 0.5, "no prediction\n(model incompatible)", ha="center", va="center", fontsize=10)

    fig.tight_layout(); plt.show()
