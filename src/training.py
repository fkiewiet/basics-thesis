# src/training.py
from __future__ import annotations
from typing import Dict, List, Tuple
import numpy as np
import torch
import torch.nn as nn
from torch.utils.data import Dataset, DataLoader, random_split

try:
    import matplotlib.pyplot as plt
except Exception:
    plt = None  # pragma: no cover

class XYDictDataset(Dataset):
    def __init__(self, data: Dict[str, np.ndarray]):
        self.x = torch.from_numpy(data["x"]).float()
        self.y = torch.from_numpy(data["y"]).float()
    def __len__(self): return self.x.shape[0]
    def __getitem__(self, i): return self.x[i], self.y[i]

def train_model(model: torch.nn.Module, data: Dict[str, np.ndarray],
                epochs=20, batch_size=8, lr=1e-3, device=None):
    device = device or torch.device("cuda" if torch.cuda.is_available() else "cpu")
    ds = XYDictDataset(data)
    n_val = max(1, int(0.2*len(ds)))
    n_tr  = len(ds) - n_val
    tr, va = random_split(ds, [n_tr, n_val])
    tl = DataLoader(tr, batch_size=batch_size, shuffle=True)
    vl = DataLoader(va, batch_size=batch_size)
    model = model.to(device)
    opt = torch.optim.Adam(model.parameters(), lr=lr)
    crit = nn.MSELoss()
    hist = {"train": [], "val": []}
    for ep in range(1, epochs+1):
        model.train(); tot=0
        for x,y in tl:
            x,y = x.to(device), y.to(device)
            opt.zero_grad(); p = model(x); loss = crit(p,y)
            loss.backward(); opt.step()
            tot += loss.item()*x.size(0)
        tr_loss = tot/len(tr)
        model.eval(); tot=0
        with torch.no_grad():
            for x,y in vl:
                x,y = x.to(device), y.to(device)
                p = model(x); tot += crit(p,y).item()*x.size(0)
        va_loss = tot/len(va)
        hist["train"].append(tr_loss); hist["val"].append(va_loss)
        print(f"Epoch {ep:02d} | train {tr_loss:.4e} | val {va_loss:.4e}")
    return model, hist

def show_example(model: torch.nn.Module, data: Dict[str, np.ndarray], index=0, title=""):
    if plt is None:
        return
    model.eval()
    x = torch.from_numpy(data["x"][index:index+1]).float()
    with torch.no_grad():
        p = model(x.to(next(model.parameters()).device)).cpu().numpy()[0]
    y = data["y"][index]
    fig, ax = plt.subplots(2, 3, figsize=(11,6))
    # Inputs: ch0, ch1, ch2
    for j in range(3):
        im = ax[0, j].imshow(x.numpy()[0, j]); ax[0, j].set_title(f"input ch{j}"); ax[0, j].axis("off"); fig.colorbar(im, ax=ax[0, j], fraction=0.046, pad=0.04)
    # Truth vs pred
    im = ax[1,0].imshow(y[0]); ax[1,0].set_title("Re(u) truth"); ax[1,0].axis("off"); fig.colorbar(im, ax=ax[1,0], fraction=0.046, pad=0.04)
    im = ax[1,1].imshow(y[1]); ax[1,1].set_title("Im(u) truth"); ax[1,1].axis("off"); fig.colorbar(im, ax=ax[1,1], fraction=0.046, pad=0.04)
    magp = np.linalg.norm(p, axis=0)
    im = ax[1,2].imshow(magp); ax[1,2].set_title("|u| pred"); ax[1,2].axis("off"); fig.colorbar(im, ax=ax[1,2], fraction=0.046, pad=0.04)
    fig.suptitle(title); fig.tight_layout(); plt.show()
