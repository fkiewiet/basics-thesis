# src/transfer.py
#from __future__ import annotations
import torch
import torch.nn as nn

class LearnedTransfer(nn.Module):
    """Local operator T; can be initialized with a LocalCNN."""
    def __init__(self, core: nn.Module):
        super().__init__()
        self.core = core
    def forward(self, u): return self.core(u)

def intertwining_penalty(T: nn.Module, u, apply_L_lo, apply_L_hi):
    """|| L_hi(T u) - T(L_lo u) ||_2"""
    Tu = T(u)
    return torch.linalg.vector_norm(apply_L_hi(Tu) - T(apply_L_lo(u)))
