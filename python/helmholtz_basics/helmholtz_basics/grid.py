"""Utilities for generating structured grids."""

from __future__ import annotations

from typing import Tuple

import numpy as np

from .config import GridSpec


def build_grid(spec: GridSpec) -> Tuple[np.ndarray, ...]:
    """Return coordinate arrays for each dimension."""

    axes = []
    for length, n in zip(spec.lengths, spec.shape):
        axes.append(np.linspace(0.0, length, n, dtype=float))
    return tuple(np.meshgrid(*axes, indexing="ij", sparse=False))