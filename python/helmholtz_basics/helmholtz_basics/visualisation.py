"""Optional Matplotlib-based plotting helpers."""

from __future__ import annotations

from typing import Iterable

import numpy as np

try:  # pragma: no cover - optional dependency
    import matplotlib.pyplot as plt
except Exception:  # pragma: no cover
    plt = None  # type: ignore

from .solvers import SolverResult


def plot_residuals(result: SolverResult, *, ax: "plt.Axes | None" = None) -> "plt.Axes | None":
    if plt is None:
        return None
    if ax is None:
        _, ax = plt.subplots()
    ax.semilogy(result.residuals, marker="o")
    ax.set_xlabel("Iteration")
    ax.set_ylabel("Residual norm")
    ax.set_title("GMRES convergence")
    return ax


def plot_field(field: np.ndarray, shape: Iterable[int], *, ax: "plt.Axes | None" = None) -> "plt.Axes | None":
    if plt is None:
        return None
    field = field.reshape(tuple(shape))
    if field.ndim == 1:
        if ax is None:
            _, ax = plt.subplots()
        ax.plot(field)
    elif field.ndim == 2:
        if ax is None:
            _, ax = plt.subplots()
        im = ax.imshow(np.real(field), origin="lower")
        plt.colorbar(im, ax=ax)
    else:
        raise ValueError("plot_field supports 1D or 2D data only")
    return ax
