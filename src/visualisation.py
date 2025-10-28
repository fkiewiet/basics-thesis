# src/visualisation.py
from __future__ import annotations

from typing import Iterable, Protocol, Sequence, TYPE_CHECKING, Tuple, Optional

import numpy as np

try:  # pragma: no cover - optional dependency
    import matplotlib.pyplot as plt
except Exception:  # pragma: no cover
    plt = None  # type: ignore


__all__ = [
    "plot_residuals",
    "plot_field",
    "compare_fields",
    "linecut",
    "save_snapshot",
    "plot_error_vs_k",
    "shade_pml",
]


# ----------------------------- Protocols / Types -----------------------------

class _SolverResultProtocol(Protocol):
    """Minimal surface required from solver outputs for plotting."""
    residuals: Sequence[float]


if TYPE_CHECKING:  # pragma: no cover - used only for static analysis
    from .solvers import SolverResult as SolverResult  # noqa: F401
else:
    SolverResult = _SolverResultProtocol


# ----------------------------- Small utilities -------------------------------

def _has_plt() -> bool:
    return plt is not None  # type: ignore[truthy-function]


def _imshow(ax, Z, title: str = "", cmap: str = "viridis"):
    """Safe imshow with colorbar; returns the image handle."""
    im = ax.imshow(Z, origin="lower", cmap=cmap)
    ax.set_title(title)
    ax.set_xticks([]); ax.set_yticks([])
    if _has_plt():
        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    return im


def _to_mag_phase(u_2ch: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    u_2ch: (2,H,W) -> (|u|, phase in [-pi,pi])
    """
    ur, ui = u_2ch[0], u_2ch[1]
    u = ur + 1j * ui
    mag = np.abs(u)
    phs = np.angle(u)
    return mag, phs


def _wrap_phase_diff(phi_pred: np.ndarray, phi_true: np.ndarray) -> np.ndarray:
    """Return wrapped phase difference in [-pi, pi]."""
    return np.angle(np.exp(1j * (phi_pred - phi_true)))


# ----------------------------- Core plots (kept API) -------------------------

def plot_residuals(
    result: SolverResult,
    *,
    ax: "plt.Axes | None" = None,
    tol: float | None = None
) -> "plt.Axes | None":
    """
    Semilogy residual history; optionally draw a horizontal tolerance line.
    """
    if not _has_plt():
        return None
    if ax is None:
        _, ax = plt.subplots()
    ax.semilogy(result.residuals)
    ax.set_xlabel("Iteration")
    ax.set_ylabel("Residual norm")
    ax.grid(True, which="both", alpha=0.25)
    if tol is not None:
        ax.axhline(y=tol, color="red", linestyle="--", linewidth=1.2,
                   label=f"tolerance = {tol:.0e}")
        ax.legend()
    return ax


def plot_field(
    field: np.ndarray,
    shape: Iterable[int],
    *,
    ax: "plt.Axes | None" = None,
    title: str = "",
    cmap: str = "viridis"
) -> "plt.Axes | None":
    """
    Plot 1D or 2D field. For 2D, shows real-part heatmap by default.
    """
    if not _has_plt():
        return None
    A = np.asarray(field).reshape(tuple(shape))
    if A.ndim == 1:
        if ax is None:
            _, ax = plt.subplots()
        ax.plot(A)
        ax.set_title(title or "Field")
        ax.set_xlabel("Index")
        ax.set_ylabel("Value")
        ax.grid(True, alpha=0.3)
    elif A.ndim == 2:
        if ax is None:
            _, ax = plt.subplots()
        im = ax.imshow(np.real(A), origin="lower", cmap=cmap)
        ax.set_title(title or "Field (Re)")
        ax.set_xticks([]); ax.set_yticks([])
        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    else:
        raise ValueError("plot_field supports 1D or 2D data only")
    return ax


# ----------------------------- PDE/ML helpers --------------------------------

def compare_fields(
    y_true_2ch: np.ndarray,
    y_pred_2ch: np.ndarray,
    *,
    title: str = "",
    show_phase: bool = True
) -> "plt.Figure | None":
    """
    Compare truth vs prediction for a complex field (Re/Im channels).
    Panels: |u|_true, |u|_pred, |Δ|u|| ; phase (optional);
            Re/Im truth vs pred; ΔRe/ΔIm heatmaps.
    """
    if not _has_plt():
        return None

    y_true_2ch = np.asarray(y_true_2ch)
    y_pred_2ch = np.asarray(y_pred_2ch)
    assert y_true_2ch.shape == y_pred_2ch.shape and y_true_2ch.shape[0] == 2, \
        "Inputs must be (2,H,W) arrays of identical shape."

    H, W = y_true_2ch.shape[1:]
    mag_t, phs_t = _to_mag_phase(y_true_2ch)
    mag_p, phs_p = _to_mag_phase(y_pred_2ch)
    dmag = np.abs(mag_p - mag_t)
    dRe = y_pred_2ch[0] - y_true_2ch[0]
    dIm = y_pred_2ch[1] - y_true_2ch[1]

    if show_phase:
        fig, axs = plt.subplots(3, 4, figsize=(14, 9))
        if title: fig.suptitle(title, y=0.98)
        _imshow(axs[0,0], mag_t, "True |u|")
        _imshow(axs[0,1], mag_p, "Pred |u|")
        _imshow(axs[0,2], dmag, "|Δ|u||")
        _imshow(axs[0,3], mag_p - mag_t, "Δ|u| (signed)", cmap="coolwarm")

        _imshow(axs[1,0], phs_t, "True phase", cmap="twilight")
        _imshow(axs[1,1], phs_p, "Pred phase", cmap="twilight")
        dph = _wrap_phase_diff(phs_p, phs_t)
        _imshow(axs[1,2], np.abs(dph), "|Δphase|", cmap="magma")
        _imshow(axs[1,3], dph, "Δphase", cmap="coolwarm")

        _imshow(axs[2,0], y_true_2ch[0], "Re(u) truth", cmap="coolwarm")
        _imshow(axs[2,1], y_pred_2ch[0], "Re(u) pred", cmap="coolwarm")
        _imshow(axs[2,2], y_true_2ch[1], "Im(u) truth", cmap="coolwarm")
        _imshow(axs[2,3], y_pred_2ch[1], "Im(u) pred", cmap="coolwarm")
        fig.tight_layout()
        return fig

    # No phase row
    fig, axs = plt.subplots(2, 4, figsize=(14, 6))
    if title: fig.suptitle(title, y=0.98)
    _imshow(axs[0,0], mag_t, "True |u|")
    _imshow(axs[0,1], mag_p, "Pred |u|")
    _imshow(axs[0,2], dmag, "|Δ|u||")
    _imshow(axs[0,3], mag_p - mag_t, "Δ|u| (signed)", cmap="coolwarm")

    _imshow(axs[1,0], dRe, "ΔRe", cmap="coolwarm")
    _imshow(axs[1,1], dIm, "ΔIm", cmap="coolwarm")
    _imshow(axs[1,2], y_true_2ch[0], "Re truth", cmap="coolwarm")
    _imshow(axs[1,3], y_true_2ch[1], "Im truth", cmap="coolwarm")
    fig.tight_layout()
    return fig


def linecut(
    y_true_2ch: np.ndarray,
    y_pred_2ch: np.ndarray,
    *,
    row: Optional[int] = None,
    col: Optional[int] = None,
    title: str = "Line cut"
) -> "plt.Figure | None":
    """
    Plot 1D slices through Re/Im along a row or column to see phase alignment.
    """
    if not _has_plt():
        return None
    y_true_2ch = np.asarray(y_true_2ch)
    y_pred_2ch = np.asarray(y_pred_2ch)
    H, W = y_true_2ch.shape[1:]

    if row is None and col is None:
        row = H // 2

    fig = plt.figure(figsize=(8, 4))
    fig.suptitle(title)

    if row is not None:
        x = np.arange(W)
        t_Re, p_Re = y_true_2ch[0, row, :], y_pred_2ch[0, row, :]
        t_Im, p_Im = y_true_2ch[1, row, :], y_pred_2ch[1, row, :]
        ax1 = fig.add_subplot(2, 1, 1)
        ax2 = fig.add_subplot(2, 1, 2, sharex=ax1)
        ax1.plot(x, t_Re, label="True Re"); ax1.plot(x, p_Re, "--", label="Pred Re")
        ax2.plot(x, t_Im, label="True Im"); ax2.plot(x, p_Im, "--", label="Pred Im")
        ax1.legend(); ax2.legend()
        ax1.set_ylabel("Re"); ax2.set_ylabel("Im"); ax2.set_xlabel(f"Column index (row={row})")
        ax1.grid(True, alpha=0.3); ax2.grid(True, alpha=0.3)
        return fig

    if col is not None:
        y = np.arange(H)
        t_Re, p_Re = y_true_2ch[0, :, col], y_pred_2ch[0, :, col]
        t_Im, p_Im = y_true_2ch[1, :, col], y_pred_2ch[1, :, col]
        ax1 = fig.add_subplot(2, 1, 1)
        ax2 = fig.add_subplot(2, 1, 2, sharex=ax1)
        ax1.plot(y, t_Re, label="True Re"); ax1.plot(y, p_Re, "--", label="Pred Re")
        ax2.plot(y, t_Im, label="True Im"); ax2.plot(y, p_Im, "--", label="Pred Im")
        ax1.legend(); ax2.legend()
        ax1.set_ylabel("Re"); ax2.set_ylabel("Im"); ax2.set_xlabel(f"Row index (col={col})")
        ax1.grid(True, alpha=0.3); ax2.grid(True, alpha=0.3)
        return fig

    return fig


def save_snapshot(
    y_true_2ch: np.ndarray,
    y_pred_2ch: np.ndarray,
    out_path: "str | None" = None,
    *,
    title: str = ""
) -> "plt.Figure | None":
    """
    Save a compact 2x3 figure (True |u|, Pred |u|, |Δ|u||, True Re, Pred Re, ΔRe).
    If out_path is None, just shows the figure.
    """
    if not _has_plt():
        return None
    y_true_2ch = np.asarray(y_true_2ch)
    y_pred_2ch = np.asarray(y_pred_2ch)

    def to_mag(u2): return np.abs(u2[0] + 1j * u2[1])
    mag_t, mag_p = to_mag(y_true_2ch), to_mag(y_pred_2ch)

    fig, axs = plt.subplots(2, 3, figsize=(10, 6))
    if title: fig.suptitle(title)
    for ax in axs.ravel():
        ax.set_xticks([]); ax.set_yticks([])

    _imshow(axs[0,0], mag_t, "True |u|")
    _imshow(axs[0,1], mag_p, "Pred |u|")
    _imshow(axs[0,2], np.abs(mag_p - mag_t), "|Δ|u||")
    _imshow(axs[1,0], y_true_2ch[0], "True Re", cmap="coolwarm")
    _imshow(axs[1,1], y_pred_2ch[0], "Pred Re", cmap="coolwarm")
    _imshow(axs[1,2], y_pred_2ch[0] - y_true_2ch[0], "ΔRe", cmap="coolwarm")

    fig.tight_layout()
    if out_path is not None:
        fig.savefig(out_path, dpi=200)
    else:
        plt.show()
    return fig


def plot_error_vs_k(
    ks: np.ndarray,
    rel_errors: np.ndarray,
    *,
    title: str = "Error vs frequency",
    ylabel: str = "relative L2 error"
) -> "plt.Axes | None":
    """
    Scatter plot of error vs. frequency k. Accepts 1D arrays of equal length.
    """
    if not _has_plt():
        return None
    ks = np.asarray(ks).ravel()
    rel = np.asarray(rel_errors).ravel()
    assert ks.shape == rel.shape, "ks and rel_errors must have the same length."
    fig, ax = plt.subplots(figsize=(6, 3.5))
    ax.scatter(ks, rel, s=14)
    ax.set_xlabel("k")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    return ax


def shade_pml(ax, H: int, W: int, pml_thickness: int) -> None:
    """
    Draw dashed rectangles indicating PML borders on a 2D image plot.
    Call this after imshow; (H,W) are the array sizes.
    """
    if not _has_plt():
        return
    t = int(pml_thickness)
    if t <= 0:
        return
    import matplotlib.patches as patches
    # top
    ax.add_patch(patches.Rectangle((0, 0), W, t, fill=False, ec="k", ls="--", lw=0.8))
    # bottom
    ax.add_patch(patches.Rectangle((0, H - t), W, t, fill=False, ec="k", ls="--", lw=0.8))
    # left
    ax.add_patch(patches.Rectangle((0, 0), t, H, fill=False, ec="k", ls="--", lw=0.8))
    # right
    ax.add_patch(patches.Rectangle((W - t, 0), t, H, fill=False, ec="k", ls="--", lw=0.8))


# ----------------------------- CLI guidance ----------------------------------

if __name__ == "__main__":  # pragma: no cover
    import textwrap
    message = textwrap.dedent(
        """
        Plotting helpers for your Helmholtz ML notebooks.

        Typical usage:
            from src.visualisation import compare_fields, linecut, plot_residuals

        - compare_fields(y_true_2ch, y_pred_2ch, title="(f,k)->u")
        - linecut(y_true_2ch, y_pred_2ch, row=H//2)
        - plot_residuals(result, tol=1e-8)

        All functions degrade gracefully if matplotlib is not installed.
        """
    ).strip()
    print(message)
