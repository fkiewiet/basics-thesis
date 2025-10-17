from __future__ import annotations

from typing import Iterable, Protocol, Sequence, TYPE_CHECKING

import numpy as np

try:  # pragma: no cover - optional dependency
    import matplotlib.pyplot as plt
except Exception:  # pragma: no cover
    plt = None  # type: ignore

class _SolverResultProtocol(Protocol):
    """Minimal surface required from solver outputs for plotting."""

    residuals: Sequence[float]


if TYPE_CHECKING:  # pragma: no cover - used only for static analysis
    from .solvers import SolverResult as SolverResult  # noqa: F401
else:
    SolverResult = _SolverResultProtocol


def plot_residuals(result: SolverResult, *, ax: "plt.Axes | None" = None, tol: float | None = None) -> "plt.Axes | None":
    if plt is None:
        return None
    if ax is None:
        _, ax = plt.subplots()
    ax.semilogy(result.residuals) #, marker="o")
    ax.set_xlabel("Iteration")
    ax.set_ylabel("Residual norm")
    if tol is not None:
        ax.axhline(y=tol, color="red", linestyle="--", linewidth=1.2,
                   label=f"tolerance = {tol:.0e}")
        ax.legend()
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


if __name__ == "__main__":  # pragma: no cover - guidance for direct execution
    import textwrap

    message = textwrap.dedent(
        """
        These plotting helpers are meant to be imported from Python (for example
        `from helmholtz_basics.visualisation import plot_residuals`).

        If you would like to take them for a quick spin from the command line,
        `cd` into the `python/` directory of the repository (or add it to your
        `PYTHONPATH`) and run `python -m helmholtz_basics.visualisation`.  You
        can also open a Python or notebook session and import the functions
        there.
        """
    ).strip()
    print(message)
