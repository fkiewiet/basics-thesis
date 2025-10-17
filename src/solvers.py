"""Iterative solvers and helper routines."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Optional

import numpy as np
from scipy.sparse.linalg import LinearOperator, gmres as scipy_gmres

try:  # pragma: no cover - optional dependency typing
    from scipy import sparse
except Exception:  # pragma: no cover
    sparse = None  # type: ignore


@dataclass
class SolverResult:
    """Container returned by solver wrappers."""

    solution: np.ndarray
    residuals: list[float]
    converged: bool
    info: int


@dataclass(frozen=True)
class GMRESOptions:
    """Basic knobs exposed for GMRES runs."""

    restart: Optional[int] = None
    tol: float = 1e-8
    maxiter: Optional[int] = None


def gmres_solve(
    matrix: "sparse.spmatrix | LinearOperator",
    rhs: np.ndarray,
    *,
    options: Optional[GMRESOptions] = None,
    callback: Optional[Callable[[np.ndarray], None]] = None,
) -> SolverResult:
    """Wrap SciPy's GMRES; works with both tol and rtol/atol SciPy versions."""
    from scipy.sparse.linalg import gmres as scipy_gmres
    import inspect

    opts = options or GMRESOptions()
    residuals: list[float] = []

    def _callback(residual: np.ndarray) -> None:
        residuals.append(float(np.linalg.norm(residual)))
        if callback is not None:
            callback(residual)

    # Build kwargs robustly
    sig = inspect.signature(scipy_gmres).parameters
    kwargs = dict(
        restart=opts.restart,
        maxiter=opts.maxiter,
        callback=_callback,
    )

    # Only pass tolerances if user provided one; else use SciPy defaults
    if opts.tol is not None:
        if "rtol" in sig:         # newer SciPy
            kwargs["rtol"] = float(opts.tol)
            kwargs["atol"] = 0.0
        elif "tol" in sig:        # older SciPy
            kwargs["tol"] = float(opts.tol)

    solution, info = scipy_gmres(matrix, rhs, **kwargs)
    return SolverResult(solution=solution, residuals=residuals, converged=(info == 0), info=info)

