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
    tol: float = 1e-6
    maxiter: Optional[int] = None


def gmres_solve(
    matrix: "sparse.spmatrix | LinearOperator",
    rhs: np.ndarray,
    *,
    options: Optional[GMRESOptions] = None,
    callback: Optional[Callable[[np.ndarray], None]] = None,
) -> SolverResult:
    """Wrap SciPy's GMRES; logs relative residuals and is compatible with tol / rtol APIs."""
    from scipy.sparse.linalg import gmres as scipy_gmres
    import inspect
    import numpy as np

    opts = options or GMRESOptions()
    residuals: list[float] = []

    # relative residual logging: ||r|| / ||b||
    b_norm = float(np.linalg.norm(rhs)) or 1.0  # guard against zero RHS

    def _callback(residual: np.ndarray) -> None:
        norm_rel = float(np.linalg.norm(residual) / b_norm)
        residuals.append(norm_rel)
        if callback is not None:
            callback(residual)

    # Build kwargs to match the installed SciPy signature
    sig = inspect.signature(scipy_gmres).parameters
    kwargs = dict(
        restart=opts.restart,
        maxiter=opts.maxiter,
        callback=_callback,
    )

    if opts.tol is not None:
        if "rtol" in sig:        # newer SciPy
            kwargs["rtol"] = float(opts.tol)   # relative tolerance
            kwargs["atol"] = 0.0               # no absolute floor
        else:                      # older SciPy
            kwargs["tol"] = float(opts.tol)    # legacy single tol

    solution, info = scipy_gmres(matrix, rhs, **kwargs)
    converged = (info == 0)
    return SolverResult(solution=solution, residuals=residuals, converged=converged, info=info)


