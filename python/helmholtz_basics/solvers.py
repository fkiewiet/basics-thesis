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
    """Wrap SciPy's GMRES with a friendlier result object."""

    opts = options or GMRESOptions()
    residuals: list[float] = []

    def _callback(residual: np.ndarray) -> None:
        norm = float(np.linalg.norm(residual))
        residuals.append(norm)
        if callback is not None:
            callback(residual)

    solution, info = scipy_gmres(
        matrix,
        rhs,
        restart=opts.restart,
        tol=opts.tol,
        maxiter=opts.maxiter,
        callback=_callback,
    )
    converged = info == 0
    return SolverResult(solution=solution, residuals=residuals, converged=converged, info=info)
