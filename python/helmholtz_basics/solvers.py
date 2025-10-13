"""Iterative solvers and helper routines."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Optional

import numpy as np
from scipy.sparse.linalg import LinearOperator, gmres as scipy_gmres


@dataclass
class SolverResult:
    solution: np.ndarray
    residuals: list[float]
    converged: bool
    info: int


def gmres_solve(
    matrix: "sparse.spmatrix | LinearOperator",
    rhs: np.ndarray,
    restart: Optional[int] = None,
    tol: float = 1e-8,
    maxiter: Optional[int] = None,
    callback: Optional[Callable[[np.ndarray], None]] = None,
) -> SolverResult:
    """Wrap SciPy's GMRES with a friendlier result object."""

    residuals: list[float] = []

    def _callback(residual: np.ndarray) -> None:
        norm = float(np.linalg.norm(residual))
        residuals.append(norm)
        if callback is not None:
            callback(residual)

    solution, info = scipy_gmres(matrix, rhs, restart=restart, tol=tol, maxiter=maxiter, callback=_callback)
    converged = info == 0
    return SolverResult(solution=solution, residuals=residuals, converged=converged, info=info)


try:  # pragma: no cover - optional dependency typing
    from scipy import sparse
except Exception:  # pragma: no cover
    sparse = None  # type: ignore
