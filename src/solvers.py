"""Iterative solvers and helper routines."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Optional, Any, Dict
from types import SimpleNamespace


import numpy as np
from scipy.sparse.linalg import LinearOperator, gmres as scipy_gmres
from scipy.sparse import csc_matrix, isspmatrix_csr, isspmatrix_csc
from scipy.sparse.linalg import spsolve, splu

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




# -- Optional settings for the direct solve
@dataclass
class DirectOptions:
    """Knobs for direct sparse linear solves."""
    method: str = "spsolve"        # "spsolve" | "splu" | "factorized"
    permc_spec: str = "COLAMD"     # ordering for splu
    diag_pivot_thresh: float = 1.0 # pivot threshold for splu
    cache: bool = False            # if you plan to reuse LU outside this call

def _to_csc(A):
    """Ensure sparse matrix is CSC-format for SuperLU."""
    if isspmatrix_csc(A):
        return A
    if isspmatrix_csr(A):
        return A.tocsc()
    return csc_matrix(A)

def direct_solve(A, b, options: Optional[DirectOptions] = None):
    """
    Direct linear solve for A x = b (complex supported).
    Returns an object with .solution, .converged, .info.
    """
    from types import SimpleNamespace
    from scipy.sparse.linalg import spsolve, splu, factorized

    options = options or DirectOptions()
    A = _to_csc(A)
    b = np.asarray(b)

    method = options.method.lower()
    info: Dict[str, Any] = {"method": method, "nnz": A.nnz, "shape": A.shape}

    if method == "splu":
        # explicit LU decomposition
        lu = splu(A, permc_spec=options.permc_spec, diag_pivot_thresh=options.diag_pivot_thresh)
        x = lu.solve(b)
        info.update({
            "permc_spec": options.permc_spec,
            "diag_pivot_thresh": options.diag_pivot_thresh,
        })

    elif method == "factorized":
        # cacheable solver (useful if same A is solved many times)
        solve = factorized(A)
        x = solve(b)

    else:
        # safest one-shot sparse solver (SciPy >=1.10)
        x = spsolve(A, b)  # no extra kwargs

    # mimic GMRES-style result for API consistency
    return SimpleNamespace(
        solution=x.astype(np.complex128, copy=False),
        converged=True,
        info=info,
        residuals=None,
    )




# === PML-backed Helmholtz helpers ===========================================
from .operators import helmholtz_operator  # local package import

def direct_solve_helmholtz(
    *,
    shape: tuple[int, ...],
    lengths: tuple[float, ...],
    k: float,
    rhs: np.ndarray,
    dtype: np.dtype = np.complex128,
    pml=None,
    options: Optional[DirectOptions] = None,
):
    """
    Assemble the Helmholtz operator A = (-Δ - k^2) with optional PML and solve A u = rhs.
    Returns a SimpleNamespace like direct_solve(...), with .solution, .converged, .info.

    Parameters
    ----------
    shape, lengths : grid definition (same as elsewhere in src)
    k              : wavenumber
    rhs            : right-hand side (1D ndarray, flattened in the same ordering as assembly)
    dtype          : complex dtype
    pml            : PMLConfig or None (Dirichlet if None)
    options        : DirectOptions (method, pivots, etc.)
    """
    A = helmholtz_operator(shape=shape, lengths=lengths, k=k, dtype=dtype, pml=pml)
    return direct_solve(A, rhs, options=options)


def solve_helmholtz_field(
    *,
    shape: tuple[int, ...],
    lengths: tuple[float, ...],
    k: float,
    rhs: np.ndarray,
    dtype: np.dtype = np.complex128,
    pml=None,
) -> np.ndarray:
    """
    Convenience wrapper that just returns the solution array (np.ndarray).
    Useful when dataset builders expect a plain vector rather than a result object.
    """
    res = direct_solve_helmholtz(shape=shape, lengths=lengths, k=k, rhs=rhs, dtype=dtype, pml=pml)
    return res.solution


def direct_solve_auto(
    *,
    shape: tuple[int, ...],
    lengths: tuple[float, ...],
    k: float,
    rhs: np.ndarray,
    use_pml: bool = False,
    pml=None,
    dtype: np.dtype = np.complex128,
    options: Optional[DirectOptions] = None,
):
    """
    One entry point you can switch in notebooks:
      - use_pml=False  -> Dirichlet
      - use_pml=True   -> PML (requires pml=PMLConfig)
    """
    pml_cfg = pml if use_pml else None
    return direct_solve_helmholtz(
        shape=shape, lengths=lengths, k=k, rhs=rhs, dtype=dtype, pml=pml_cfg, options=options
    )
# ============================================================================
