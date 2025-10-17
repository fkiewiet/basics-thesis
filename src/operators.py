"""Operators for the Helmholtz equation.

This module assembles sparse matrices for the Helmholtz operator

    (Δ + k^2) u = f

on structured grids described by `GridSpec`. It is designed for quick
experimentation in notebooks: you can easily swap boundary conditions,
wavenumber k, dtype, and (later) discretisation schemes.

Example
-------
>>> from src.grid import GridSpec
>>> from src.operators import helmholtz_operator, laplacian_operator
>>> grid = GridSpec(dims=2, shape=(60, 60), lengths=(1.0, 1.0))
>>> H = helmholtz_operator(grid, k=30.0, bc="periodic")  # (Δ + 30^2 I)
>>> L = laplacian_operator(grid, bc="dirichlet")
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from functools import reduce
from typing import Iterable, Protocol, Sequence

import numpy as np

# Optional dependency: scipy.sparse
try:  # pragma: no cover
    from scipy import sparse
except Exception as exc:  # pragma: no cover
    sparse = None  # type: ignore[assignment]
    _sparse_import_error = exc
else:  # pragma: no cover
    _sparse_import_error = None

from .config import GridSpec


# --------------------------------------------------------------------------- #
# Public API types
# --------------------------------------------------------------------------- #

class BC(str, Enum):
    """Supported boundary conditions."""
    DIRICHLET = "dirichlet"
    NEUMANN   = "neumann"
    PERIODIC  = "periodic"

    @classmethod
    def parse(cls, value: str | BC) -> "BC":
        if isinstance(value, BC):
            return value
        v = str(value).lower()
        try:
            return BC(v)
        except ValueError as e:
            opts = ", ".join(b.value for b in BC)
            raise ValueError(f"Unsupported boundary condition: {value!r}. Choose from {{{opts}}}.") from e


@dataclass(frozen=True)
class FDConfig:
    """Finite-difference configuration (second order central differences).

    Parameters
    ----------
    bc : {"dirichlet","neumann","periodic"}
        Boundary condition applied on all faces.
    dtype : str | np.dtype
        Matrix dtype ("float64" or "complex128" are typical). Helmholtz with k>0
        often uses complex dtype when you combine with complex-valued sources.
    """
    bc: BC = BC.DIRICHLET
    dtype: np.dtype | str = np.float64


class Discretisation(Protocol):
    """How to assemble the Helmholtz operator (Δ + k² I) on a grid."""
    name: str
    def assemble(self, grid: GridSpec, k: float) -> "sparse.spmatrix":
        ...


# --------------------------------------------------------------------------- #
# High-level functions (recommended for notebooks)
# --------------------------------------------------------------------------- #

def helmholtz_operator(
    grid: GridSpec,
    k: float,
    bc: str | BC = "dirichlet",
    *,
    dtype: np.dtype | str = np.complex128,
) -> "sparse.spmatrix":
    """Assemble the Helmholtz operator (Δ + k² I).

    Parameters
    ----------
    grid : GridSpec
        Structured grid descriptor.
    k : float
        Wavenumber (|k| in the conventional Δ + k² formulation).
    bc : {"dirichlet","neumann","periodic"}
        Boundary condition.
    dtype : numpy dtype or string
        Matrix dtype; use complex if you plan complex RHS/solutions.

    Returns
    -------
    scipy.sparse.csr_matrix
    """
    _ensure_sparse()
    cfg = FDConfig(bc=BC.parse(bc), dtype=np.dtype(dtype))
    L = _assemble_laplacian_fd(grid, cfg)
    return (L + (k ** 2) * sparse.eye(grid.size, dtype=cfg.dtype, format="csr")).tocsr()


def laplacian_operator(
    grid: GridSpec,
    bc: str | BC = "dirichlet",
    *,
    dtype: np.dtype | str = np.float64,
) -> "sparse.spmatrix":
    """Assemble the Laplacian with the given BC."""
    _ensure_sparse()
    cfg = FDConfig(bc=BC.parse(bc), dtype=np.dtype(dtype))
    return _assemble_laplacian_fd(grid, cfg)

def assemble_operator(
    discretisation: Discretisation,
    grid: GridSpec,
    wavenumber: float
) -> "sparse.spmatrix":
    """Legacy helper for notebooks expecting assemble_operator().

    Simply calls discretisation.assemble(grid, wavenumber).
    """
    _ensure_sparse()
    return discretisation.assemble(grid, wavenumber)



# --------------------------------------------------------------------------- #
# Class form (kept for convenience + backward compatibility)
# --------------------------------------------------------------------------- #

@dataclass(frozen=True)
class FiniteDifference:
    """Second-order central finite-difference Helmholtz.

    You can set defaults on the object and reuse it across calls.
    """
    config: FDConfig = FDConfig()
    # Optional default k to use if none is passed at call time
    wavenumber: float | None = None
    name: str = "finite_difference"

    def _resolve_k(self, k: float | None) -> float:
        return float(self.wavenumber if k is None and self.wavenumber is not None else (k or 0.0))

    def assemble(self, grid: GridSpec, k: float | None = None) -> "sparse.spmatrix":
        _ensure_sparse()
        L = _assemble_laplacian_fd(grid, self.config)
        kk = self._resolve_k(k)
        return (L + (kk ** 2) * sparse.eye(grid.size, dtype=self.config.dtype, format="csr")).tocsr()

    # --- Backwards-compatibility aliases used in earlier notebooks -----------
    def build_operator(self, grid: GridSpec, wavenumber: float | None = None) -> "sparse.spmatrix":
        """Alias for assemble(); k defaults to self.wavenumber or 0."""
        return self.assemble(grid, wavenumber)

    def build_matrix(self, grid: GridSpec, wavenumber: float | None = None) -> "sparse.spmatrix":
        """Alias for assemble(); kept for legacy code."""
        return self.assemble(grid, wavenumber)


# --------------------------------------------------------------------------- #
# Internals: FD Laplacian (2nd order) + BC handling
# --------------------------------------------------------------------------- #

def _assemble_laplacian_fd(grid: GridSpec, cfg: FDConfig) -> "sparse.spmatrix":
    """Kronecker-sum Laplacian in arbitrary dimensions with 2nd-order stencils."""
    shape = grid.shape
    spacing = grid.spacing
    dtype = cfg.dtype

    # Identity caches
    I = [sparse.eye(n, dtype=dtype, format="csr") for n in shape]
    comps = []
    for axis, (n, h) in enumerate(zip(shape, spacing)):
        D2 = _second_derivative_1d(n, h, dtype)
        kron = list(I)
        kron[axis] = D2
        comps.append(_kron_all(kron))
    L = reduce(lambda acc, term: acc + term, comps).tocsr()

    # Apply BCs
    if cfg.bc == BC.DIRICHLET:
        return _apply_dirichlet(L, shape, dtype)
    if cfg.bc == BC.NEUMANN:
        return _apply_neumann(L, shape, spacing).astype(dtype, copy=False).tocsr()
    if cfg.bc == BC.PERIODIC:
        return _apply_periodic(L, shape, spacing).astype(dtype, copy=False).tocsr()
    raise ValueError(f"Unsupported BC: {cfg.bc}")


def _second_derivative_1d(n: int, h: float, dtype) -> "sparse.spmatrix":
    main = np.full(n, -2.0 / h**2, dtype=dtype)
    off  = np.full(n - 1,  1.0 / h**2, dtype=dtype)
    return sparse.diags([off, main, off], offsets=[-1, 0, 1], dtype=dtype, format="csr")


def _kron_all(mats: Iterable["sparse.spmatrix"]) -> "sparse.spmatrix":
    mats = list(mats)
    return reduce(lambda a, b: sparse.kron(a, b, format="csr"), mats)


# ------------------------------ boundary ops --------------------------------

def _apply_dirichlet(M: "sparse.spmatrix", shape: Sequence[int], dtype) -> "sparse.spmatrix":
    """Homogeneous Dirichlet: overwrite boundary rows with identity (u=0 on ∂Ω)."""
    A = M.tolil()
    for idx in np.ndindex(tuple(shape)):
        if _on_boundary(idx, shape):
            r = _flatten_index(idx, shape)
            A.rows[r] = [r]
            A.data[r] = [np.array(1.0, dtype=dtype)]
    return A.tocsr()


def _apply_neumann(M: "sparse.spmatrix", shape: Sequence[int], spacing: Sequence[float]) -> "sparse.spmatrix":
    """
    Homogeneous Neumann: mirrored ghost-point adjustment at faces.
    Implementation: add +1/h^2 to the diagonal at nodes touching a boundary
    face orthogonal to that axis (keeps row-sum zero in 1D).
    """
    A = M.tolil()
    for axis, h in enumerate(spacing):
        h2 = h * h
        for idx in np.ndindex(tuple(shape)):
            if idx[axis] == 0 or idx[axis] == shape[axis] - 1:
                r = _flatten_index(idx, shape)
                A[r, r] = A[r, r] + 1.0 / h2
    return A.tocsr()


def _apply_periodic(M: "sparse.spmatrix", shape: Sequence[int], spacing: Sequence[float]) -> "sparse.spmatrix":
    """Periodic: connect opposite faces (wrap-around) along each axis."""
    A = M.tolil()
    for axis, (n, h) in enumerate(zip(shape, spacing)):
        h2 = h * h
        for idx in np.ndindex(tuple(shape)):
            r = _flatten_index(idx, shape)
            i = list(idx)

            # neighbor at -1 (wrap)
            i[axis] = (i[axis] - 1) % n
            A[r, _flatten_index(tuple(i), shape)] += 1.0 / h2

            # neighbor at +1 (wrap)
            i[axis] = (i[axis] + 2) % n  # move from i-1 to i+1
            A[r, _flatten_index(tuple(i), shape)] += 1.0 / h2
    return A.tocsr()


# ----------------------------- indexing helpers -----------------------------

def _on_boundary(idx: Sequence[int], shape: Sequence[int]) -> bool:
    return any(i == 0 or i == dim - 1 for i, dim in zip(idx, shape))


def _flatten_index(idx: Sequence[int], shape: Sequence[int]) -> int:
    """Row-major (C-order) flattening of a multi-index."""
    linear = 0
    stride = 1
    for i, dim in zip(reversed(idx), reversed(shape)):
        linear += i * stride
        stride *= dim
    return int(linear)


# --------------------------------------------------------------------------- #
# Housekeeping
# --------------------------------------------------------------------------- #

def _ensure_sparse() -> None:
    if sparse is None:  # pragma: no cover
        raise ImportError("scipy.sparse is required for operator assembly") from _sparse_import_error
