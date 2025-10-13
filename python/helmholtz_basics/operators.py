"""Discretisation strategies for the Helmholtz operator."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol, Tuple

import numpy as np

try:  # pragma: no cover - optional dependency
    from scipy import sparse
except Exception as exc:  # pragma: no cover
    sparse = None  # type: ignore[assignment]
    _sparse_import_error = exc
else:
    _sparse_import_error = None

from .config import GridSpec


class Discretisation(Protocol):
    """Protocol describing how to assemble the Helmholtz operator."""

    name: str

    def assemble(self, grid: GridSpec, wavenumber: float) -> sparse.spmatrix:
        """Return the discrete Helmholtz operator."""


@dataclass(frozen=True)
class FiniteDifference:
    """Simple second-order finite difference discretisation."""

    order: int = 2
    stencil: str = "standard"
    name: str = "finite_difference"

    def assemble(self, grid: GridSpec, wavenumber: float) -> sparse.spmatrix:
        _ensure_sparse()
        if grid.dims == 1:
            return _assemble_1d(grid, wavenumber)
        if grid.dims == 2:
            return _assemble_2d(grid, wavenumber)
        if grid.dims == 3:
            return _assemble_3d(grid, wavenumber)
        raise ValueError("Unsupported dimension")


def assemble_operator(discretisation: Discretisation, grid: GridSpec, wavenumber: float) -> sparse.spmatrix:
    """Public helper mirroring the Julia API."""

    _ensure_sparse()
    return discretisation.assemble(grid, wavenumber)


def _assemble_1d(grid: GridSpec, wavenumber: float) -> sparse.spmatrix:
    n = grid.shape[0]
    h = grid.spacing[0]
    main = (2.0 / h**2 + wavenumber**2) * np.ones(n)
    off = (-1.0 / h**2) * np.ones(n - 1)
    A = sparse.diags([off, main, off], offsets=[-1, 0, 1], format="csr")
    A = A.tolil()
    A[0, :] = 0.0
    A[0, 0] = 1.0
    A[-1, :] = 0.0
    A[-1, -1] = 1.0
    return A.tocsr()


def _assemble_2d(grid: GridSpec, wavenumber: float) -> sparse.spmatrix:
    nx, ny = grid.shape
    hx, hy = grid.spacing
    main = 2.0 / hx**2 + 2.0 / hy**2 + wavenumber**2
    off_x = -1.0 / hx**2
    off_y = -1.0 / hy**2

    Ix = sparse.eye(nx)
    Iy = sparse.eye(ny)
    Dx = sparse.diags([off_x, main, off_x], offsets=[-1, 0, 1], shape=(nx, nx))
    Dy = sparse.diags([off_y, 0.0, off_y], offsets=[-1, 0, 1], shape=(ny, ny))

    laplace = sparse.kron(Dx, Iy) + sparse.kron(Ix, Dy)
    A = laplace + wavenumber**2 * sparse.eye(nx * ny)

    return _apply_dirichlet(A, grid.shape)


def _assemble_3d(grid: GridSpec, wavenumber: float) -> sparse.spmatrix:
    nx, ny, nz = grid.shape
    hx, hy, hz = grid.spacing
    off_x = -1.0 / hx**2
    off_y = -1.0 / hy**2
    off_z = -1.0 / hz**2
    main = 2.0 * (1.0 / hx**2 + 1.0 / hy**2 + 1.0 / hz**2) + wavenumber**2

    Ix = sparse.eye(nx)
    Iy = sparse.eye(ny)
    Iz = sparse.eye(nz)

    Dx = sparse.diags([off_x, main, off_x], offsets=[-1, 0, 1], shape=(nx, nx))
    Dy = sparse.diags([off_y, 0.0, off_y], offsets=[-1, 0, 1], shape=(ny, ny))
    Dz = sparse.diags([off_z, 0.0, off_z], offsets=[-1, 0, 1], shape=(nz, nz))

    laplace = (
        sparse.kron(sparse.kron(Dx, Iy), Iz)
        + sparse.kron(sparse.kron(Ix, Dy), Iz)
        + sparse.kron(sparse.kron(Ix, Iy), Dz)
    )
    A = laplace + wavenumber**2 * sparse.eye(nx * ny * nz)
    return _apply_dirichlet(A, grid.shape)


def _apply_dirichlet(matrix: sparse.spmatrix, shape: Tuple[int, ...]) -> sparse.spmatrix:
    matrix = matrix.tolil()
    def flatten_index(idx):
        stride = 1
        linear = 0
        for size, i in zip(reversed(shape), reversed(idx)):
            linear += i * stride
            stride *= size
        return linear

    iterator = np.ndindex(shape)
    for idx in iterator:
        if any(i == 0 or i == dim - 1 for i, dim in zip(idx, shape)):
            row = flatten_index(idx)
            matrix.rows[row] = [row]
            matrix.data[row] = [1.0]
    return matrix.tocsr()


def _ensure_sparse() -> None:
    if sparse is None:
        raise ImportError(
            "scipy.sparse is required for operator assembly"  # noqa: EM101
        ) from _sparse_import_error
