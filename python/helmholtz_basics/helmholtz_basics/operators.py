"""Discretisation strategies for the Helmholtz operator."""

from __future__ import annotations

from dataclasses import dataclass
from functools import reduce
from typing import Iterable, Protocol, Sequence

import numpy as np

try:  # pragma: no cover - optional dependency
    from scipy import sparse
except Exception as exc:  # pragma: no cover
    sparse = None  # type: ignore[assignment]
    _sparse_import_error = exc
else:  # pragma: no cover
    _sparse_import_error = None

from .config import GridSpec


class Discretisation(Protocol):
    """Protocol describing how to assemble the Helmholtz operator."""

    name: str

    def assemble(self, grid: GridSpec, wavenumber: float) -> "sparse.spmatrix":
        """Return the discrete Helmholtz operator."""


@dataclass(frozen=True)
class FiniteDifference:
    """Second-order central finite-difference discretisation."""

    name: str = "finite_difference"

    def assemble(self, grid: GridSpec, wavenumber: float) -> "sparse.spmatrix":
        _ensure_sparse()
        laplace = _laplacian(grid.shape, grid.spacing)
        helmholtz = laplace + (wavenumber**2) * sparse.eye(grid.size, format="csr")
        return _apply_dirichlet(helmholtz, grid.shape)


def assemble_operator(
    discretisation: Discretisation, grid: GridSpec, wavenumber: float
) -> "sparse.spmatrix":
    """Public helper to match the package style."""

    _ensure_sparse()
    return discretisation.assemble(grid, wavenumber)


def _laplacian(shape: Sequence[int], spacing: Sequence[float]) -> "sparse.spmatrix":
    components = []
    identity_cache = [sparse.eye(n, format="csr") for n in shape]
    for axis, (n, h) in enumerate(zip(shape, spacing)):
        second = _second_derivative_matrix(n, h)
        kron_factors = list(identity_cache)
        kron_factors[axis] = second
        components.append(_kron_all(kron_factors))
    laplace = reduce(lambda acc, term: acc + term, components)
    return laplace.tocsr()


def _second_derivative_matrix(n: int, h: float) -> "sparse.spmatrix":
    main = np.full(n, -2.0 / h**2)
    off = np.full(n - 1, 1.0 / h**2)
    return sparse.diags([off, main, off], offsets=[-1, 0, 1], format="csr")


def _kron_all(matrices: Iterable["sparse.spmatrix"]) -> "sparse.spmatrix":
    matrices = list(matrices)
    return reduce(lambda a, b: sparse.kron(a, b, format="csr"), matrices)


def _apply_dirichlet(matrix: "sparse.spmatrix", shape: Sequence[int]) -> "sparse.spmatrix":
    matrix = matrix.tolil()
    iterator = np.ndindex(tuple(shape))
    for idx in iterator:
        if _on_boundary(idx, shape):
            row = _flatten_index(idx, shape)
            matrix.rows[row] = [row]
            matrix.data[row] = [1.0]
    return matrix.tocsr()


def _on_boundary(idx: Sequence[int], shape: Sequence[int]) -> bool:
    return any(i == 0 or i == dim - 1 for i, dim in zip(idx, shape))


def _flatten_index(idx: Sequence[int], shape: Sequence[int]) -> int:
    linear = 0
    stride = 1
    for i, dim in zip(reversed(idx), reversed(shape)):
        linear += i * stride
        stride *= dim
    return int(linear)


def _ensure_sparse() -> None:
    if sparse is None:  # pragma: no cover - handled above
        raise ImportError(
            "scipy.sparse is required for operator assembly"  # noqa: EM101
        ) from _sparse_import_error