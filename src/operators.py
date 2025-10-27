"""
operators.py — Finite-difference operators for Helmholtz/Laplacian
with Dirichlet BC and optional PML via complex coordinate stretching.

Sign convention: we assemble the discrete operator for
    (-Δ - k^2) u = f
on a rectangular domain with uniform spacing per axis.

PML (optional):
We implement a diagonal complex stretching s_j(x_j) per axis j.
In PML mode, the operator is
    Σ_j ∂/∂x_j [ (1 / s_j) ∂u/∂x_j ]  +  k^2 S0 u,
where S0 = Π_j s_j. We discretize each axis term in divergence form
using variable-coefficient second differences with Dirichlet at the
physical boundary (the outer grid nodes are still unknowns; the PML
attenuation handles absorption). If you want hard Dirichlet walls
without PML, keep pml=None.

Notes
-----
- s_j(x_j) := 1 + i * σ_j(x_j) / k   (dimensionless)
- σ_j ramps from 0 inside the domain to σ_max at the outer edge of the
  PML layer with polynomial order m (default m=3).
- If your physical model requires a different scaling (e.g., σ/ω),
  adjust `stretch_profile()` accordingly.

Public API
----------
- laplacian_operator(shape, lengths, dtype)
- helmholtz_operator(shape, lengths, k, dtype, pml=None)
- FiniteDifference(...).assemble(kind="laplacian"|"helmholtz", k=None, pml=None)

Minimal dependencies: NumPy + SciPy (sparse).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Sequence, Tuple, List

import numpy as np
from scipy.sparse import diags, kron, eye, csc_matrix, spdiags
from scipy.sparse import issparse
from numpy.typing import NDArray


# ----------------------------
# Configuration & PML settings
# ----------------------------

@dataclass(frozen=True)
class FDConfig:
    """Finite-difference configuration."""
    dtype: np.dtype = np.complex128  # complex dtype recommended for Helmholtz/PML


@dataclass(frozen=True)
class PMLConfig:
    """
    Perfectly Matched Layer configuration.

    thickness: number of grid cells per axis used as PML (at *each* boundary).
    m: polynomial order of the damping ramp (>= 1).
    sigma_max: maximum damping at the outer boundary of PML.
    """
    thickness: int = 10
    m: int = 3
    sigma_max: float = 50.0  # increase for stronger absorption


# ----------------------------
# Public functional API
# ----------------------------

def laplacian_operator(
    shape: Sequence[int],
    lengths: Sequence[float],
    dtype: np.dtype = np.complex128,
) -> csc_matrix:
    """Assemble -Δ with Dirichlet BC (no PML)."""
    fd = FiniteDifference(FDConfig(dtype=dtype))
    return fd.assemble(kind="laplacian", shape=shape, lengths=lengths)


def helmholtz_operator(
    shape: Sequence[int],
    lengths: Sequence[float],
    k: float,
    dtype: np.dtype = np.complex128,
    pml: Optional[PMLConfig] = None,
) -> csc_matrix:
    """Assemble (-Δ - k^2) with Dirichlet BC, or PML if `pml` is given."""
    fd = FiniteDifference(FDConfig(dtype=dtype))
    return fd.assemble(kind="helmholtz", shape=shape, lengths=lengths, k=k, pml=pml)


# ----------------------------
# Core class
# ----------------------------

class FiniteDifference:
    """
    Finite-difference assembly on a tensor-product grid with uniform spacing per axis.

    Methods:
      - assemble(kind="laplacian"|"helmholtz", k=None, pml=None)
    """

    def __init__(self, config: FDConfig = FDConfig()):
        self.cfg = config

    # ---------- public ----------

    def assemble(
        self,
        *,
        kind: str,
        shape: Sequence[int],
        lengths: Sequence[float],
        k: Optional[float] = None,
        pml: Optional[PMLConfig] = None,
    ) -> csc_matrix:
        shape = tuple(int(n) for n in shape)
        lengths = tuple(float(L) for L in lengths)
        assert len(shape) == len(lengths), "shape and lengths must have same dimensionality"
        assert all(n >= 2 for n in shape), "each axis must have at least 2 points"

        if kind not in {"laplacian", "helmholtz"}:
            raise ValueError("kind must be 'laplacian' or 'helmholtz'")

        if pml is None:
            L = self._assemble_laplacian_dirichlet(shape, lengths).astype(self.cfg.dtype)
            if kind == "laplacian":
                return L.tocsc()
            if k is None:
                raise ValueError("helmholtz assembly requires k")
            N = np.prod(shape)
            M = eye(N, dtype=self.cfg.dtype, format="csc")
            return (L - (k**2) * M).tocsc()

        # PML path
        if kind == "laplacian":
            # interpret as Σ_j ∂/∂x_j [(1/s_j) ∂/∂x_j] with k only used to scale s_j (use k=1)
            k_eff = 1.0
            A = self._assemble_pml_helmholtz_like(shape, lengths, k=k_eff, pml=pml)
            return A.tocsc()
        else:
            if k is None:
                raise ValueError("helmholtz assembly requires k")
            A = self._assemble_pml_helmholtz_like(shape, lengths, k=k, pml=pml)
            return A.tocsc()

    # ---------- Dirichlet core (no PML) ----------

    def _assemble_laplacian_dirichlet(
        self, shape: Tuple[int, ...], lengths: Tuple[float, ...]
    ) -> csc_matrix:
        """Kronecker sum of 1D Dirichlet second-difference operators (with minus sign)."""
        dims = len(shape)
        mats: List[csc_matrix] = []
        for ax in range(dims):
            n = shape[ax]
            h = lengths[ax] / (n - 1)
            L1 = self._second_derivative_1d_dirichlet(n, h)  # returns -d2/dx2
            mats.append(L1.astype(self.cfg.dtype))

        # Kronecker-sum: L = Lx ⊕ Ly (⊕ Lz ...) = kron(I, Ly) + kron(Lx, I) + ...
        L = None
        for ax, Lax in enumerate(mats):
            Ia = 1
            for b in range(dims):
                if b == ax:
                    continue
                Ib = eye(shape[b], dtype=self.cfg.dtype, format="csc")
                Ia = kron(Ia, Ib, format="csc") if not isinstance(Ia, csc_matrix) else kron(Ia, Ib, format="csc")
            term = kron(Ia, Lax, format="csc") if dims == 1 else self._place_axis(Lax, shape, ax)
            L = term if L is None else (L + term)
        return L.tocsc()

    @staticmethod
    def _second_derivative_1d_dirichlet(n: int, h: float) -> csc_matrix:
        """
        Return the 1D operator for -d2/dx2 with Dirichlet BC on both ends.
        Shape: (n x n). Interior 3-point stencil; boundaries have the usual 2-point rows.
        """
        main = 2.0 * np.ones(n)
        off = -1.0 * np.ones(n - 1)
        D2 = diags([off, main, off], offsets=[-1, 0, 1], shape=(n, n), format="csc")
        return (D2 / (h**2)).astype(np.complex128)

    def _place_axis(self, Lax: csc_matrix, shape: Tuple[int, ...], ax: int) -> csc_matrix:
        """
        Embed a 1D operator Lax (size n_ax) into the full ND space via Kronecker products,
        acting along axis `ax`.
        """
        dims = len(shape)
        result = None
        for d in range(dims):
            Id = eye(shape[d], dtype=self.cfg.dtype, format="csc") if d != ax else Lax
            result = Id if result is None else kron(result, Id, format="csc")
        return result.tocsc()

    # ---------- PML assembly ----------

    def _assemble_pml_helmholtz_like(
        self,
        shape: Tuple[int, ...],
        lengths: Tuple[float, ...],
        *,
        k: float,
        pml: PMLConfig,
    ) -> csc_matrix:
        """
        Assemble Σ_j ∂/∂x_j [(1/s_j) ∂u/∂x_j]  +  k^2 S0 u  with Dirichlet-style
        stencil inside, complex stretching s_j in the PML collar.

        Discretization per axis j:
          L_j u ≈ -(1/h_j^2) * [ a_{i+1/2} (u_{i+1} - u_i) - a_{i-1/2} (u_i - u_{i-1}) ],
        where a = 1 / s_j and s_j varies only along axis j. Boundaries use the
        usual reduced stencil (Dirichlet style row).
        """
        dims = len(shape)
        # Build stretching arrays per axis
        s_axes: List[NDArray[np.complex128]] = []
        for ax in range(dims):
            s_axes.append(self._stretch_axis(shape[ax], lengths[ax], k=k, pml=pml))

        # Axis operators (variable-coefficient 1D) placed via Kronecker
        L_sum = None
        for ax in range(dims):
            n = shape[ax]
            h = lengths[ax] / (n - 1)
            s = s_axes[ax]  # shape (n,)
            a = 1.0 / s     # coefficient for divergence-form diffusion

            L1 = self._varcoef_second_derivative_dirichlet(a, h)  # csc (n x n)

            # Place along axis ax
            Term_ax = self._place_axis(L1.astype(self.cfg.dtype), shape, ax)
            L_sum = Term_ax if L_sum is None else (L_sum + Term_ax)

        # Mass term: k^2 * S0 * I, where S0 = Π_j s_j  (pointwise).
        # Realistically, S0 varies across the grid; implement as diagonal mass.
        # Construct S0 as Kronecker sum of per-axis diag(s_j) expanded across other axes.
        S0_diag = self._tensor_product_pointwise_product(s_axes)  # (N,) complex
        M = spdiags(S0_diag.astype(self.cfg.dtype), 0, S0_diag.size, S0_diag.size).tocsc()

        return (L_sum + (k**2) * M).tocsc()

    @staticmethod
    def _varcoef_second_derivative_dirichlet(a: NDArray[np.complex128], h: float) -> csc_matrix:
        """
        Variable-coefficient -d/dx ( a(x) d/dx ) with Dirichlet rows at boundaries.
        a is sampled at **nodes**; we form face-centered coefficients by arithmetic mean.
        """
        n = a.size
        # Face-centered a_{i+1/2}
        a_face = 0.5 * (a[:-1] + a[1:])  # size n-1

        # Tridiagonal entries
        main = np.zeros(n, dtype=np.complex128)
        upper = np.zeros(n - 1, dtype=np.complex128)
        lower = np.zeros(n - 1, dtype=np.complex128)

        # Interior i=1..n-2
        main[1:-1] = (a_face[:-1] + a_face[1:])
        upper[:-1] = -a_face[:-1]
        lower[:-1] = -a_face[:-1]
        # last upper/lower entry
        upper[-1] = -a_face[-1]
        lower[-1] = -a_face[-1]

        # Dirichlet rows (i=0, i=n-1) -> standard 1D Lap row with a_face at the single interior face
        main[0] = a_face[0]
        main[-1] = a_face[-1]

        D = diags([lower, main, upper], offsets=[-1, 0, 1], shape=(n, n), format="csc")
        return (D / (h**2)).astype(np.complex128)

    def _stretch_axis(
        self, n: int, L: float, *, k: float, pml: PMLConfig
    ) -> NDArray[np.complex128]:
        """
        Build complex stretching s(x) = 1 + i * σ(x) / k over a 1D axis.
        σ ramps to σ_max within PML thickness at both ends; zero in interior.
        """
        # grid coordinates in [0, L], uniform spacing
        x = np.linspace(0.0, L, n)
        h = x[1] - x[0]

        t = pml.thickness
        m = pml.m
        sig_max = pml.sigma_max

        sigma = np.zeros(n, dtype=float)

        # Left PML: nodes 0..t-1
        if t > 0:
            xi = np.clip((t - 1 - np.arange(t)) / max(t - 1, 1), 0.0, 1.0)  # 1 at boundary, 0 at interior edge
            ramp = (xi ** m)
            sigma[:t] = sig_max * ramp

            # Right PML: nodes n-t..n-1
            xi_r = np.clip((np.arange(t)) / max(t - 1, 1), 0.0, 1.0)
            ramp_r = (xi_r ** m)
            sigma[-t:] = np.maximum(sigma[-t:], sig_max * ramp_r)

        s = 1.0 + 1j * (sigma / max(k, 1e-12))
        return s.astype(np.complex128)

    # ---------- helpers ----------

    def _tensor_product_pointwise_product(self, axes_arrays: List[NDArray[np.complex128]]) -> NDArray[np.complex128]:
        """
        Given per-axis arrays a_j of shapes (n_j,), return the flattened (N,) array of
        Π_j a_j on the tensor grid in C-order (row-major) consistent with our Kronecker layout.
        """
        # Meshgrid in 'ij' indexing then flatten in C-order
        grids = np.meshgrid(*axes_arrays, indexing='ij')
        prod = np.ones_like(grids[0], dtype=np.complex128)
        for g in grids:
            prod = prod * g
        return prod.ravel(order='C')


# ----------------------------
# End of module
# ----------------------------


# ----------------------------
# Backward-compatibility shim
# ----------------------------
from enum import Enum

class BC(Enum):
    """Legacy placeholder (Dirichlet only). Kept for backward compatibility."""
    DIRICHLET = "dirichlet"

class Discretisation:
    """Legacy no-op placeholder so old imports don't crash."""
    pass

def assemble_operator(
    # legacy signature variants supported:
    disc=None,
    grid=None,
    k: float | None = None,
    *,
    kind: str = "helmholtz",
    shape: Sequence[int] | None = None,
    lengths: Sequence[float] | None = None,
    dtype: np.dtype = np.complex128,
    bc: BC | None = None,            # ignored; Dirichlet is default in new path
    pml: PMLConfig | None = None,
):
    """
    Back-compat wrapper that maps old calls to the new API.
    Accepted call styles:
        assemble_operator(disc, grid, k, kind="helmholtz", ...)
        assemble_operator(kind="laplacian", shape=..., lengths=..., dtype=..., pml=...)
    """
    # Extract shape/lengths from grid if provided
    if grid is not None:
        shape = getattr(grid, "shape", shape)
        lengths = getattr(grid, "lengths", lengths)

    if shape is None or lengths is None:
        raise ValueError("assemble_operator: 'shape' and 'lengths' must be provided (directly or via 'grid').")

    fd = FiniteDifference(FDConfig(dtype=dtype))

    kind = kind.lower()
    if kind == "laplacian":
        # PML path allowed: builds Σ_j ∂x_j(1/s_j ∂x_j) when pml is not None
        return fd.assemble(kind="laplacian", shape=shape, lengths=lengths, pml=pml)
    elif kind == "helmholtz":
        if k is None:
            raise ValueError("assemble_operator(kind='helmholtz') requires k.")
        return fd.assemble(kind="helmholtz", shape=shape, lengths=lengths, k=k, pml=pml)
    else:
        raise ValueError("assemble_operator: kind must be 'laplacian' or 'helmholtz'.")
