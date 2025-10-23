"""Factory functions for Helmholtz right-hand sides."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Protocol, Tuple, Union, Any

import numpy as np

from .config import GridSpec


class Load(Protocol):
    """Protocol for building Helmholtz right-hand sides."""

    name: str

    def build(self, grid: GridSpec) -> np.ndarray:
        """Return a flattened array matching the grid size."""


@dataclass(frozen=True)
class PointSource:
    """Load representing a single impulse on the grid.

    Parameters
    ----------
    location : Tuple[float, ...] | str
        Physical coordinates (one per dimension) OR a keyword:
        - "centre"/"center"/"middle"
        - "origin"/"zero"
        - "random" (uses `self.seed` if present for reproducible RNG)
    amplitude : float
        Magnitude of the point source.
    phase : float
        Phase in radians (the complex value written is amplitude * exp(1j * phase)).
    name : str
        Label used in experiment summaries/plots.
    """

    location: Union[Tuple[float, ...], str]
    amplitude: float = 1.0
    phase: float = 0.0      # <-- NEW
    name: str = "point_source"

    def build(self, grid: GridSpec) -> np.ndarray:
        """
        Return a RHS vector with a unit point load at the requested location
        (scaled by `self.amplitude`).

        `self.location` may be:
        - an iterable of physical coordinates (len == grid.dims), or
        - a keyword string: "centre"/"center"/"middle", "origin"/"zero", "random".
            For "random", a reproducible RNG is used if `self.seed` exists.
        """
        # --- resolve physical coordinates ---
        if isinstance(self.location, str):
            key = self.location.lower()
            if key in ("centre", "center", "middle"):
                physical = tuple(L * 0.5 for L in grid.lengths)
            elif key in ("origin", "zero"):
                physical = tuple(0.0 for _ in range(grid.dims))
            elif key == "random":
                rng = np.random.default_rng(getattr(self, "seed", None))
                physical = tuple(float(rng.uniform(0.0, L)) for L in grid.lengths)
            else:
                raise ValueError(f"Unknown location keyword: {self.location!r}")
        else:
            # Assume iterable of floats
            try:
                physical = tuple(self.location[: grid.dims])  # supports lists/tuples/np arrays
            except Exception as exc:  # noqa: BLE001
                raise TypeError("PointSource.location must be a string or an iterable of floats") from exc
            if len(physical) < grid.dims:
                raise ValueError("PointSource.location must supply one coordinate per dimension")

        # --- map physical coords -> lattice indices (clamped) ---
        idx = []
        for coord, h, n in zip(physical, grid.spacing, grid.shape):
            lattice = float(coord) / float(h)
            idx.append(int(round(np.clip(lattice, 0.0, n - 1))))

        # --- assemble RHS vector ---
        rhs = np.zeros(grid.size, dtype=np.complex128)
        lin = _flatten_index(tuple(idx), grid.shape)
        rhs[lin] = float(self.amplitude) * np.exp(1j * float(self.phase))
        return rhs


@dataclass
class RandomPointSource:
    amplitude_range: Tuple[float, float] = (0.5, 1.5)
    phase_range: Tuple[float, float] = (0.0, 2*np.pi)
    interior_only: bool = True
    seed: int | None = None
    name: str = "random_point_source"   # <-- add this

    _rng: Any = field(default=None, init=False, repr=False)
    _last: dict | None = field(default=None, init=False, repr=False)

    def __post_init__(self):
        self._rng = np.random.default_rng(self.seed)

    def _sample_location(self, shape: Tuple[int, ...]) -> Tuple[int, ...]:
        idx = []
        for n in shape:
            if self.interior_only and n > 2:
                idx.append(int(self._rng.integers(1, n-1)))
            else:
                idx.append(int(self._rng.integers(0, n)))
        return tuple(idx)

    def build(self, grid) -> np.ndarray:
        idx = self._sample_location(grid.shape)
        A = self._rng.uniform(*self.amplitude_range)
        phase = self._rng.uniform(*self.phase_range)
        value = A * np.exp(1j * phase)

        rhs = np.zeros(grid.size, dtype=np.complex128)
        lin = np.ravel_multi_index(idx, grid.shape)
        rhs[lin] = value
        self._last = {"A": A, "phase": phase, "lattice_idx": tuple(map(int, idx))}
        return rhs

@dataclass
class RandomContinuousPointSource:
    """
    Random point source at a continuous (x,y), deposited bilinearly to the grid.
    Keeps amplitude/phase sampling like RandomPointSource.
    """
    amplitude_range: Tuple[float, float] = (0.5, 1.5)
    phase_range: Tuple[float, float] = (0.0, 2 * np.pi)
    interior_margin_cells: int = 1    # keep at least this many cells from each boundary
    seed: int | None = None
    name: str = "random_continuous_point_source"

    _rng: Any = field(default=None, init=False, repr=False)
    _last: dict | None = field(default=None, init=False, repr=False)

    def __post_init__(self):
        self._rng = np.random.default_rng(self.seed)

    def build(self, grid: GridSpec) -> np.ndarray:
        if grid.dims != 2:
            raise ValueError("RandomContinuousPointSource currently supports 2D grids only")
        nx, ny = grid.shape
        hx, hy = grid.spacing
        Lx, Ly = grid.lengths

        # interior margin in *physical* units
        mx = max(0, self.interior_margin_cells) * hx
        my = max(0, self.interior_margin_cells) * hy

        x = float(self._rng.uniform(mx, Lx - mx))
        y = float(self._rng.uniform(my, Ly - my))

        A = float(self._rng.uniform(*self.amplitude_range))
        phase = float(self._rng.uniform(*self.phase_range))
        value = A * np.exp(1j * phase)

        rhs = np.zeros(grid.size, dtype=np.complex128)
        _bilinear_deposit(rhs, grid, x, y, value)

        self._last = {"A": A, "phase": phase, "xy": (x, y)}
        return rhs

def _bilinear_deposit(rhs: np.ndarray, grid: GridSpec, x: float, y: float, value: complex) -> None:
    """Deposit `value` at physical (x,y) onto the 2×2 neighboring nodes with bilinear weights."""
    nx, ny = grid.shape
    hx, hy = grid.spacing

    # convert to lattice coordinates (node i at i*hx)
    u = np.clip(x / hx, 0.0, nx - 1.0)
    v = np.clip(y / hy, 0.0, ny - 1.0)

    i0 = int(np.floor(u));  j0 = int(np.floor(v))
    i1 = min(i0 + 1, nx - 1); j1 = min(j0 + 1, ny - 1)

    tx = float(u - i0); ty = float(v - j0)

    w00 = (1 - tx) * (1 - ty)
    w10 = tx * (1 - ty)
    w01 = (1 - tx) * ty
    w11 = tx * ty

    rhs[_flatten_index((i0, j0), grid.shape)] += w00 * value
    rhs[_flatten_index((i1, j0), grid.shape)] += w10 * value
    rhs[_flatten_index((i0, j1), grid.shape)] += w01 * value
    rhs[_flatten_index((i1, j1), grid.shape)] += w11 * value



def build_load(load, grid):
    """
    Build a right-hand side vector b compatible with `grid`.

    Accepts any object with a `.build(grid)` method that returns either:
      - a 1D array-like of length prod(grid.shape), or
      - a 2D array-like with shape == grid.shape

    Returns
    -------
    np.ndarray
        Flattened, C-contiguous RHS (dtype complex128), length == grid.size.
    """
    if not (hasattr(load, "build") and callable(load.build)):
        raise TypeError(
            f"Unsupported load type: {type(load).__name__}. "
            "Provide an object with a `.build(grid)` method."
        )

    b = load.build(grid)
    b = np.asarray(b)  # zero-copy when possible

    # Normalize shape: accept either flat or grid-shaped inputs
    expected_shape = tuple(getattr(grid, "shape", ()))
    expected_size = int(np.prod(expected_shape)) if expected_shape else b.size

    if b.ndim == 2 and b.shape == expected_shape:
        # grid-shaped -> flatten
        b = b.ravel(order="C")
    elif b.ndim == 1 and b.size == expected_size:
        # already flat -> keep
        pass
    elif b.size == expected_size:
        # different ndim but total size matches -> reshape then flatten
        b = np.reshape(b, expected_shape).ravel(order="C")
    else:
        raise ValueError(
            f"Load produced array with shape {b.shape} (size {b.size}), "
            f"but grid expects size {expected_size} (shape {expected_shape})."
        )

    # Ensure complex dtype without touching non-portable dtypes (like np.complex256)
    if not np.issubdtype(b.dtype, np.complexfloating):
        b = b.astype(np.complex128, copy=False)

    # Finalize memory layout (helps sparse ops / GMRES)
    b = np.ascontiguousarray(b)

    # Optional safety: check for NaN/Inf (comment out if you prefer silent pass-through)
    if not np.all(np.isfinite(b)):
        n_bad = np.size(b) - np.count_nonzero(np.isfinite(b))
        raise ValueError(f"RHS contains {n_bad} non-finite entries (NaN/Inf).")

    return b


def _flatten_index(idx: Tuple[int, ...], shape: Tuple[int, ...]) -> int:
    linear = 0
    stride = 1
    for i, size in zip(reversed(idx), reversed(shape)):
        linear += int(i) * stride
        stride *= size
    return linear




@dataclass(frozen=True)
class PlaneWaveSource:
    """Complex plane wave forcing term."""

    direction: Tuple[float, ...]
    phase: float = 0.0
    name: str = "plane_wave"

    def build(self, grid: GridSpec) -> np.ndarray:
        axes = [np.linspace(0.0, L, n, dtype=float) for L, n in zip(grid.lengths, grid.shape)]
        mesh = np.meshgrid(*axes, indexing="ij")
        if len(self.direction) < grid.dims:
            raise ValueError("PlaneWaveSource.direction must have one entry per dimension")
        direction = np.array(self.direction[: grid.dims], dtype=float)
        norm = np.linalg.norm(direction)
        if norm == 0.0:
            raise ValueError("PlaneWaveSource.direction must be non-zero")
        direction /= norm
        phase_field = sum(axis * component for axis, component in zip(mesh, direction))
        values = np.exp(1j * (phase_field + self.phase))
        return values.reshape(-1)