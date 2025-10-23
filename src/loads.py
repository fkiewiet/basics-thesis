"""Factory functions for Helmholtz right-hand sides."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol, Tuple

import numpy as np

from .config import GridSpec


from typing import Tuple, Any


class Load(Protocol):
    """Protocol for building Helmholtz right-hand sides."""

    name: str

    def build(self, grid: GridSpec) -> np.ndarray:
        """Return a flattened array matching the grid size."""


@dataclass(frozen=True)
class PointSource:
    """Load representing a single impulse on the grid."""

    location: Tuple[float, ...]
    amplitude: float = 1.0
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
        rhs = np.zeros(grid.size, dtype=float)
        rhs[_flatten_index(tuple(idx), grid.shape)] = float(self.amplitude)
        return rhs


@dataclass
class RandomPointSource:
    """
    Randomized point source with random amplitude, phase, and lattice location.

    Parameters
    ----------
    amplitude_range : (float, float)
        Uniform range for |A| (inclusive of endpoints in NumPy sense).
    phase_range : (float, float)
        Uniform range for phase in radians.
    interior_only : bool
        If True and a dimension has size > 2, sample indices in [1, n-2] to avoid boundaries.
    seed : int | None
        RNG seed for reproducibility. If None, uses a non-deterministic seed.
    """
    amplitude_range: Tuple[float, float] = (0.5, 1.5)
    phase_range: Tuple[float, float] = (0.0, 2*np.pi)
    interior_only: bool = True
    seed: int | None = None

    # Internal state (not in signature)
    _rng: Any = None
    _last: dict | None = None

    def __post_init__(self):
        self._rng = np.random.default_rng(self.seed)

    def _sample_location(self, shape: Tuple[int, ...]) -> Tuple[int, ...]:
        idx = []
        for n in shape:
            if self.interior_only and n > 2:
                i = int(self._rng.integers(1, n - 1))
            else:
                i = int(self._rng.integers(0, n))
            idx.append(i)
        return tuple(idx)

    def build(self, grid) -> np.ndarray:
        """
        Build a flattened complex RHS vector b with a single complex spike at a random index.
        b[idx] = A * exp(1j * phase), with A, phase, and idx sampled randomly.
        """
        # Prefer grid.shape; fall back to grid.size as 1D
        shape = getattr(grid, "shape", None)
        if shape is None:
            size = getattr(grid, "size", None)
            if size is None:
                raise ValueError("Grid must expose .shape or .size")
            shape = (int(size),)

        size = int(getattr(grid, "size", int(np.prod(shape))))
        b = np.zeros(size, dtype=np.complex128)

        # Draw amplitude and phase
        A = float(self._rng.uniform(*self.amplitude_range))
        phase = float(self._rng.uniform(*self.phase_range))

        # Pick random lattice index and map to linear
        lattice_idx = self._sample_location(tuple(int(s) for s in shape))
        lin_idx = int(np.ravel_multi_index(lattice_idx, shape))

        b[lin_idx] = A * np.exp(1j * phase)

        # Keep last sample for debugging/labels
        self._last = {
            "A": A,
            "phase": phase,
            "lattice_idx": lattice_idx,
            "lin_idx": lin_idx,
        }
        return b


# --- OPTIONAL: if you have a build_load helper, ensure it accepts any object with .build ---
def build_load(load, grid):
    """
    Build a right-hand side vector `b` compatible with `grid`.

    Accepts:
      - Known load classes (e.g., PointSource, PlaneWaveSource, RandomPointSource, ...)
      - Any object with a `.build(grid)` method returning a flat array of length grid.size

    Returns
    -------
    np.ndarray
        Flattened RHS vector (complex128 for safety).
    """
    if hasattr(load, "build") and callable(load.build):
        b = load.build(grid)
    else:
        raise TypeError(
            f"Unsupported load type: {type(load).__name__}. "
            "Provide an object with a `.build(grid)` method."
        )
    b = np.asarray(b)
    # Ensure complex dtype to support complex loads (plane waves, randomized phase)
    if b.dtype not in (np.complex64, np.complex128, np.complex256):
        b = b.astype(np.complex128, copy=False)
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