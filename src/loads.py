"""Factory functions for Helmholtz right-hand sides."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol, Tuple

import numpy as np

from .config import GridSpec


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


@dataclass(frozen=True)
class RandomSource:
    """Gaussian random forcing for stress-testing solvers."""

    seed: int | None = None
    name: str = "random"

    def build(self, grid: GridSpec) -> np.ndarray:
        rng = np.random.default_rng(self.seed)
        return rng.standard_normal(grid.size)


def build_load(load: Load, grid: GridSpec) -> np.ndarray:
    """Helper mirroring the protocol."""

    return load.build(grid)


def _flatten_index(idx: Tuple[int, ...], shape: Tuple[int, ...]) -> int:
    linear = 0
    stride = 1
    for i, size in zip(reversed(idx), reversed(shape)):
        linear += int(i) * stride
        stride *= size
    return linear