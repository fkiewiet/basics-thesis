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
        spacing = grid.spacing
        if len(self.location) < grid.dims:
            raise ValueError("PointSource.location must supply one coordinate per dimension")
        physical = self.location[: grid.dims]
        idx = []
        for coord, h, n in zip(physical, spacing, grid.shape):
            lattice = coord / h
            idx.append(int(round(np.clip(lattice, 0.0, n - 1))))
        rhs = np.zeros(grid.size, dtype=float)
        rhs[_flatten_index(tuple(idx), grid.shape)] = self.amplitude
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