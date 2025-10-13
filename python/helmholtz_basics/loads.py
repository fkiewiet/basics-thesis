"""Factory functions for Helmholtz right-hand sides."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol, Tuple

import numpy as np

from .config import GridSpec


class Load(Protocol):
    name: str

    def build(self, grid: GridSpec) -> np.ndarray:
        """Return a vectorised load with one entry per grid node."""


@dataclass(frozen=True)
class PointSource:
    location: Tuple[float, ...]
    amplitude: float = 1.0
    name: str = "point_source"

    def build(self, grid: GridSpec) -> np.ndarray:
        idx = [int(round(x / (L / (n - 1)))) for x, L, n in zip(self.location, grid.lengths, grid.shape)]
        idx = [min(max(i, 0), n - 1) for i, n in zip(idx, grid.shape)]
        rhs = np.zeros(np.prod(grid.shape), dtype=float)
        rhs[_flatten_index(tuple(idx), grid.shape)] = self.amplitude
        return rhs


@dataclass(frozen=True)
class PlaneWaveSource:
    direction: Tuple[float, ...]
    phase: float = 0.0
    name: str = "plane_wave"

    def build(self, grid: GridSpec) -> np.ndarray:
        coords = np.stack(
            np.meshgrid(
                *[np.linspace(0.0, L, n) for L, n in zip(grid.lengths, grid.shape)],
                indexing="ij",
            ),
            axis=-1,
        )
        direction = np.array(self.direction[: grid.dims], dtype=float)
        direction /= np.linalg.norm(direction)
        phase_field = np.tensordot(coords[..., : grid.dims], direction, axes=([coords.ndim - 1], [0]))
        return np.exp(1j * (phase_field + self.phase)).reshape(-1)


@dataclass(frozen=True)
class RandomSource:
    seed: int | None = None
    name: str = "random"

    def build(self, grid: GridSpec) -> np.ndarray:
        rng = np.random.default_rng(self.seed)
        rhs = rng.standard_normal(np.prod(grid.shape))
        return rhs


def build_load(load: Load, grid: GridSpec) -> np.ndarray:
    return load.build(grid)


def _flatten_index(idx: Tuple[int, ...], shape: Tuple[int, ...]) -> int:
    linear = 0
    stride = 1
    for i, size in zip(reversed(idx), reversed(shape)):
        linear += i * stride
        stride *= size
    return linear
