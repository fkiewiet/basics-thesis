"""Configuration containers for Helmholtz simulations."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Iterable, Mapping, MutableMapping, Optional, Sequence, Tuple

if False:  # pragma: no cover - imported for typing only
    from .loads import Load
    from .operators import Discretisation
    from .solvers import SolverResult


@dataclass(frozen=True)
class GridSpec:
    """Describe the computational domain and resolution."""

    dims: int
    shape: Tuple[int, ...]
    lengths: Tuple[float, ...]

    def __post_init__(self) -> None:
        if not (1 <= self.dims <= 3):
            raise ValueError("Only 1D, 2D, or 3D grids are supported")
        if len(self.shape) != self.dims:
            raise ValueError("shape must match the number of dimensions")
        if len(self.lengths) != self.dims:
            raise ValueError("lengths must match the number of dimensions")
        if any(n <= 1 for n in self.shape):
            raise ValueError("each dimension needs at least two nodes")
        if any(L <= 0 for L in self.lengths):
            raise ValueError("domain lengths must be positive")

    @property
    def spacing(self) -> Tuple[float, ...]:
        """Return the mesh spacing for each dimension."""

        return tuple(L / (n - 1) for L, n in zip(self.lengths, self.shape))

    @property
    def size(self) -> int:
        """Total number of grid points."""

        total = 1
        for n in self.shape:
            total *= n
        return total


@dataclass(frozen=True)
class HelmholtzConfig:
    """Collect the ingredients for a single solve."""

    wavenumber: float
    grid: GridSpec
    load: "Load"
    discretisation: "Discretisation"
    label: Optional[str] = None
    metadata: Mapping[str, object] = field(default_factory=dict)


@dataclass(frozen=True)
class SweepConfig:
    """Describe a Cartesian product of parameters for experiments."""

    grids: Sequence[GridSpec]
    wavenumbers: Sequence[float]
    loads: Sequence["Load"]
    discretisations: Sequence["Discretisation"]
    hook: Optional[Callable[[HelmholtzConfig, "SolverResult"], None]] = None
    annotations: MutableMapping[str, object] = field(default_factory=dict)

    def iter_configs(self) -> Iterable[HelmholtzConfig]:
        for grid in self.grids:
            for k in self.wavenumbers:
                for load in self.loads:
                    for disc in self.discretisations:
                        yield HelmholtzConfig(
                            wavenumber=k,
                            grid=grid,
                            load=load,
                            discretisation=disc,
                            metadata=self.annotations,
                        )