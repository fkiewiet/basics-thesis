"""Configuration containers for Helmholtz simulations."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Dict, Iterable, Optional, Sequence, Tuple


@dataclass(frozen=True)
class GridSpec:
    """Describe the computational domain and resolution."""

    dims: int
    shape: Tuple[int, ...]
    lengths: Tuple[float, ...]

    def __post_init__(self) -> None:
        if self.dims < 1 or self.dims > 3:
            raise ValueError("Only 1D, 2D, or 3D grids are supported")
        if len(self.shape) != self.dims:
            raise ValueError("shape must provide one entry per dimension")
        if len(self.lengths) != self.dims:
            raise ValueError("lengths must provide one entry per dimension")
        if any(n <= 1 for n in self.shape):
            raise ValueError("each grid dimension must have more than one node")
        if any(L <= 0 for L in self.lengths):
            raise ValueError("domain lengths must be positive")

    @property
    def spacing(self) -> Tuple[float, ...]:
        return tuple(L / (n - 1) for L, n in zip(self.lengths, self.shape))


@dataclass(frozen=True)
class BoundarySpec:
    """Placeholder for boundary conditions."""

    kind: str = "dirichlet"
    value: float = 0.0


@dataclass(frozen=True)
class HelmholtzConfig:
    """Collects the ingredients for a single solve."""

    wavenumber: float
    grid: GridSpec
    load: "Load"
    discretisation: "Discretisation"
    boundary: BoundarySpec = BoundarySpec()


@dataclass(frozen=True)
class SweepConfig:
    """Describe parameter sweeps for experiments."""

    grids: Sequence[GridSpec]
    wavenumbers: Sequence[float]
    loads: Sequence["Load"]
    discretisations: Sequence["Discretisation"]
    metadata: Dict[str, object] = field(default_factory=dict)
    hook: Optional[Callable[["HelmholtzConfig", "SolverResult"], None]] = None

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
                            boundary=BoundarySpec(),
                        )


# Type checking support via forward references
if False:  # pragma: no cover
    from .loads import Load
    from .operators import Discretisation
    from .solvers import SolverResult
