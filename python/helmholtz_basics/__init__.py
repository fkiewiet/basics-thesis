"""Python companion package for modular Helmholtz GMRES experiments."""

from .config import GridSpec, HelmholtzConfig, SweepConfig
from .grid import build_grid
from .operators import (
    Discretisation,
    FiniteDifference,
    assemble_operator,
)
from .loads import Load, PointSource, PlaneWaveSource, RandomSource, build_load
from .solvers import GMRESOptions, SolverResult, gmres_solve
from .experiments import run_experiment

__all__ = [
    "GridSpec",
    "HelmholtzConfig",
    "SweepConfig",
    "build_grid",
    "Discretisation",
    "FiniteDifference",
    "assemble_operator",
    "Load",
    "PointSource",
    "PlaneWaveSource",
    "RandomSource",
    "build_load",
    "SolverResult",
    "GMRESOptions",
    "gmres_solve",
    "run_experiment",
]
