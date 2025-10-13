"""Python companion package for modular Helmholtz GMRES experiments."""

# Optional version exposure (works if the package is installed; falls back locally)
try:
    from importlib import metadata as _metadata  # Python 3.8+
    __version__ = _metadata.version("python-helmholtz-basics")
except Exception:  # not installed / local dev
    __version__ = "0+local"

# Public API (relative imports because this file lives inside the package)
from .config import GridSpec, HelmholtzConfig, SweepConfig
from .grid import build_grid
from .operators import Discretisation, FiniteDifference, assemble_operator
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

# Helpful guard if someone tries `python path/to/__init__.py`
if __name__ == "__main__":
    raise SystemExit(
        "This package cannot be executed directly. "
        "Use it from Python, e.g.: "
        "python -c \"from python_helmholtz_basics import run_experiment; run_experiment(...)\""
    )
