"""Python companion package for modular Helmholtz GMRES experiments.

The imports below prefer relative resolution (normal package usage) but fall back
to absolute imports when the module is executed as a script. The latter guards
against the confusing ``ParserError`` messages that appear if someone runs the
file directly in a shell instead of importing it from Python.
"""

try:  # pragma: no cover - import plumbing
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
except ImportError:  # pragma: no cover - direct execution fallback
    from helmholtz_basics.config import GridSpec, HelmholtzConfig, SweepConfig
    from helmholtz_basics.grid import build_grid
    from helmholtz_basics.operators import (
        Discretisation,
        FiniteDifference,
        assemble_operator,
    )
    from helmholtz_basics.loads import (
        Load,
        PointSource,
        PlaneWaveSource,
        RandomSource,
        build_load,
    )
    from helmholtz_basics.solvers import GMRESOptions, SolverResult, gmres_solve
    from helmholtz_basics.experiments import run_experiment

if __name__ == "__main__":  # pragma: no cover - friendly guidance
    import textwrap

    print(
        textwrap.dedent(
            """
            This file re-exports the main helpers from ``helmholtz_basics``.

            To explore the toolkit interactively, start a Python session and import
            the package, for example:

                python -m pip install -e .  # optional, if you want editable installs
                python
                >>> import helmholtz_basics as hb
                >>> hb.GridSpec(dim=2, shape=(64, 64), spacing=0.02)

            Running ``python python/helmholtz_basics/__init__.py`` is safe but not
            especially useful; the imports above simply prepare the package
            namespace. When used from the shell (PowerShell, bash, etc.) without
            starting Python first, commands such as ``from .config import ...`` will
            raise syntax errors because the shell is not a Python interpreter.
            """
        ).strip()
    )

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
