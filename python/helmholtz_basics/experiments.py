"""Experiment runners that sweep parameters."""

from __future__ import annotations

from typing import Dict, Iterable

from .config import SweepConfig
from .loads import build_load
from .operators import assemble_operator
from .solvers import gmres_solve


def run_experiment(config: SweepConfig) -> Iterable[Dict[str, object]]:
    """Yield a record for each combination of parameters."""

    for case in config.iter_configs():
        matrix = assemble_operator(case.discretisation, case.grid, case.wavenumber)
        rhs = build_load(case.load, case.grid)
        result = gmres_solve(matrix, rhs)
        record = {
            "grid": case.grid,
            "wavenumber": case.wavenumber,
            "load": case.load,
            "discretisation": case.discretisation,
            "solver_result": result,
        }
        if config.hook is not None:
            config.hook(case, result)
        yield record
