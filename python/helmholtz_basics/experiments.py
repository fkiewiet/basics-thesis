"""Experiment runners that sweep parameters."""

from __future__ import annotations

from typing import Dict, Iterable, Optional

from .config import SweepConfig
from .loads import build_load
from .operators import assemble_operator
from .solvers import GMRESOptions, gmres_solve


def run_experiment(
    config: SweepConfig,
    *,
    solver_options: Optional[GMRESOptions] = None,
) -> Iterable[Dict[str, object]]:
    """Yield a record for each combination of parameters."""

    for case in config.iter_configs():
        matrix = assemble_operator(case.discretisation, case.grid, case.wavenumber)
        rhs = build_load(case.load, case.grid)
        result = gmres_solve(matrix, rhs, options=solver_options)
        record: Dict[str, object] = {
            "grid": case.grid,
            "wavenumber": case.wavenumber,
            "load": case.load,
            "discretisation": case.discretisation,
            "solver_result": result,
        }
        if case.label is not None:
            record["label"] = case.label
        if case.metadata:
            record["metadata"] = dict(case.metadata)
        if config.hook is not None:
            config.hook(case, result)
        yield record
