# Python-first layout for Helmholtz + GMRES experiments

The repository used to experiment with a Julia module split, but the active
codebase is now entirely Python.  Everything lives under
`python/helmholtz_basics/`, and the goal is to give notebooks and command-line
scripts a consistent set of helpers for exploring different spatial dimensions,
grid sizes, wavenumbers, loads, and discretisations.  The Jupyter side mirrors
that modularity: the `notebooks/` directory contains five short notebooks that
each focus on one step of the workflow (see [`docs/NOTEBOOKS.md`](NOTEBOOKS.md)
for a quick index).

## Package map

```
python/helmholtz_basics/
  __init__.py        # convenience exports for interactive use
  config.py          # dataclasses describing grids, problems, and sweeps
  grid.py            # mesh/coordinate generation utilities
  loads.py           # right-hand side factories (point, plane wave, random)
  operators.py       # finite-difference Helmholtz assembly in 1/2/3D
  solvers.py         # GMRES wrapper + solver options
  experiments.py     # small driver that loops through sweep combinations
  visualisation.py   # optional Matplotlib helpers for results
```

You can import the package in a notebook or script and reach for the pieces you
need:

```python
from helmholtz_basics import (
    GridSpec,
    FiniteDifference,
    PointSource,
    gmres_solve,
)
```

The public API prefers plain dataclasses and protocols so it is straightforward
to add your own load type or discretisation without editing existing files.

## Handling different dimensions and grids

`GridSpec` records both the number of dimensions and the shape/lengths per axis.
It provides convenience properties such as `spacing` and `size`, so assembly
routines only need to look at `spec.dims` entries.  The same `FiniteDifference`
class handles 1D, 2D, and 3D because helper functions iterate over the axes.

## Loads and discretisations

Loads implement a light `Load` protocol with a `build(grid)` method.  The built
ins cover point sources, phase-controlled plane waves, and random draws.  New
sources only need to honour the protocol.  Discretisations do the same via a
`Discretisation` protocol—`FiniteDifference` ships with a standard second-order
stencil, and alternate schemes (higher order, absorbing layers, FEM wrappers)
can be dropped in beside it.

## Sweep orchestration

`SweepConfig` holds the Cartesian product of parameters you want to explore and
an optional callback hook for logging.  `run_experiment` assembles the operator,
builds the right-hand side, calls the solver, and yields a record per case.  The
records are plain dictionaries so downstream code can serialise them, compute
statistics, or plot without ceremony.

## Suggested workflow

1. Start in a notebook or script and `from helmholtz_basics import *` the pieces
   you need.
2. Create a `GridSpec`, choose a `FiniteDifference` discretisation, and pick one
   of the load factories (or implement your own).
3. Run individual solves with `gmres_solve` while prototyping, then switch to a
   `SweepConfig` and `run_experiment` once you want to study trends.
4. Use the optional plotting helpers in `visualisation.py` to keep quick-look
   figures separate from the solver core.
