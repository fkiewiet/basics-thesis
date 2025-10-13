# Python companion outline

I sketched a small Python package in `python/helmholtz_basics/` so the Julia
notebook has a twin that feels familiar.  The folder mirrors the Julia modules:
configuration dataclasses live in `config.py`, grids in `grid.py`, operators in
`operators.py`, and so on.  Everything is wired together through the package's
`__init__` so you can write imports such as:

```python
from helmholtz_basics import GridSpec, FiniteDifference, PointSource, gmres_solve
```

The stubs already cover the things we talked about tweaking—dimensions, grid
sizes, wavenumbers, and different loads.  Right now the defaults use a simple
second-order finite difference scheme, but you can add new discretisations by
creating another class with an `assemble` method.  Experiments are handled by a
small sweep helper that loops over every combination of parameters and hands
results to an optional callback for logging or plotting.

Once SciPy and Matplotlib are available, the code is ready for real runs.  The
plan is to port pieces from the notebook into these modules one at a time so we
keep the story in the write-up clear while the heavy lifting happens in tidy,
testable Python files.
