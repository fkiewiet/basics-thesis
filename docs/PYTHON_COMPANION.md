# Python toolkit in plain words

The working code now lives entirely in `python/helmholtz_basics/`.  Think of it
as a tray of Lego bricks: grids, operators, loads, solvers, and a tiny
experiment runner.  You only grab the bricks you need for the question at hand.

The package exports the usual suspects through `__init__.py`, so quick-start
imports stay short.  Open a Python interpreter (``python`` on macOS/Linux,
``py`` on Windows) _before_ running the import; shells like PowerShell will throw
syntax errors if you paste Python directly into them.

```python
from helmholtz_basics import GridSpec, FiniteDifference, PointSource, gmres_solve
```

- `config.py` keeps the small dataclasses that describe grids and parameter
  sweeps.
- `grid.py` turns those specs into NumPy coordinate arrays when you want to
  inspect or plot a field.
- `operators.py` currently offers a second-order finite-difference
  discretisation that works in 1D, 2D, or 3D.  Drop a new class in the same file
  if you fancy higher order stencils or absorbing layers.
- `loads.py` holds the right-hand side factories—point sources, plane waves, and
  a quick random generator for stress tests.
- `solvers.py` wraps SciPy's GMRES so you always get back the solution together
  with the residual history.
- `experiments.py` wires the pieces into a simple sweep loop and optionally lets
  you log or plot each run via a callback.

Once SciPy and Matplotlib are installed, you can run small scripts straight away
or keep working from notebooks.  The aim is to keep the notebook narrative light
and move the heavy lifting into these modules where it is easier to maintain and
extend.

Prefer to tinker in Jupyter?  Head to the `notebooks/` folder.  It now contains
an introduction plus four short companions that mirror the steps above: grid
building, single solves, parameter sweeps, and plotting.  Copy a cell, tweak a
parameter, and keep iterating without scrolling through a monolithic notebook.
