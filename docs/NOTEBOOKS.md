# Notebook overview

The original notebook grew into a long, monolithic walkthrough.  The new layout
splits the story into focused pieces so you can jump directly to the stage you
care about:

- `00_introduction.ipynb` gives a short tour and points to the other notebooks.
- `01_grid_and_operator.ipynb` builds grids and assembles Helmholtz operators.
- `02_single_run.ipynb` performs one GMRES solve from start to finish.
- `03_parameter_sweeps.ipynb` demonstrates the sweep driver for parameter
  studies.
- `04_visualisation.ipynb` collects plotting snippets for fields and residuals.
- `legacy/Helmholtz_GMRES_Julia_MSc.ipynb` preserves the original narrative for
  reference.

Each notebook leans on the shared `helmholtz_basics` Python package so code
snippets stay short and reusable.  Feel free to copy the cells into your own
experiments or add new notebooks beside them for specialised studies (e.g.
preconditioning trials, 3D benchmarks, or load comparisons).
