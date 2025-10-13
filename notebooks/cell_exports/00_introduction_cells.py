# %% [markdown]
"""
Cells extracted from 00_introduction.ipynb
"""

# %% [markdown]
"""
# Helmholtz basics notebook tour

This short notebook introduces the Python toolkit that drives the experiments.
Each subsequent notebook focuses on one idea so you can mix and match pieces.
"""

# %%
from helmholtz_basics import GridSpec, FiniteDifference, PointSource, gmres_solve, SweepConfig, run_experiment
print('Toolkit pieces imported successfully.')

# %% [markdown]
"""
## Suggested reading order

1. **01_grid_and_operator** – build grids and assemble the Helmholtz matrix.
2. **02_single_run** – run GMRES for a single configuration and plot diagnostics.
3. **03_parameter_sweeps** – explore parameter variations with the experiment driver.
4. **04_visualisation** – collect helper snippets for plotting residuals and fields.

The legacy Julia walkthrough now lives in `notebooks/legacy/` for reference.
"""
