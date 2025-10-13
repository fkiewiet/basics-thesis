# %% [markdown]
"""
Cells extracted from 04_visualisation.ipynb
"""

# %% [markdown]
"""
# 04 · Visualisation snippets

Collect helper functions and examples for visualising solutions and residuals.
"""

# %%
from helmholtz_basics import GridSpec, FiniteDifference, PointSource, gmres_solve
from helmholtz_basics.visualisation import plot_field, plot_residuals

%matplotlib inline

grid = GridSpec(dims=2, shape=(50, 50), lengths=(1.0, 1.0))
disc = FiniteDifference(wavenumber=25.0)
A = disc.build_operator(grid)
load = PointSource(location='centre')
b = load.build(grid)
result = gmres_solve(A, b, tol=1e-6, maxiter=200)

field = result.solution.reshape(grid.shape)
plot_field(grid, field, title='Field amplitude')
plot_residuals(result.residuals)

# %% [markdown]
"""
### Suggested follow-ups

- Swap `plot_field` for your own Matplotlib code if you need custom styling.
- Capture plots to disk inside a results directory for reproducible studies.
- Add widgets (e.g. ipywidgets) to build interactive dashboards.
"""
