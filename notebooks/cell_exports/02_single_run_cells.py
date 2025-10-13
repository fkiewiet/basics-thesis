# %% [markdown]
"""
Cells extracted from 02_single_run.ipynb
"""

# %% [markdown]
"""
# 02 · Single GMRES run

Drive the GMRES solver for a single configuration and track residuals.
"""

# %%
from helmholtz_basics import GridSpec, FiniteDifference, PointSource, gmres_solve

# Build grid and operator
grid = GridSpec(dims=2, shape=(60, 60), lengths=(1.0, 1.0))
disc = FiniteDifference(wavenumber=30.0)
A = disc.build_operator(grid)

# Point-source load at the centre
load = PointSource(location='centre', amplitude=1.0)
b = load.build(grid)

# Run GMRES
result = gmres_solve(A, b, tol=1e-6, maxiter=200)
print(f"Converged: {result.converged} in {result.iterations} iterations")

# %% [markdown]
"""
### Suggested follow-ups

- Plot `result.residuals` to inspect convergence.
- Experiment with different loads (e.g. `PlaneWave`, `RandomLoad`).
- Adjust `wavenumber` or `grid` to see how the solver reacts.
"""
