# %% [markdown]
"""
Cells extracted from 01_grid_and_operator.ipynb
"""

# %% [markdown]
"""
# 01 · Grid and operator assembly

Set up grid specifications and assemble discrete Helmholtz operators in 1D, 2D, or 3D.
"""

# %%
from helmholtz_basics import GridSpec, FiniteDifference

# 2D grid example
grid = GridSpec(dims=2, shape=(50, 50), lengths=(1.0, 1.0))
print(grid)

# Assemble the standard second-order finite-difference operator
discretisation = FiniteDifference(wavenumber=25.0)
A = discretisation.build_operator(grid)
print(A.shape, 'non-zeros:', A.nnz)

# %% [markdown]
"""
### What to try next

- Change `dims`, `shape`, or `lengths` to see how the grid reacts.
- Inspect `A` (e.g. `A.diagonal()[:5]`) to understand the stencil weights.
- Swap in your own discretisation class if you want higher-order stencils.
"""
