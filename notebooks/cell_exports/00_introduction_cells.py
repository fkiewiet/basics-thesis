# %% [markdown]
"""
Cells extracted from 00_introduction.ipynb
"""

# %% [markdown]
"""
# Helmholtz basics tour

Welcome! This short intro gives you the landmarks for the Python toolkit so you can jump straight to the bit you need. Each follow-up notebook zooms in on one task—building grids, running a single solve, sweeping parameters, or plotting diagnostics.
"""

# %% [markdown]
"""
## How to use this tour

* Run the quick import cell below to make sure the core pieces load.
* Skim the map to decide which focused notebook you want to open next.
* Jump to the docs whenever you want a deeper dive or a script-friendly view.
"""

# %%
# A tiny smoke test so you know the toolkit is wired up.
from helmholtz_basics.config import GridSpec
from helmholtz_basics.loads import PointSource

grid = GridSpec(dims=2, shape=(32, 32), lengths=(1.0, 1.0))
source = PointSource(location=(0.5, 0.5))
print(f"Grid spacing: {grid.spacing}")
print("Sample load entry:", source.build(grid).nonzero()[0][:5])

try:
    from helmholtz_basics.operators import FiniteDifference
except ImportError as exc:  # noqa: F841
    print("FiniteDifference needs SciPy. Install it when you want to assemble matrices.")
else:
    disc = FiniteDifference()
    print(f"Discretisation ready: {disc.name}")

# %% [markdown]
"""
## Map of the follow-up notebooks

1. **01_grid_and_operator** – set up 1D/2D/3D grids and assemble Helmholtz matrices.
2. **02_single_run** – solve one configuration, inspect GMRES residuals, and sanity-check the field.
3. **03_parameter_sweeps** – script repeatable runs across dimensions, wavenumbers, or loads.
4. **04_visualisation** – keep plotting snippets handy for residual histories and field slices.
"""

# %% [markdown]
"""
## Handy extras

* **Docs:** `docs/ARCHITECTURE.md` sketches how the package fits together; `docs/PYTHON_COMPANION.md` adds plain-language notes.
* **Script exports:** every notebook has a twin under `notebooks/cell_exports/` with `# %%` markers for quick reuse.
* **Legacy material:** the original long-form walkthrough now lives in `notebooks/legacy/` if you ever want the full story.
"""
