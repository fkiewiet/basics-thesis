# %% [markdown]
"""
Cells extracted from 03_parameter_sweeps.ipynb
"""

# %% [markdown]
"""
# 03 · Parameter sweeps

Use the experiment runner to explore multiple configurations automatically.
"""

# %%
from helmholtz_basics import GridSpec, FiniteDifference, PointSource, PlaneWave, SweepConfig, run_experiment

config = SweepConfig(
    grids=[
        GridSpec(dims=2, shape=(40, 40), lengths=(1.0, 1.0)),
        GridSpec(dims=2, shape=(80, 80), lengths=(1.0, 1.0)),
    ],
    discretisations=[
        FiniteDifference(wavenumber=20.0),
        FiniteDifference(wavenumber=35.0),
    ],
    loads=[
        PointSource(location='centre'),
        PlaneWave(direction=(1.0, 0.0), phase=0.0),
    ],
)

for record in run_experiment(config):
    summary = (
        record['grid'].shape,
        record['discretisation'].wavenumber,
        type(record['load']).__name__,
        record['solver'].converged,
        record['solver'].iterations,
    )
    print(summary)

# %% [markdown]
"""
### Suggested follow-ups

- Extend `SweepConfig` with more grids, loads, or solvers.
- Attach a callback to log or save results as they stream.
- Convert the loop into a Pandas DataFrame for analysis.
"""
