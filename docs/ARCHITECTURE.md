# Modularising the Helmholtz GMRES project

This repository currently contains the exploratory notebook
`Helmholtz_GMRES_Julia_MSc.ipynb`.  To support experiments with different
spatial dimensions, grid sizes, frequencies, right-hand sides, and
alternative discretisations, treat the notebook as a consumer of a small
library housed in `src/`.

## Proposed layout

```
src/
  BasicsHelmholtz.jl     # umbrella module re-exporting public API
  types.jl               # declarative data containers (grid, boundary, load)
  grids.jl               # coordinate generation utilities
  operators.jl           # discretisation definitions and matrix assembly
  loads.jl               # right-hand side factories
  solvers.jl             # GMRES and related Krylov solvers
  experiments.jl         # orchestration helpers for parameter sweeps
  visualisation.jl       # plotting hooks kept optional
```

Notebooks and scripts should `include("src/BasicsHelmholtz.jl")` and then call
high-level helpers such as `solve_helmholtz` or `run_experiment!`.  All heavy
lifting happens inside the library modules so that functionality is shareable
across multiple front-ends.

## Extending by dimension

The `GridSpec` object stores dimension, lengths, and discrete shape in a single
place.  This allows you to create configurations such as:

```julia
GridSpec(dims = 1, lengths = (1.0, 0.0, 0.0), shape = (512, 1, 1))
GridSpec(dims = 2, lengths = (1.0, 1.0, 0.0), shape = (256, 256, 1))
GridSpec(dims = 3, lengths = (1.0, 1.0, 1.0), shape = (64, 64, 64))
```

Grid-dependent logic should use the first `spec.dims` entries in each tuple so
functions naturally degrade from 3D to 2D and 1D without rewriting code.

## Managing discretisations

Implement each spatial discretisation as a subtype of `Discretisation`.  For
finite differences, use small concrete structs such as `FivePointStencil` and
specialise `assemble_operator`.  Alternative discretisations—finite elements,
discontinuous Galerkin, or PML-augmented operators—can register their own
subtypes and overloads without modifying existing files.

## Right-hand sides

Similarly, declare each forcing term as a subtype of `Load`.  The `build_load`
dispatch allows custom behaviour per load (e.g. mapping a physical point source
to the closest grid node, drawing random fields, or injecting phase-controlled
plane waves).  Keeping loads distinct makes it easy to compose experiments that
sweep through source types.

## Experiment orchestration

The `ExperimentConfig` container captures the cross-product of parameters.  A
single call to `run_experiment!` can sweep through grids, discretisations,
frequencies, and loads, handing the results to a callback for logging, plotting,
or saving to disk.  This isolates scientific questions from plumbing code and
encourages reproducibility.

## Next steps

1. Port the existing notebook functions into the appropriate modules.
2. Replace placeholder implementations (identity operator, zero solution) with
   working finite-difference assembly and GMRES logic.
3. Write lightweight tests that exercise the modules in `src/` so future
   refactors remain safe.
4. Gradually move plotting helpers into `visualisation.jl`, allowing the
   notebook to focus on narrative and analysis rather than implementation.
