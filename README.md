# basics-thesis

This repository houses exploratory work on Helmholtz solvers and GMRES
convergence studies.  Everything has now been refactored around a single
Python-first toolkit so experiments, scripts, and notebooks can share the same
building blocks.

- `python/helmholtz_basics/` contains the reusable pieces: grid builders,
  operator assembly, right-hand sides, GMRES wrappers, and a small experiment
  runner.  Import what you need and stay focused on the physics you want to
  poke at.
- [`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md) summarises how the package is
  organised and the kind of variations it already supports (dimension, grid
  size, frequency, loads, and discretisations).
- [`docs/PYTHON_COMPANION.md`](docs/PYTHON_COMPANION.md) is a short note that
  reads like a teammate explaining how to get started.

The original Julia notebook is still in the repo for reference, but new work is
expected to grow inside the Python modules.
