# basics-thesis

This repository houses exploratory work on Helmholtz solvers and GMRES
convergence studies.  The `src/` directory now contains a light-weight module
structure that can be imported from notebooks or scripts to organise solver
components.

If you prefer working in Python, check out the companion scaffolding in
`python/helmholtz_basics/`, which mirrors the Julia layout.  A short tour is
available in [`docs/PYTHON_COMPANION.md`](docs/PYTHON_COMPANION.md).

See [`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md) for the proposed modular
layout and migration steps on the Julia side.
