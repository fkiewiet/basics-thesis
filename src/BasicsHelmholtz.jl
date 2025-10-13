module BasicsHelmholtz

using Base: @kwdef

"""
    using BasicsHelmholtz

High-level module that wires together mesh generation, discretisation, right-hand
side definitions, and linear solvers for Helmholtz experiments.  Users are expected
to construct a `GridSpec` describing the dimensionality and resolution, choose a
`Discretisation` implementation, and supply a `Load` before calling into the solver
layer.
"""

include("types.jl")
include("grids.jl")
include("operators.jl")
include("loads.jl")
include("solvers.jl")
include("experiments.jl")
include("visualisation.jl")

export GridSpec, CoordinateAxes, BoundarySpec, CartesianGrid, spacings,
       Discretisation, FivePointStencil, SevenPointStencil, NinePointStencil,
       Load, PointSource, PlaneWave, RandomLoad,
       assemble_operator, build_grid, build_load,
       solve_helmholtz, run_experiment!

end # module
