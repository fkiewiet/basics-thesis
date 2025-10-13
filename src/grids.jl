"""
Construct coordinate axes for a Cartesian grid given a `GridSpec`.  Supports 1D,
2D, and 3D grids by truncating unused directions.
"""
function build_grid(spec::GridSpec)
    axes = ntuple(i -> range(0, spec.lengths[i]; length = spec.shape[i]), 3)
    CoordinateAxes(x = axes[1], y = axes[2], z = axes[3])
end

struct CartesianGrid
    axes::CoordinateAxes
end

"""
Convenience constructor that produces both the coordinate axes and a helper
`CartesianGrid` wrapper.
"""
CartesianGrid(spec::GridSpec) = CartesianGrid(build_grid(spec))

"""
Return the grid spacing for each axis; useful for defining discretisations.
"""
function spacings(spec::GridSpec)
    ntuple(i -> spec.lengths[i] / (spec.shape[i] - 1), 3)
end
