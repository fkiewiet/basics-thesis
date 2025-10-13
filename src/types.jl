"""
    struct GridSpec
Describes the dimensionality, coordinate extents, and resolution of the
computational domain.  This object is intentionally lightweight so it can be
serialised from notebooks or configuration files.
"""
Base.@kwdef mutable struct GridSpec
    dims::Int
    lengths::NTuple{3,Float64}
    shape::NTuple{3,Int}
    boundary::Symbol = :dirichlet
end

"""
    struct CoordinateAxes
Convenience container for the grid coordinates in each spatial direction.
"""
Base.@kwdef struct CoordinateAxes
    x::AbstractVector
    y::AbstractVector
    z::AbstractVector
end

"""
Abstract discretisation strategy, e.g. finite differences or finite elements.
"""
abstract type Discretisation end

"""
Abstract type for the forcing term in the Helmholtz equation.
"""
abstract type Load end

"""
    struct BoundarySpec
Encodes the boundary conditions applied to the Helmholtz problem.  Extend this
with additional fields (e.g. impedance parameters) as you explore new scenarios.
"""
Base.@kwdef mutable struct BoundarySpec
    kind::Symbol = :dirichlet
    value::Float64 = 0.0
end
