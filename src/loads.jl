using Random

struct PointSource <: Load
    location::NTuple{3,Float64}
    amplitude::ComplexF64
end

struct PlaneWave <: Load
    direction::NTuple{3,Float64}
    amplitude::ComplexF64
    wavenumber::Float64
end

struct RandomLoad <: Load
    seed::Int
    amplitude::Float64
end

"""
Return a discretised load vector matching the grid and discretisation choices.
Each load type can specialise this function.
"""
function build_load(spec::GridSpec, load::PointSource)
    n = prod(spec.shape[1:spec.dims])
    rhs = zeros(ComplexF64, n)
    rhs[1] = load.amplitude
    rhs
end

function build_load(spec::GridSpec, load::PlaneWave)
    n = prod(spec.shape[1:spec.dims])
    fill(load.amplitude, n)
end

function build_load(spec::GridSpec, load::RandomLoad)
    rng = MersenneTwister(load.seed)
    n = prod(spec.shape[1:spec.dims])
    load.amplitude .* randn(rng, n)
end

build_load(::GridSpec, load::Load) = error("Load not implemented: $(typeof(load))")
