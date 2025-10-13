mutable struct ExperimentConfig
    specs::Vector{GridSpec}
    discretisations::Vector{Discretisation}
    loads::Vector{Load}
    frequencies::Vector{Float64}
end

"""
Iterate over combinations of grid specs, discretisations, loads, and frequencies
running the solver for each case.  Records metadata so notebooks can focus on
analysis and visualisation.
"""
function run_experiment!(callback, config::ExperimentConfig; solver_kwargs...)
    for spec in config.specs, disc in config.discretisations, load in config.loads, ω in config.frequencies
        solution, history = solve_helmholtz(spec, disc, load; ω, solver_kwargs...)
        callback((; spec, disc, load, ω, solution, history))
    end
end
