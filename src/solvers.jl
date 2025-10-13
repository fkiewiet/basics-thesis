using LinearAlgebra

"""
Basic GMRES solver interface that accepts a pre-assembled operator.  Extend this
with preconditioner options or alternative Krylov solvers.
"""
function solve_helmholtz(A::AbstractMatrix, b; restart::Union{Nothing,Int}=nothing, tol=1e-8, maxiter=size(A,1))
    # Placeholder returning zero solution.  Replace with GMRES implementation.
    zeros(eltype(b), length(b)), fill(0.0, 0)
end

"""
High-level helper that assembles the operator and load from declarative inputs.
"""
function solve_helmholtz(spec::GridSpec, disc::Discretisation, load::Load; ω, solver_kwargs...)
    A = assemble_operator(spec, disc; ω)
    b = build_load(spec, load)
    solve_helmholtz(A, b; solver_kwargs...)
end
