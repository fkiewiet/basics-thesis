using SparseArrays
using LinearAlgebra: I

"""
    struct FiniteDifferenceStencil
Encapsulates a dimension-aware finite-difference stencil.  Use concrete subtypes
for common patterns, or define new ones for alternative discretisations.
"""
abstract type FiniteDifferenceStencil <: Discretisation end

struct FivePointStencil <: FiniteDifferenceStencil end
struct SevenPointStencil <: FiniteDifferenceStencil end
struct NinePointStencil <: FiniteDifferenceStencil end

"""
Construct a sparse Helmholtz operator for the given grid and discretisation.
This function only sketches the API; fill in the implementation as you migrate
code out of the notebook.
"""
function assemble_operator(spec::GridSpec, disc::FiniteDifferenceStencil; ω::Float64)
    n = prod(spec.shape[1:spec.dims])
    # placeholder identity; replace with real stencil assembly
    SparseArrays.spdiagm(0 => fill(ω^2, n)) + I
end

"""
Hook to support alternative discretisation families (e.g. FEM, DG).  Overload for
your own `Discretisation` subtypes.
"""
assemble_operator(spec::GridSpec, disc::Discretisation; kwargs...) = error("Discretisation not implemented: $(typeof(disc))")
