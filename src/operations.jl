# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    ThinningMethod

A method for thinning spatial point processes and patterns.
"""
abstract type ThinningMethod end

"""
    thin(process, method)

Thin spatial point `process` with thinning `method`.
"""
function thin end

#-----------------
# IMPLEMENTATIONS
#-----------------

include("thinning/random.jl")

#---------------
# PROCESS UNION
#---------------

"""
    p₁ ∪ p₂

Return the union of point processes `p₁` and `p₂`.
"""
Base.union(p₁::PointProcess, p₂::PointProcess) = UnionProcess(p₁, p₂)

Base.union(p₁::PoissonProcess{<:Real}, p₂::PoissonProcess{<:Real}) = PoissonProcess(p₁.λ + p₂.λ)

Base.union(p₁::PoissonProcess{<:Function}, p₂::PoissonProcess{<:Function}) = PoissonProcess(x -> p₁.λ(x) + p₂.λ(x))

Base.union(p₁::PoissonProcess{<:Real}, p₂::PoissonProcess{<:Function}) = PoissonProcess(x -> p₁.λ + p₂.λ(x))

Base.union(p₁::PoissonProcess{<:Function}, p₂::PoissonProcess{<:Real}) = PoissonProcess(x -> p₁.λ(x) + p₂.λ)

Base.union(p₁::PoissonProcess{<:AbstractVector}, p₂::PoissonProcess{<:AbstractVector}) =
  PoissonProcess(x -> p₁.λ + p₂.λ)
