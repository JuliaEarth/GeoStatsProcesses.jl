# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    PoissonProcess(λ)

A Poisson process with intensity `λ`. For a homogeneous process,
define `λ` as a constant real value, while for an inhomogeneous process,
define `λ` as a function or vector of values. If `λ` is a vector, it is
assumed that the process is associated with a `Domain` with the same
number of elements as `λ`.
"""
struct PoissonProcess{L<:Union{Real,Function,AbstractVector}} <: PointProcess
  λ::L
end

ishomogeneous(p::PoissonProcess{<:Real}) = true

ishomogeneous(p::PoissonProcess{<:Function}) = false

ishomogeneous(p::PoissonProcess{<:AbstractVector}) = false
