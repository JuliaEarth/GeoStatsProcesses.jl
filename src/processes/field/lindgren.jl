# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    LindgrenProcess(range=1.0, sill=1.0)

Lindgren process with given `range` and `sill`.

The process relies on the discretization of the Laplace-Beltrami
operator on meshes and is associated with a specific SPDE.

## Examples

```julia
# set range
LindgrenProcess(20.0)

# set range and sill
LindgrenProcess(20.0, 2.0)
```

## References

* Lindgren et al. 2011. [An explicit link between Gaussian fields and
  Gaussian Markov random fields: the stochastic partial differential
  equation approach](https://rss.onlinelibrary.wiley.com/doi/10.1111/j.1467-9868.2011.00777.x)

### Notes

The process is particularly useful in highly curved domains (e.g. surfaces)
given that it approximates geodesics as opposed to naive Euclidean distances.

It is also known as Gaussian Markov Random Field (GMRF) in the literature.

"""
struct LindgrenProcess{ℒ<:Len,V} <: FieldProcess
  range::ℒ
  sill::V
  function LindgrenProcess(range::ℒ, sill::V) where {ℒ<:Len,V}
    @assert range > zero(range) "range must be positive"
    @assert sill > zero(sill) "sill must be positive"
    new{float(ℒ),float(V)}(range, sill)
  end
end

LindgrenProcess(range=1.0u"m", sill=1.0) = LindgrenProcess(aslen(range), sill)

iscontinuous(process::LindgrenProcess) = true

isanalytical(process::LindgrenProcess) = true

"""
    precisionmatrix(process::LindgrenProcess, domain::Domain)

Compute the Matérn precision matrix of the Lindgren process at
the vertices of the given `domain`.
"""
function precisionmatrix(process::LindgrenProcess, domain)
  # process parameters
  𝓁 = process.range
  σ² = process.sill

  # Laplace-Beltrami operator
  W = laplacematrix(domain)
  M = measurematrix(domain)
  Δ = inv(M) * W

  # retrieve parametric dimension
  d = paramdim(domain)

  # LHS of SPDE (κ² - Δ)Z = τW with Δ = M⁻¹W
  α = 2
  ν = α - d / 2
  κ = 1 / 𝓁
  A = κ^2 * I - Δ

  # Matérn precision matrix
  τ² = σ² * κ^(2ν) * (4π)^(d / 2) * gamma(α) / gamma(ν)
  ustrip.(A'A / τ²)
end
