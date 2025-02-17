# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    LindgrenProcess(range=1.0, sill=1.0; init=NearestInit())

Lindgren process with given `range` (correlation length)
and `sill` (total variance) as described in Lindgren 2011.

The process relies relies on a discretization of the Laplace-Beltrami 
operator on meshes and is adequate for highly curved domains (e.g. surfaces).

## References

* Lindgren et al. 2011. [An explicit link between Gaussian fields and
  Gaussian Markov random fields: the stochastic partial differential
  equation approach](https://rss.onlinelibrary.wiley.com/doi/10.1111/j.1467-9868.2011.00777.x)
"""
struct LindgrenProcess{ℒ<:Len,V} <: FieldProcess
  range::ℒ
  sill::V
  LindgrenProcess(range::ℒ, sill::V) where {ℒ<:Len,V} = new{float(ℒ),float(V)}(range, sill)
end

LindgrenProcess(range=1.0u"m", sill=1.0) = LindgrenProcess(addunit(range, u"m"), sill)
