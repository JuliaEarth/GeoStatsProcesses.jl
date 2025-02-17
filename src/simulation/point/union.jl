# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function randsingle(rng::AbstractRNG, p::UnionProcess, g)
  pp₁ = rand(rng, p.p₁, g)
  pp₂ = rand(rng, p.p₂, g)
  PointSet([collect(pp₁); collect(pp₂)])
end
