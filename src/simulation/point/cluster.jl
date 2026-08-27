# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function randsingle(rng::AbstractRNG, p::ClusterProcess, g)
  # generate parents
  parents = rand(rng, p.proc, g)

  # initialize result
  result = eltype(parents)[]

  # append offsprings
  for point in parents
    pset = p.ofun(point)
    isnothing(pset) && continue
    pint = pset ∩ g
    isnothing(pint) && continue
    append!(result, _flat(pint))
  end

  isempty(result) ? nothing : PointSet(result)
end

_flat(p::Point) = [p]
_flat(m::Multi) = parent(m)
