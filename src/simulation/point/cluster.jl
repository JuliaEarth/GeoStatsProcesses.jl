# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function randsingle(rng::AbstractRNG, p::ClusterProcess, g)
  # generate parents
  parents = rand(rng, p.proc, g)

  # generate offsprings
  offsprings = filter(!isnothing, p.ofun.(parents))

  # intersect with geometry
  intersects = filter(!isnothing, [o ∩ g for o in offsprings])

  # combine offsprings into single set
  isempty(intersects) ? nothing : PointSet(mapreduce(collect, vcat, intersects))
end
