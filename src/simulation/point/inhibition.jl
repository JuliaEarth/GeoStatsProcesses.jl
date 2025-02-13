# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

randsingle(rng::AbstractRNG, p::InhibitionProcess, g) = PointSet(sample(rng, g, MinDistanceSampling(p.δ)))
