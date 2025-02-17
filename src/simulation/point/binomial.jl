# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

randsingle(rng::AbstractRNG, p::BinomialProcess, g) = PointSet(sample(rng, g, HomogeneousSampling(p.n)))
