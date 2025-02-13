# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

#------------------
# HOMOGENEOUS CASE
#------------------

function randsingle(rng::AbstractRNG, p::PoissonProcess{<:Real}, g)
  # simulate number of points
  λ = p.λ
  V = measure(g)
  n = rand(rng, Poisson(ustrip(λ * V)))

  # simulate n points
  iszero(n) ? nothing : PointSet(sample(rng, g, HomogeneousSampling(n)))
end

#--------------------
# INHOMOGENEOUS CASE
#--------------------

function randsingle(rng::AbstractRNG, p::PoissonProcess{<:Function}, g)
  # upper bound for intensity
  λmax = _maxintensity(p, g)

  # simulate a homogeneous process
  pset = randsingle(rng, PoissonProcess(ustrip(λmax)), g)

  # thin point pattern
  isnothing(pset) ? nothing : thin(pset, RandomThinning(x -> p.λ(x) / λmax))
end

function randsingle(rng::AbstractRNG, p::PoissonProcess{<:AbstractVector}, d::Domain)
  # simulate number of points
  λ = p.λ .* measure.(d)
  n = rand(rng, Poisson(ustrip(sum(λ))))

  # simulate point pattern
  iszero(n) ? nothing : PointSet(sample(rng, d, HomogeneousSampling(n, λ)))
end

function _maxintensity(p::PoissonProcess{<:Function}, g)
  points = sample(g, HomogeneousSampling(10000))
  λmin, λmax = extrema(p.λ, points)
  λmax + 0.05 * (λmax - λmin)
end
