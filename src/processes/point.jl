# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    PointProcess <: GeoStatsProcess

Parent type of all point processes.
"""
abstract type PointProcess <: GeoStatsProcess end

"""
    ishomogeneous(process::PointProcess)

Tells whether or not the spatial point process `process` is homogeneous.
"""
ishomogeneous(process::PointProcess) = false

"""
    rand([rng], process::PointProcess, geometry, [nreals])
    rand([rng], process::PointProcess, domain, [nreals])

Generate one or `nreals` realizations of the point `process` inside
`geometry` or `domain`. Optionally specify the random number generator
`rng`.
"""
Base.rand(process::PointProcess, geomdom) = rand(Random.default_rng(), process, geomdom)

Base.rand(process::PointProcess, geomdom, nreals::Int) = rand(Random.default_rng(), process, geomdom, nreals)

Base.rand(rng::AbstractRNG, process::PointProcess, geomdom) = randsingle(rng, process, geomdom)

Base.rand(rng::AbstractRNG, process::PointProcess, geomdom, nreals::Int) = [randsingle(rng, process, geomdom) for _ in 1:nreals]

#-----------------
# IMPLEMENTATIONS
#-----------------

include("point/binomial.jl")
include("point/poisson.jl")
include("point/inhibition.jl")
include("point/cluster.jl")
include("point/union.jl")
