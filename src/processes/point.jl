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

#-----------------
# IMPLEMENTATIONS
#-----------------

include("point/binomial.jl")
include("point/poisson.jl")
include("point/inhibition.jl")
include("point/cluster.jl")
include("point/union.jl")
