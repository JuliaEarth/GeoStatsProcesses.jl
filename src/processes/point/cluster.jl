# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    ClusterProcess(proc, ofun)

A cluster process with parent process `proc` and offsprings generated
with `ofun`. It is a function that takes a parent point and returns
a point pattern from another point process.

    ClusterProcess(proc, offs, gfun)

Alternatively, specify the parent process `proc`, the offspring process
`offs` and the geometry function `gfun`. It is a function that takes a
parent point and returns a geometry or domain for the offspring process.
"""
struct ClusterProcess{P<:PointProcess,F<:Function} <: PointProcess
  proc::P
  ofun::F
end

ClusterProcess(p::PointProcess, o::PointProcess, gfun::Function) = ClusterProcess(p, parent -> rand(o, gfun(parent)))
