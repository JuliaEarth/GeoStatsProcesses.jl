# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    GeoStatsProcess

A geostatistical process that can be simulated with the `rand` function.
"""
abstract type GeoStatsProcess end

"""
    RandMethod

A `rand` method for geostatistical processes.
"""
abstract type RandMethod end

"""
    DefaultRandMethod()

Default `rand` method used by some geostatistical processes.
"""
struct DefaultRandMethod <: RandMethod end

"""
    defaultmethod(process, setup) -> RandMethod

Returns the default method for the `process` and `setup`.
"""
defaultmethod(process, setup) = DefaultRandMethod()

# ----------------
# IMPLEMENTATIONS
# ----------------

include("processes/point.jl")
include("processes/field.jl")

# -----------
# IO METHODS
# -----------

Base.summary(io::IO, process::GeoStatsProcess) = print(io, prettyname(process))

function Base.show(io::IO, process::GeoStatsProcess)
  name = prettyname(process)
  ioctx = IOContext(io, :compact => true)
  print(io, "$name(")
  printfields(ioctx, process, singleline=true)
  print(io, ")")
end

function Base.show(io::IO, ::MIME"text/plain", process::GeoStatsProcess)
  summary(io, process)
  printfields(io, process)
end
