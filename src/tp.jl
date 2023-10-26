# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

@kwdef struct TP{P,B,E,I} <: GeoStatsProcess
  params::P = nothing
  blur::B = nothing
  edge::E = nothing
  iter::I = 100
end
