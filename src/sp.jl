# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

@kwdef struct SP{E,S,ST,N,FB,FT} <: GeoStatsProcess
  environment::E
  state::S = nothing
  stack::ST = :erosional
  nepochs::N = 10
  fillbase::FB = NaN
  filltop::FT = NaN
end
