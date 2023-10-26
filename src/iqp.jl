# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

@kwdef struct IQP{TR,TS,O,P,IN,S,T,I} <: FieldProcess
  trainimg::TR
  tilesize::TS
  overlap::O = nothing
  path::P = :raster
  inactive::IN = nothing
  soft::S = nothing
  tol::T = 0.1
  init::I = NearestInit()
end
