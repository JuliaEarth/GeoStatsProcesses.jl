# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

const Len{T} = Quantity{T,u"𝐋"}

addunit(x::Number, u) = x * u
addunit(x::Quantity, _) = x

function _pairwise(f::GeoStatsFunction, dom₁, dom₂)
  s = ustrip(sill(f))
  F = GeoStatsFunctions.pairwise(f, dom₁, dom₂)
  isbanded(f) || (F .= s .- F)
  F
end

function _pairwise(f::GeoStatsFunction, dom)
  s = ustrip(sill(f))
  F = GeoStatsFunctions.pairwise(f, dom)
  isbanded(f) || (F .= s .- F)
  F
end
