# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

# ------
# UNITS
# ------

const Len{T} = Quantity{T,u"𝐋"}

addunit(x::Number, u) = x * u
addunit(x::Quantity, _) = x

prettyname(obj) = prettyname(typeof(obj))
function prettyname(T::Type)
  name = string(T)
  name = replace(name, r"{.*" => "")
  replace(name, r".*\." => "")
end

# ---
# IO
# ---

printfields(io, obj; kwargs...) = printfields(io, obj, fieldnames(typeof(obj)); kwargs...)
function printfields(io, obj, fnames; singleline=false)
  if singleline
    vals = map(enumerate(fnames)) do (i, field)
      val = getfield(obj, i)
      str = repr(val, context=io)
      "$field: $str"
    end
    join(io, vals, ", ")
  else
    len = length(fnames)
    for (i, field) in enumerate(fnames)
      div = i == len ? "\n└─ " : "\n├─ "
      val = getfield(obj, i)
      str = repr(val, context=io)
      print(io, "$div$field: $str")
    end
  end
end

# -----------
# SIMULATION
# -----------

function _pairwise(fun, dom₁, dom₂)
  s = ustrip(sill(fun))
  F = GeoStatsFunctions.pairwise(fun, dom₁, dom₂)
  isbanded(fun) || (F .= s .- F)
  F
end

function _pairwise(fun, dom)
  s = ustrip(sill(fun))
  F = GeoStatsFunctions.pairwise(fun, dom)
  isbanded(fun) || (F .= s .- F)
  F
end

function _scale(dom, dat, fun)
  α₁ = _scalefactor(dom)
  α₂ = isnothing(dat) ? 1 : _scalefactor(domain(dat))
  α₃ = _scalefactor(fun)
  α = inv(max(α₁, α₂, α₃))

  sdom = dom |> Scale(α)
  sdat = isnothing(dat) ? nothing : (dat |> Scale(α))
  sfun = GeoStatsFunctions.scale(fun, α)

  sdom, sdat, sfun, α
end

function _scalefactor(domain::Domain)
  pmin, pmax = extrema(boundingbox(domain))
  cmin = abs.(to(pmin))
  cmax = abs.(to(pmax))
  ustrip(max(cmin..., cmax...))
end

_scalefactor(fun::GeoStatsFunction) = ustrip(range(fun))
