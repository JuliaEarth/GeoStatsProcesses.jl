# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

# ------
# UNITS
# ------

const Len{T} = Quantity{T,u"𝐋"}

addunit(x::Number, u) = x * u
addunit(x::Quantity, _) = x

# ---
# IO
# ---

prettyname(obj) = prettyname(typeof(obj))
function prettyname(T::Type)
  name = string(T)
  name = replace(name, r"{.*" => "")
  replace(name, r".*\." => "")
end

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

# ---------
# PAIRWISE
# ---------

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
