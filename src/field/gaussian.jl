# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

@kwdef struct GaussianProcess{V,T} <: FieldProcess
  variogram::V = GaussianVariogram()
  mean::T = 0.0
end

#---------
# METHODS
#---------

include("gaussian/lu.jl")
include("gaussian/fft.jl")
include("gaussian/seq.jl")

function defaultmethod(::GaussianProcess, setup::RandSetup)
  dom = setup.domain
  udom, _ = unview(dom)
  if udom isa Grid
    FFTMethod()
  elseif nelements(dom) < 10000
    LUMethod()
  else
    SEQMethod()
  end
end
