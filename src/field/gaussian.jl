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
  N = 10_000
  domain, _ = unview(setup.domain)
  if domain isa Grid
    FFTMethod()
  elseif nelements(domain) < N
    LUMethod()
  else
    SEQMethod()
  end
end
