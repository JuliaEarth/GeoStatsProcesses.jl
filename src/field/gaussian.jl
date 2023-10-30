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

defaultmethod(::GaussianProcess, ::RandSetup) = LUMethod()
