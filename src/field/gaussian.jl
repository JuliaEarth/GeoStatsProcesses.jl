# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    GaussianProcess([parameters])

Gaussian process with given variogram and mean.

## Parameters

* `variogram` - Theoretical variogram model (default to `GaussianVariogram()`)
* `mean`      - Mean field value (default to `0.0`)

"""
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
