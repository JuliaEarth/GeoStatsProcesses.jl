# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    GaussianProcess(variogram, mean)
    GaussianProcess(variogram; mean)
    GaussianProcess(; variogram, mean)

Gaussian process with given variogram and mean.

## Parameters

* `variogram` - Theoretical variogram model (default to `GaussianVariogram()`)
* `mean`      - Mean field value (default to `0.0`)

"""
struct GaussianProcess{V,T} <: FieldProcess
  variogram::V
  mean::T
end

GaussianProcess(variogram; mean=0.0) = GaussianProcess(variogram, mean)
GaussianProcess(; variogram=GaussianVariogram(), mean=0.0) = GaussianProcess(variogram, mean)

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
