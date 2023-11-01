# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    GaussianProcess(variogram=GaussianVariogram(), mean=0.0)

Gaussian process with given theoretical `variogram` and global `mean`.
"""
struct GaussianProcess{V,T} <: FieldProcess
  variogram::V
  mean::T
end

GaussianProcess(variogram) = GaussianProcess(variogram, 0.0)
GaussianProcess() = GaussianProcess(GaussianVariogram(), 0.0)

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
