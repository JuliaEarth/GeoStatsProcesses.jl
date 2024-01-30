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

function defaultmethod(process::GaussianProcess, setup::RandSetup)
  γ = process.variogram
  d = setup.domain
  p = parent(d)
  b = boundingbox(p)
  if p isa Grid && range(γ) ≤ minimum(sides(b)) / 3
    FFTMethod()
  elseif nelements(d) < 100 * 100
    LUMethod()
  else
    SEQMethod()
  end
end
