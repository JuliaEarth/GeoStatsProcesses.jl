# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    GaussianProcess(function, mean=0.0)

Gaussian process with given geostatistical `function` and global `mean`.

## Examples

```
GaussianProcess(GaussianVariogram())
GaussianProcess(SphericalCovariance(), 0.5)
```
"""
struct GaussianProcess{F,T} <: FieldProcess
  func::F
  mean::T
end

GaussianProcess(func) = GaussianProcess(func, 0.0)
GaussianProcess() = GaussianProcess(GaussianVariogram(), 0.0)

#---------
# METHODS
#---------

include("gaussian/lu.jl")
include("gaussian/fft.jl")
include("gaussian/seq.jl")

function defaultmethod(process::GaussianProcess, setup::RandSetup)
  f = process.func
  d = setup.domain
  p = parent(d)
  b = boundingbox(p)
  if p isa Grid && range(f) ≤ minimum(sides(b)) / 3
    FFTMethod()
  elseif nelements(d) < 100 * 100
    LUMethod()
  else
    SEQMethod()
  end
end
