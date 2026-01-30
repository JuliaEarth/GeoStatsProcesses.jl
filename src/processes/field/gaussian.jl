# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    GaussianProcess(func, mean=0)

Gaussian process with given geostatistical `func`tion
(e.g. variogram) and `mean`.

## Examples

```julia
# univariate processes
GaussianProcess(GaussianVariogram())
GaussianProcess(SphericalCovariance(), 0.5)

# multivariate processes
GaussianProcess(LinearTransiogram(), [0.5, 0.5])
```
"""
struct GaussianProcess{F,M} <: FieldProcess
  func::F
  mean::M

  function GaussianProcess{F,M}(func, mean) where {F,M}
    nm = length(mean)
    nf = nvariates(func)
    @assert nm == nf "mean must have $nf components, received $nm"
    new(func, mean)
  end
end

GaussianProcess(func, mean) = GaussianProcess{typeof(func),typeof(mean)}(func, mean)
GaussianProcess(func) = GaussianProcess(func, _zeros(nvariates(func)))

_zeros(n) = n > 1 ? zeros(n) : 0.0

iscontinuous(process::GaussianProcess) = true

isanalytical(process::GaussianProcess) = true

function defaultschema(process::GaussianProcess)
  nvars = nvariates(process.func)
  names = nvars > 1 ? ntuple(i -> Symbol(:field, i), nvars) : (:field,)
  types = ntuple(i -> typeof(process.mean[i]), nvars)
  Tables.Schema(names, types)
end
