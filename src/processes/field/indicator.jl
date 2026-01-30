# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    IndicatorProcess(func, prob)

Indicator process with multivariate geostatistical `func`tion
(e.g., multivariate covariance) and a priori `prob`abilities.

    IndicatorProcess(transiogram, prob=proportions(transiogram))

Alternatively, construct indicator process from `transiogram`.

## Examples

```julia
# covariance and explicit probabilities
IndicatorProcess([1.0 0.5; 0.5 1.0] * SphericalCovariance(), [0.8, 0.2])

# transiogram and implicit probabilities
IndicatorProcess(LinearTransiogram())

# transiogram and explicit probabilities
IndicatorProcess(ExponentialTransiogram(), [0.8, 0.2])
```
"""
struct IndicatorProcess{F,P} <: FieldProcess
  func::F
  prob::P

  function IndicatorProcess{F,P}(func, prob) where {F,P}
    np = length(prob)
    nf = nvariates(func)
    @assert nf > 1 "indicator process requires multivariate function"
    @assert nf == np "probabilities must have $nf components, received $np"
    @assert all(p -> 0 ≤ p ≤ 1, prob) "probabilities must be in [0, 1]"
    @assert sum(prob) ≈ 1 "probabilities must sum up to one"
    new(func, prob)
  end
end

IndicatorProcess(func, prob) = IndicatorProcess{typeof(func),typeof(prob)}(func, prob)

IndicatorProcess(func::Transiogram) = IndicatorProcess(func, collect(proportions(func)))

function defaultschema(process::IndicatorProcess)
  nvars = nvariates(process.func)
  names = ntuple(i -> Symbol(:field, i), nvars)
  types = ntuple(i -> Bool, nvars)
  Tables.Schema(names, types)
end
