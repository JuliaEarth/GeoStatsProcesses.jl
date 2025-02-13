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
