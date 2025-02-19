# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    FieldSimulationMethod

A method for simulating field processes.
"""
abstract type FieldSimulationMethod end

"""
    LUSIM(; [options])

The LU simulation method introduced by Alabert 1987.

The full covariance matrix is built to include all locations
of the data, and samples from the multivariate Gaussian are
drawn via a lower-upper (LU) matrix factorization.

## Options

`factorization` - Factorization function from LinearAlgebra (default to `cholesky`)

## References

* Alabert 1987. [The practice of fast conditional simulations
  through the LU decomposition of the covariance matrix]
  (https://link.springer.com/article/10.1007/BF00897191)

* Oliver 2003. [Gaussian cosimulation: modeling of the cross-covariance]
  (https://link.springer.com/article/10.1023%2FB%3AMATG.0000002984.56637.ef)

### Notes

The method is only adequate for domains with relatively small
number of elements (e.g. 100x100 grids) where it is feasible to
factorize the full covariance.

For larger domains (e.g. 3D grids), other methods are preferred
such as [`SEQSIM`](@ref) and [`FFTSIM`](@ref).
"""
@kwdef struct LUSIM{F} <: FieldSimulationMethod
  factorization::F = cholesky
end

"""
    SEQSIM(; [options])

The sequential simulation method introduced by Gomez-Hernandez 1993.

The method traverses all locations of the geospatial domain according to
a path, approximates the conditional distribution within a neighborhood
using a geostatistical model, and assigns a value to the center of the
neighborhood by sampling from this distribution.

## Options

* `path`         - Process path (default to `LinearPath()`)
* `minneighbors` - Minimum number of neighbors (default to `1`)
* `maxneighbors` - Maximum number of neighbors (default to `26`)
* `neighborhood` - Search neighborhood (default to `:range`)
* `distance`     - Distance used to find nearest neighbors (default to `Euclidean()`)

For each location in the process `path`, a maximum number of
neighbors `maxneighbors` is used to fit the conditional Gaussian
distribution. The neighbors are searched according to a `neighborhood`.

The `neighborhood` can be a `MetricBall`, the symbol `:range` or `nothing`.
The symbol `:range` is converted to `MetricBall(range(f))` where `f` is the
geostatistical function of the Gaussian process. If `neighborhood` is `nothing`,
the nearest neighbors are used without additional constraints.

## References

* Gomez-Hernandez & Journel 1993. [Joint Sequential Simulation of
  MultiGaussian Fields](https://link.springer.com/chapter/10.1007/978-94-011-1739-5_8)

### Notes

This method is very sensitive to the neighbor search options.
Care must be taken to make sure that enough neighbors are used
in the underlying geostatistical model.
"""
@kwdef struct SEQSIM{P,N,D} <: FieldSimulationMethod
  path::P = LinearPath()
  minneighbors::Int = 1
  maxneighbors::Int = 26
  neighborhood::N = :range
  distance::D = Euclidean()
end

"""
    FFTSIM(; [options])

The FFT simulation method introduced by Gutjahr 1997.

The covariance function is perturbed in the frequency domain
after a fast Fourier transform. White noise is added to the
phase of the spectrum, and a realization is produced by an
inverse Fourier transform.

## Options

* `minneighbors` - Minimum number of neighbors (default to `1`)
* `maxneighbors` - Maximum number of neighbors (default to `26`)
* `neighborhood` - Search neighborhood (default to `nothing`)
* `distance`     - Distance used to find nearest neighbors (default to `Euclidean()`)

## References

* Gutjahr 1997. [General joint conditional simulations using a fast
  Fourier transform method](https://link.springer.com/article/10.1007/BF02769641)

* Gómez-Hernández, J. & Srivastava, R. 2021. [One Step at a Time: The Origins
  of Sequential Simulation and Beyond](https://link.springer.com/article/10.1007/s11004-021-09926-0)

### Notes

The method is limited to simulations on regular grids, and care must be
taken to make sure that the correlation length is small enough compared to
the grid size. As a general rule of thumb, avoid correlation lengths greater
than 1/3 of the grid.

Visual artifacts can appear near the boundaries of the grid if the correlation
length is large compared to the grid itself.
"""
@kwdef struct FFTSIM{N,D} <: FieldSimulationMethod
  minneighbors::Int = 1
  maxneighbors::Int = 26
  neighborhood::N = nothing
  distance::D = Euclidean()
end
