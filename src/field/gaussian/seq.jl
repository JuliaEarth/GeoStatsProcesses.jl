# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    SEQMethod(; [paramaters])

The sequential process method introduced by Gomez-Hernandez 1993.
It traverses all locations of the geospatial domain according to a path,
approximates the conditional Gaussian distribution within a neighborhood
using Kriging, and assigns a value to the center of the neighborhood by
sampling from this distribution.

## Parameters

* `path`         - Process path (default to `LinearPath()`)
* `minneighbors` - Minimum number of neighbors (default to `1`)
* `maxneighbors` - Maximum number of neighbors (default to `10`)
* `neighborhood` - Search neighborhood (default to `nothing`)
* `distance`     - Distance used to find nearest neighbors (default to `Euclidean()`)
* `init`         - Data initialization method (default to `NearestInit()`)

For each location in the process `path`, a maximum number of
neighbors `maxneighbors` is used to fit the conditional Gaussian
distribution. The neighbors are searched according to a `neighborhood`.

### References

* Gomez-Hernandez & Journel 1993. [Joint Sequential Simulation of
  MultiGaussian Fields](https://link.springer.com/chapter/10.1007/978-94-011-1739-5_8)

### Notes

* This method is very sensitive to the process path and number of
  samples. Care must be taken to make sure that neighborhoods have
  enough samples for the geostatistical model (e.g. Kriging).
"""
@kwdef struct SEQMethod{P,N,D,I} <: RandMethod
  path::P = LinearPath()
  minneighbors::Int = 1
  maxneighbors::Int = 10
  neighborhood::N = nothing
  distance::D = Euclidean()
  init::I = NearestInit()
end

function randprep(::AbstractRNG, process::GaussianProcess, method::SEQMethod, setup::RandSetup)
  (; variogram, mean) = process
  (; path, minneighbors, maxneighbors, neighborhood, distance, init) = method
  probmodel = GeoStatsModels.SimpleKriging(variogram, mean)
  marginal = Normal(mean, √sill(variogram))
  SequentialProcess(probmodel, marginal; path, minneighbors, maxneighbors, neighborhood, distance, init)
end

function randsingle(rng::AbstractRNG, ::GaussianProcess, ::SEQMethod, setup::RandSetup, prep)
  seq = prep
  seqmethod = DefaultRandMethod()
  seqprep = randprep(rng, seq, seqmethod, setup)
  randsingle(rng, seq, seqmethod, setup, seqprep)
end
