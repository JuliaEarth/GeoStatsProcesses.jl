# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------


"""
    SEQ([paramaters])

A sequential simulation solver.

For each location in the simulation `path`, a maximum number
of neighbors `maxneighbors` is used to fit a distribution.
The neighbors are searched according to a `neighborhood`,
and in case there are none, use a `marginal` distribution.

## Parameters

* `estimator`    - CDF estimator
* `marginal`     - Marginal distribution
* `path`         - Simulation path (default to `LinearPath()`)
* `minneighbors` - Minimum number of neighbors (default to `1`)
* `maxneighbors` - Maximum number of neighbors (default to `10`)
* `neighborhood` - Search neighborhood (default to `nothing`)
* `distance`     - Distance used to find nearest neighbors (default to `Euclidean()`)
* `init`         - Data initialization method (default to `NearestInit()`)
"""
@kwdef struct SEQ{M,D,P,N,DT,I} <: FieldProcess
  probmodel::M
  marginal::D
  path::P = LinearPath()
  minneighbors::Int = 1
  maxneighbors::Int = 10
  neighborhood::N = nothing
  distance::DT = Euclidean()
  init::I = NearestInit()
end

"""
    SGP([paramaters])

The sequential Gaussian simulation solver introduced by Gomez-Hernandez 1993.
It traverses all locations of the geospatial domain according to a path,
approximates the conditional Gaussian distribution within a neighborhood
using Kriging, and assigns a value to the center of the neighborhood by
sampling from this distribution.

## Parameters

* `variogram`    - Variogram model (default to `GaussianVariogram()`)
* `mean`         - mean for simple Kriging (default to `0.0`)
* `path`         - Simulation path (default to `LinearPath()`)
* `minneighbors` - Minimum number of neighbors (default to `1`)
* `maxneighbors` - Maximum number of neighbors (default to `10`)
* `neighborhood` - Search neighborhood (default to `nothing`)
* `distance`     - Distance used to find nearest neighbors (default to `Euclidean()`)
* `init`         - Data initialization method (default to `NearestInit()`)

For each location in the simulation `path`, a maximum number of
neighbors `maxneighbors` is used to fit the conditional Gaussian
distribution. The neighbors are searched according to a `neighborhood`.

### References

* Gomez-Hernandez & Journel 1993. [Joint Sequential Simulation of
  MultiGaussian Fields](https://link.springer.com/chapter/10.1007/978-94-011-1739-5_8)

### Notes

* This solver is very sensitive to the simulation path and number of
  samples. Care must be taken to make sure that neighborhoods have
  enough samples for the Kriging estimator.
"""
function SGP(; variogram=GaussianVariogram(), mean=0.0, kwargs...)
  probmodel = GeoStatsModels.SimpleKriging(variogram, mean)
  marginal = Normal(mean, √sill(variogram))
  SEQ(; probmodel, marginal, kwargs...)
end

function randprep(::AbstractRNG, process::SEQ, setup::RandSetup)
  # retrieve domain info
  domain = setup.domain
  # retrieve process paramaters
  (; minneighbors, maxneighbors, neighborhood, distance) = process

  nelems = nelements(domain)
  if maxneighbors > nelems || maxneighbors < 1
    @warn "Invalid maximum number of neighbors. Adjusting to $nelems..."
    maxneighbors = nelems
  end

  if minneighbors > maxneighbors || minneighbors < 1
    @warn "Invalid minimum number of neighbors. Adjusting to 1..."
    minneighbors = 1
  end

  # determine search method
  searcher = if isnothing(neighborhood)
    # nearest neighbor search with a metric
    KNearestSearch(domain, maxneighbors; metric=distance)
  else
    # neighbor search with ball neighborhood
    KBallSearch(domain, maxneighbors, neighborhood)
  end

  (; minneighbors, maxneighbors, searcher)
end

function randsingle(rng::AbstractRNG, process::SEQ, setup::RandSetup, prep)
  # retrieve parameters
  (; probmodel, marginal, path, init) = process
  (; domain, geotable, varnames, vartypes) = setup
  (; minneighbors, maxneighbors, searcher) = prep

  # initialize buffers for realization and simulation mask
  vars = Dict(zip(varnames, vartypes))
  buff, mask = initbuff(domain, vars, init, data=geotable)

  # consider point set with centroids for now
  pset = PointSet([centroid(domain, ind) for ind in 1:nelements(domain)])

  varreals = map(varnames) do var
    # pre-allocate memory for neighbors
    neighbors = Vector{Int}(undef, maxneighbors)

    # retrieve realization and mask for variable
    realization = buff[var]
    simulated = mask[var]

    # simulation loop
    for ind in traverse(domain, path)
      if !simulated[ind]
        center = pset[ind]
        # search neighbors with simulated data
        nneigh = search!(neighbors, center, searcher, mask=simulated)

        if nneigh < minneighbors
          # draw from marginal
          realization[ind] = rand(rng, marginal)
        else
          # neighborhood with data
          neigh = let
            ninds = view(neighbors, 1:nneigh)
            dom = view(pset, ninds)
            val = view(realization, ninds)
            tab = (; var => val)
            georef(tab, dom)
          end

          # fit distribution probmodel
          fitted = GeoStatsModels.fit(probmodel, neigh)

          # draw from conditional or marginal
          distribution = if GeoStatsModels.status(fitted)
            GeoStatsModels.predictprob(fitted, var, center)
          else
            marginal
          end
          realization[ind] = rand(rng, distribution)
        end

        # mark location as simulated and continue
        simulated[ind] = true
      end
    end

    var => realization
  end

  Dict(varreals)
end
