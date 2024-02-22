# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    SEQMethod(; [paramaters])

The sequential process method introduced by Gomez-Hernandez 1993.
It traverses all locations of the geospatial domain according to a path,
approximates the conditional Gaussian distribution within a neighborhood
using simple Kriging, and assigns a value to the center of the neighborhood
by sampling from this distribution.

## Parameters

* `path`         - Process path (default to `LinearPath()`)
* `minneighbors` - Minimum number of neighbors (default to `1`)
* `maxneighbors` - Maximum number of neighbors (default to `36`)
* `neighborhood` - Search neighborhood (default to `:range`)
* `distance`     - Distance used to find nearest neighbors (default to `Euclidean()`)
* `init`         - Data initialization method (default to `NearestInit()`)

For each location in the process `path`, a maximum number of
neighbors `maxneighbors` is used to fit the conditional Gaussian
distribution. The neighbors are searched according to a `neighborhood`.

The `neighborhood` can be a `MetricBall`, the symbol `:range` or `nothing`.
The symbol `:range` is converted to `MetricBall(range(γ))` where `γ` is the
variogram of the Gaussian process. If `neighborhood` is `nothing`, the nearest
neighbors are used without additional constraints.

## References

* Gomez-Hernandez & Journel 1993. [Joint Sequential Simulation of
  MultiGaussian Fields](https://link.springer.com/chapter/10.1007/978-94-011-1739-5_8)

### Notes

* This method is very sensitive to the various parameters.
  Care must be taken to make sure that enough neighbors
  are used in the underlying Kriging model.
"""
@kwdef struct SEQMethod{P,N,D,I} <: RandMethod
  path::P = LinearPath()
  minneighbors::Int = 1
  maxneighbors::Int = 36 # 6x6 grid cells
  neighborhood::N = :range
  distance::D = Euclidean()
  init::I = NearestInit()
end

function randprep(::AbstractRNG, process::GaussianProcess, method::SEQMethod, setup::RandSetup)
  # retrieve paramaters
  (; variogram, mean) = process
  (; minneighbors, maxneighbors, neighborhood, distance) = method

  # scale domains for numerical stability
  finv, dom, data = _scaledomains(setup)

  # scale variogram model accordingly
  gamma = GeoStatsFunctions.scale(variogram, finv)

  # scale neighborhood accordingly
  neigh = if neighborhood isa MetricBall
    finv * neighborhood
  elseif neighborhood == :range
    MetricBall(range(gamma))
  else
    nothing
  end

  # determine probability model
  probmodel = GeoStatsModels.SimpleKriging(gamma, mean)
  marginal = Normal(mean, √sill(gamma))

  # adjust min/max neighbors
  nobs = nelements(dom)
  if maxneighbors > nobs || maxneighbors < 1
    maxneighbors = nobs
  end
  if minneighbors > maxneighbors || minneighbors < 1
    minneighbors = 1
  end

  # determine search method
  searcher = if isnothing(neigh)
    # nearest neighbor search with a metric
    KNearestSearch(dom, maxneighbors; metric=distance)
  else
    # neighbor search with ball neighborhood
    KBallSearch(dom, maxneighbors, neigh)
  end

  (; dom, data, probmodel, marginal, minneighbors, maxneighbors, searcher)
end

function randsingle(rng::AbstractRNG, ::GaussianProcess, method::SEQMethod, setup::RandSetup, prep)
  # retrieve parameters
  (; path, init) = method
  (; varnames, vartypes) = setup
  (; dom, data, probmodel, marginal, minneighbors, maxneighbors, searcher) = prep

  # initialize buffers for realization and simulation mask
  vars = Dict(zip(varnames, vartypes))
  buff, mask = initbuff(dom, vars, init, data=data)

  # consider point set with centroids for now
  pointset = PointSet([centroid(dom, ind) for ind in 1:nelements(dom)])

  varreals = map(varnames) do var
    # pre-allocate memory for neighbors
    neighbors = Vector{Int}(undef, maxneighbors)

    # retrieve realization and mask for variable
    realization = buff[var]
    simulated = mask[var]

    # simulation loop
    for ind in traverse(dom, path)
      if !simulated[ind]
        center = pointset[ind]
        # search neighbors with simulated data
        nneigh = search!(neighbors, center, searcher, mask=simulated)

        if nneigh < minneighbors
          # draw from marginal
          realization[ind] = rand(rng, marginal)
        else
          # neighborhood with data
          neigh = let
            ninds = view(neighbors, 1:nneigh)
            dom = view(pointset, ninds)
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

function _scaledomains(setup)
  sdom = setup.domain
  sdat = setup.geotable
  fdom = _scalefactor(sdom)
  fdat = isnothing(sdat) ? 1 : _scalefactor(domain(sdat))
  fmax = max(fdom, fdat)
  finv = 1 / fmax
  func = Scale(finv)

  dom = func(sdom)
  dat = isnothing(sdat) ? nothing : func(sdat)

  finv, dom, dat
end

function _scalefactor(domain)
  pmin, pmax = extrema(boundingbox(domain))
  cmin = abs.(coordinates(pmin))
  cmax = abs.(coordinates(pmax))
  max(cmin..., cmax...)
end
