# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    SequentialProcess(probmodel, marginal; [paramaters])

A sequential simulation process.

For each location in the simulation `path`, a maximum number
of neighbors `maxneighbors` is used to fit a distribution.
The neighbors are searched according to a `neighborhood`,
and in case there are none, use a `marginal` distribution.

## Parameters

### Required

* `probmodel`    - Conditional distribution model from GeoStatsModels.jl
* `marginal`     - Marginal distribution from Distributions.jl

### Optional

* `path`         - Simulation path (default to `LinearPath()`)
* `minneighbors` - Minimum number of neighbors (default to `1`)
* `maxneighbors` - Maximum number of neighbors (default to `10`)
* `neighborhood` - Search neighborhood (default to `nothing`)
* `distance`     - Distance used to find nearest neighbors (default to `Euclidean()`)
* `init`         - Data initialization method (default to `NearestInit()`)
"""
struct SequentialProcess{M,MD,P,N,D,I} <: FieldProcess
  probmodel::M
  marginal::MD
  path::P
  minneighbors::Int
  maxneighbors::Int
  neighborhood::N
  distance::D
  init::I
end

SequentialProcess(
  probmodel,
  marginal;
  path=LinearPath(),
  minneighbors=1,
  maxneighbors=10,
  neighborhood=nothing,
  distance=Euclidean(),
  init=NearestInit()
) = SequentialProcess(probmodel, marginal, path, minneighbors, maxneighbors, neighborhood, distance, init)

function randprep(::AbstractRNG, process::SequentialProcess, ::DefaultRandMethod, setup::RandSetup)
  # retrieve domain info
  domain = setup.domain
  # retrieve process paramaters
  (; minneighbors, maxneighbors, neighborhood, distance) = process

  nobs = nelements(domain)
  if maxneighbors > nobs || maxneighbors < 1
    maxneighbors = nobs
  end
  if minneighbors > maxneighbors || minneighbors < 1
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

function randsingle(rng::AbstractRNG, process::SequentialProcess, ::DefaultRandMethod, setup::RandSetup, prep)
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
