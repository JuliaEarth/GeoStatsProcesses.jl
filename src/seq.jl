# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

@kwdef struct SEQ{M,D,P,N,DT,I} <: GeoStatsProcess
  probmodel::M
  marginal::D
  path::P = LinearPath()
  minneighbors::Int = 1
  maxneighbors::Int = 10
  neighborhood::N = nothing
  distance::DT = Euclidean()
  init::I = NearestInit()
end

function SGP(; variogram=GaussianVariogram(), mean=0.0, kwargs...)
  probmodel = GeoStatsprobmodels.SimpleKriging(variogram, mean)
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
