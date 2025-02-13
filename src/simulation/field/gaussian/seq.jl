# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function preprocess(::AbstractRNG, process::GaussianProcess, method::SEQSIM, domain, data)
  # function and mean
  f = process.func
  μ = process.mean

  # scale domains for numerical stability
  finv, dom, data = _scaledomains(domain, data)

  # scale function accordingly
  f̄ = GeoStatsFunctions.scale(f, finv)

  # scale neighborhood accordingly
  neigh = if neighborhood isa MetricBall
    finv * neighborhood
  elseif neighborhood == :range
    MetricBall(range(f̄))
  else
    nothing
  end

  # determine probability model
  probmodel = Kriging(f̄, μ)
  marginal = Normal(μ, √sill(f̄))

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

function randsingle(rng::AbstractRNG, ::GaussianProcess, method::SEQSIM, domain, data, preproc)
  # retrieve parameters
  (; path, init) = method
  (; dom, data, probmodel, marginal, minneighbors, maxneighbors, searcher) = preproc

  # initialize buffers for realization and simulation mask
  vars = Dict(zip(varnames, vartypes))
  buff, mask = initbuff(dom, vars, init, data=data)

  # consider point set with centroids for now
  pointset = PointSet([centroid(dom, ind) for ind in 1:nelements(dom)])

  pairs = map(varnames) do var
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

  (; pairs...)
end

function _scaledomains(domain, data)
  sdom = domain
  sdat = data
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
  cmin = abs.(to(pmin))
  cmax = abs.(to(pmax))
  ustrip(max(cmin..., cmax...))
end
