# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function preprocess(::AbstractRNG, process::GaussianProcess, method::SEQSIM, init, domain, data)
  # function and mean
  func = process.func
  mean = process.mean

  # method options
  path = method.path
  minneighbors = method.minneighbors
  maxneighbors = method.maxneighbors
  neighborhood = method.neighborhood
  distance = method.distance

  # scale domains and function for numerical stability
  dom, dat, f̄, α = _scale(domain, data, func)

  # scale neighborhood accordingly
  neigh = if neighborhood isa MetricBall
    α * neighborhood
  elseif neighborhood == :range
    MetricBall(range(f̄))
  else
    nothing
  end

  # adjust min/max neighbors
  nelm = nelements(dom)
  if maxneighbors > nelm || maxneighbors < 1
    maxneighbors = nelm
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

  # store adjusted parameters
  params = (; path, minneighbors, maxneighbors, searcher)

  # determine probability model
  μ̄, Σ̄ = mean, ustrip.(sill(f̄))
  model = Kriging(f̄, μ̄)
  prior = nvariates(f̄) > 1 ? MvNormal(μ̄, Σ̄) : Normal(μ̄, √Σ̄)

  # additional transformations
  sdom = PointSet([centroid(dom, ind) for ind in 1:nelements(dom)])
  sdat = dat

  (; params, model, prior, init, sdom, sdat)
end

function randsingle(rng::AbstractRNG, process::GaussianProcess, ::SEQSIM, domain, data, preproc)
  # retrieve preprocessing results
  (; params, model, prior, init, sdom, sdat) = preproc

  # retrieve search parameters
  (; path, minneighbors, maxneighbors, searcher) = params

  # initialize realization and mask
  real, mask = randinit(process, domain, data, init)

  # retrieve variable names
  vars = keys(real)

  # variables passed to probability model
  mvars = length(vars) == 1 ? first(vars) : vars

  # pre-allocate memory for neighbors
  neighbors = Vector{Int}(undef, maxneighbors)

  # matrix buffer for realization
  realization = stack(Tables.rowtable(real))

  # locations with all variables already simulated
  simulated = map(all, Tables.rowtable(mask))

  # simulation loop
  @inbounds for ind in traverse(domain, path)
    if !simulated[ind]
      # center of target location
      center = centroid(sdom, ind)

      # buffer at target location
      buffer = view(realization, :, ind)

      # search neighbors with simulated data
      nneigh = search!(neighbors, center, searcher, mask=simulated)

      if nneigh < minneighbors
        # draw from prior
        _draw!(rng, prior, buffer)
      else
        # neighborhood with data
        neigh = let
          ninds = view(neighbors, 1:nneigh)
          ndom = view(sdom, ninds)
          nmat = view(realization, :, ninds)
          ntab = (; zip(vars, eachrow(nmat))...)
          georef(ntab, ndom)
        end

        # fit probabilistic model
        fitted = GeoStatsModels.fit(model, neigh)

        # draw from conditional
        conditional = if GeoStatsModels.status(fitted)
          GeoStatsModels.predictprob(fitted, mvars, center)
        else
          prior
        end
        _draw!(rng, conditional, buffer)
      end

      # mark location as simulated and continue
      simulated[ind] = true
    end
  end

  # convert back to table format
  @inbounds for (i, var) in enumerate(vars)
    real[var] .= realization[i,:]
  end

  real
end
