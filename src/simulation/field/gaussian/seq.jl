# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function preprocess(::AbstractRNG, process::GaussianProcess, method::SEQSIM, init, domain, data)
  # function and mean
  f = process.func
  μ = process.mean

  # method options
  minneighbors = method.minneighbors
  maxneighbors = method.maxneighbors
  neighborhood = method.neighborhood
  distance = method.distance

  # scale domains for numerical stability
  finv, dom, dat = _scaledomains(domain, data)

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
  μ̄, σ̄ = μ, _std(f̄)
  probmodel = Kriging(f̄, μ̄)
  gaussians = map((μ̄ᵢ, σ̄ᵢ) -> Normal(μ̄ᵢ, σ̄ᵢ), μ̄, σ̄)
  marginal = product_distribution(gaussians)

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

  (; dom, dat, init, probmodel, marginal, minneighbors, maxneighbors, searcher)
end

function randsingle(rng::AbstractRNG, process::GaussianProcess, method::SEQSIM, domain, data, preproc)
  # retrieve preprocessing results
  (; dom, dat, init, probmodel, marginal, minneighbors, maxneighbors, searcher) = preproc

  # process parameters
  func = process.func

  # method options
  path = method.path

  # initialize realization and mask
  real, mask = randinit(process, dom, dat, init)

  # retrieve variable names
  vars = keys(real)

  # sanity checks
  if length(vars) != nvariates(func)
    throw(ArgumentError("incompatible number of variables for geostatistical function"))
  end

  # consider point set with centroids for now
  pointset = PointSet([centroid(dom, ind) for ind in 1:nelements(dom)])

  # pre-allocate memory for neighbors
  neighbors = Vector{Int}(undef, maxneighbors)

  # matrix buffer for realization
  realization = stack(Tables.rowtable(real))

  # locations with all variables already simulated
  simulated = map(all, Tables.rowtable(mask))

  # simulation loop
  @inbounds for ind in traverse(dom, path)
    if !simulated[ind]
      # center of target location
      center = pointset[ind]

      # search neighbors with simulated data
      nneigh = search!(neighbors, center, searcher, mask=simulated)

      if nneigh < minneighbors
        # draw from marginal
        realization[:,ind] .= rand(rng, marginal)
      else
        # neighborhood with data
        neigh = let
          ninds = view(neighbors, 1:nneigh)
          ndom = view(pointset, ninds)
          nmat = view(realization, :, ninds)
          ntab = (; zip(vars, eachrow(nmat))...)
          georef(ntab, ndom)
        end

        # fit probabilistic model
        fitted = GeoStatsModels.fit(probmodel, neigh)

        # draw from conditional
        conditional = if GeoStatsModels.status(fitted)
          dists = GeoStatsModels.predictprob(fitted, vars, center)
          product_distribution(dists)
        else
          marginal
        end
        realization[:,ind] .= rand(rng, conditional)
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

# -----------------
# HELPER FUNCTIONS
# -----------------

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

_std(f) = ustrip.(sqrt.(_diag(sill(f))))
_diag(M::Matrix) = diag(M)
_diag(M) = M
