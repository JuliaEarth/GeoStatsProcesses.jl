# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function preprocess(::AbstractRNG, process, method::SEQSIM, init, domain, data)
  # function and mean
  func = process.func
  mean = process.mean

  # scale objects for numerical stability
  dom, dat, fun, neigh = _scale(domain, data, func, method.neighborhood)

  # determine search method and min/max neighbors
  path, searcher, nmin, nmax = _search(dom, neigh, method)

  # store adjusted parameters
  params = (; path, searcher, nmin, nmax)

  # determine probability model
  model, prior = _probmodel(process, fun, mean)

  # additional transformations
  sdom = PointSet([centroid(dom, ind) for ind in 1:nelements(dom)])
  sdat = dat

  (; params, model, prior, init, sdom, sdat)
end

function randsingle(rng::AbstractRNG, process, ::SEQSIM, domain, data, preproc)
  # retrieve preprocessing results
  (; params, model, prior, init, sdom, sdat) = preproc

  # retrieve search parameters
  (; path, searcher, nmin, nmax) = params

  # initialize realization and mask
  real, mask = randinit(process, domain, data, init)

  # retrieve variable names
  vars = keys(real)

  # variables passed to probability model
  mvars = length(vars) == 1 ? first(vars) : vars

  # pre-allocate memory for neighbors
  neighbors = Vector{Int}(undef, nmax)

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
      n = search!(neighbors, center, searcher, mask=simulated)

      if n < nmin
        # draw from prior
        _draw!(rng, prior, buffer)
      else
        # neighborhood with data
        neigh = let
          inds = view(neighbors, 1:n)
          ndom = view(sdom, inds)
          nmat = view(realization, :, inds)
          ntab = (; zip(vars, eachrow(nmat))...)
          georef(ntab, ndom)
        end

        # fit probability model
        fitted = GeoStatsModels.fit(model, neigh)

        # draw from conditional
        conditional = if GeoStatsModels.status(fitted)
          _conditional(process, fitted, mvars, center)
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

# --------
# SCALING
# --------

function _scale(dom, dat, fun, neigh)
  α₁ = _scalefactor(dom)
  α₂ = isnothing(dat) ? 1 : _scalefactor(domain(dat))
  α₃ = _scalefactor(fun)
  α = inv(max(α₁, α₂, α₃))

  sdom = dom |> Scale(α)
  sdat = isnothing(dat) ? nothing : (dat |> Scale(α))
  sfun = GeoStatsFunctions.scale(fun, α)
  sneigh = if neigh isa MetricBall
    α * neigh
  elseif neigh == :range
    MetricBall(range(sfun))
  else
    nothing
  end

  sdom, sdat, sfun, sneigh
end

function _scalefactor(domain::Domain)
  pmin, pmax = extrema(boundingbox(domain))
  cmin = abs.(to(pmin))
  cmax = abs.(to(pmax))
  ustrip(max(cmin..., cmax...))
end

_scalefactor(fun::GeoStatsFunction) = ustrip(range(fun))

# ----------
# SEARCHING
# ----------

function _search(dom, neigh, method)
  path = method.path
  nmin = method.minneighbors
  nmax = method.maxneighbors
  dist = method.distance

  # adjust min/max neighbors
  nelm = nelements(dom)
  if nmax > nelm || nmax < 1
    nmax = nelm
  end
  if nmin > nmax || nmin < 1
    nmin = 1
  end

  # determine search method
  searcher = if isnothing(neigh)
    # nearest neighbor search with a metric
    KNearestSearch(dom, nmax; metric=dist)
  else
    # neighbor search with ball neighborhood
    KBallSearch(dom, nmax, neigh)
  end

  path, searcher, nmin, nmax
end

# ------------------
# PROBABILITY MODEL
# ------------------

function _probmodel(::GaussianProcess, fun, mean)
  f̄ = fun
  μ̄ = mean
  Σ̄ = ustrip.(sill(f̄))

  model = Kriging(f̄, μ̄)
  prior = nvariates(f̄) > 1 ? MvNormal(μ̄, Σ̄) : Normal(μ̄, √Σ̄)

  model, prior
end

function _conditional(::GaussianProcess, fitted, vars, geom)
  GeoStatsModels.predictprob(fitted, vars, geom)
end

# --------
# DRAWING
# --------

function _draw!(rng, dist::Distribution, buffer)
  val = rand(rng, dist)
  @inbounds for i in eachindex(buffer)
    buffer[i] = val[i]
  end
end

function _draw!(rng, dist::Categorical, buffer)
  j = rand(rng, dist)
  @inbounds for i in eachindex(buffer)
    buffer[i] = (i == j)
  end
end
