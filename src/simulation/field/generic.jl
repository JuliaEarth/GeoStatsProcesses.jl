# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function preprocess(::AbstractRNG, process, method::SEQSIM, init, domain, data)
  # geostatistical function
  func = process.func

  # scale objects for numerical stability
  dom, dat, fun, neigh = _scale(domain, data, func, method.neighborhood)

  # determine search method and min/max neighbors
  path, searcher, nmin, nmax = _search(dom, neigh, method)

  # store adjusted parameters
  params = (; path, searcher, nmin, nmax)

  # determine probability model
  model, prior = _probmodel(process, fun)

  # transform process and data
  sproc, sdom, sdat, cache = _transform(process, dom, dat)

  (; params, model, prior, sproc, sdom, sdat, cache, init)
end

function randsingle(rng::AbstractRNG, process, ::SEQSIM, domain, data, preproc)
  # retrieve preprocessing results
  (; params, model, prior, sproc, sdom, sdat, cache, init) = preproc

  # retrieve search parameters
  (; path, searcher, nmin, nmax) = params

  # initialize realization and mask
  real, mask = randinit(sproc, sdom, sdat, init)

  # realization in matrix form for efficient updates
  realization = stack(Tables.rowtable(real))

  # locations with all variables already simulated
  simulated = map(all, Tables.rowtable(mask))

  # pre-allocate memory for neighbors
  neighbors = Vector{Int}(undef, nmax)

  # retrieve variable names
  vars = keys(real)

  # variables passed to probability model
  mvars = length(vars) == 1 ? first(vars) : vars

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

  # undo data transformations
  rdat = _bwdtransform(process, georef(real, sdom), cache)

  # return realization values
  values(rdat)
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

# ----------------
# TRANSFORMATIONS
# ----------------

function _transform(process, dom, dat)
  # transformed process for sequential simulation
  sproc = _transformed(process)

  # consider point set of centroids for now
  sdom = PointSet([centroid(dom, ind) for ind in 1:nelements(dom)])

  # transform data and always produce valid cache
  sdat, cache = if isnothing(dat)
    dat, _cache(process)
  else
    _fwdtransform(process, dat)
  end

  sproc, sdom, sdat, cache
end

_transformed(process::GaussianProcess) = process

_transformed(process::IndicatorProcess) = GaussianProcess(process.func, process.prob)

_fwdtransform(::GaussianProcess, dat) = apply(Identity(), dat)

_fwdtransform(::IndicatorProcess, dat) = apply(OneHot(1), dat)

_bwdtransform(::GaussianProcess, rdat, cache) = revert(Identity(), rdat, cache)

_bwdtransform(::IndicatorProcess, rdat, cache) = revert(OneHot(1), rdat, cache)

_cache(::GaussianProcess) = nothing

function _cache(process::IndicatorProcess)
  f = process.func
  n = nvariates(f)
  t = (field = 1:n,)
  apply(OneHot(1), t) |> last
end

# ------------------
# PROBABILITY MODEL
# ------------------

function _probmodel(process::GaussianProcess, func)
  f = func
  μ = process.mean
  Σ = ustrip.(sill(f))

  model = Kriging(f, μ)
  prior = nvariates(f) > 1 ? MvNormal(μ, Σ) : Normal(μ, √Σ)

  model, prior
end

function _probmodel(process::IndicatorProcess, func)
  f = func
  p = process.prob
  μ = zeros(length(p))

  model = Kriging(f, μ)
  prior = Categorical(p)

  model, prior
end

function _conditional(::GaussianProcess, fitted, vars, geom)
  GeoStatsModels.predictprob(fitted, vars, geom)
end

function _conditional(::IndicatorProcess, fitted, vars, geom)
  p = GeoStatsModels.predict(fitted, vars, geom)
  p̂ = normalize(clamp.(p, 0, 1), 1)
  Categorical(p̂)
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
