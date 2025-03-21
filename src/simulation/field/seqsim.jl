# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    SEQSIM(; [options])

The sequential simulation method introduced by Gomez-Hernandez 1993.

The method traverses all locations of the geospatial domain according to
a path, approximates the conditional distribution within a neighborhood
using a geostatistical model, and assigns a value to the center of the
neighborhood by sampling from this distribution.

## Options

* `path`         - Process path (default to `LinearPath()`)
* `minneighbors` - Minimum number of neighbors (default to `1`)
* `maxneighbors` - Maximum number of neighbors (default to `26`)
* `neighborhood` - Search neighborhood (default to `:range`)
* `distance`     - Distance used to find nearest neighbors (default to `Euclidean()`)

For each location in the process `path`, a maximum number of
neighbors `maxneighbors` is used to fit the conditional Gaussian
distribution. The neighbors are searched according to a `neighborhood`.

The `neighborhood` can be a `MetricBall`, the symbol `:range` or `nothing`.
The symbol `:range` is converted to `MetricBall(range(f))` where `f` is the
geostatistical function of the Gaussian process. If `neighborhood` is `nothing`,
the nearest neighbors are used without additional constraints.

## References

* Gomez-Hernandez & Journel 1993. [Joint Sequential Simulation of
  MultiGaussian Fields](https://link.springer.com/chapter/10.1007/978-94-011-1739-5_8)

### Notes

This method is very sensitive to neighbor search options and
simulation path. Care must be taken to make sure that enough
neighbors are used in the underlying geostatistical model.
"""
@kwdef struct SEQSIM{P,N,D} <: FieldSimulationMethod
  path::P = LinearPath()
  minneighbors::Int = 1
  maxneighbors::Int = 26
  neighborhood::N = :range
  distance::D = Euclidean()
end

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
  sdom, sdat, cache = _transform(process, dom, dat)

  (; params, model, prior, sdom, sdat, cache, init)
end

function randsingle(rng::AbstractRNG, process, ::SEQSIM, domain, data, preproc)
  # retrieve preprocessing results
  (; params, model, prior, sdom, sdat, cache, init) = preproc

  # retrieve search parameters
  (; path, searcher, nmin, nmax) = params

  # initialize realization and mask
  real, mask = randinit(process, sdom, sdat, init)

  # realization in matrix form for efficient updates
  realization = ustrip.(stack(Tables.rowtable(real)))

  # save units of columns to restore later
  units = isnothing(sdat) ? _units(process) : _units(sdat)

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
    real[var] .= realization[i,:] * units[i]
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
  # consider point set of centroids for now
  sdom = PointSet([centroid(dom, ind) for ind in 1:nelements(dom)])

  # transform data and always produce valid cache
  sdat, cache = if isnothing(dat)
    dat, _cache(process)
  else
    _fwdtransform(process, dat)
  end

  sdom, sdat, cache
end

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

  # https://github.com/JuliaStats/Distributions.jl/issues/1413
  u = unit(eltype(μ))
  μᵤ = ustrip.(u, μ)
  Σᵤ = ustrip.(u^2, sill(f))

  model = Kriging(f, μᵤ)
  prior = nvariates(f) > 1 ? MvNormal(μᵤ, Σᵤ) : Normal(μᵤ, √Σᵤ)

  model, prior
end

function _probmodel(process::IndicatorProcess, func)
  f = func
  p = process.prob

  model = Kriging(f, p)
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

# ---------------
# PHYSICAL UNITS
# ---------------

function _units(process::GaussianProcess)
  μ = process.mean
  ntuple(i -> unit(μ[i]), length(μ))
end

function _units(process::IndicatorProcess)
  p = process.prob
  ntuple(i -> NoUnits, length(p))
end

function _units(data)
  types = Tables.schema(values(data)).types
  ntuple(i -> unit(types[i]), length(types))
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
