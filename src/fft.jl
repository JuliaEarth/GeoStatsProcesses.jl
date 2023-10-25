# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

@kwdef struct FFTGP{V,N,D} <: GeoStatsProcess
  variogram::V = GaussianVariogram()
  mean::Float64 = 0.0
  minneighbors::Int = 1
  maxneighbors::Int = 10
  neighborhood::N = nothing
  distance::D = Euclidean()
end

function randprep(::AbstractRNG, process::FFTGP, setup::RandSetup)
  # retrive variogram model and mean
  γ = process.variogram
  μ = process.mean

  # check stationarity
  if !isstationary(γ)
    throw(ArgumentError("variogram model must be stationary"))
  end

  dom = setup.domain
  data = setup.geotable
  grid, _ = unview(dom)
  dims = size(grid)
  nelms = nelements(grid)
  center = CartesianIndex(dims .÷ 2)
  cindex = LinearIndices(dims)[center]

  # number of threads in FFTW
  FFTW.set_num_threads(setu.threads)

  pairs = map(setup.vartypes, setup.varnames) do V, var
    # compute covariances between centroid and all points
    𝒟c = [centroid(grid, cindex)]
    𝒟p = [centroid(grid, i) for i in 1:nelms]
    cs = sill(γ) .- Variography.pairwise(γ, 𝒟c, 𝒟p)
    C = reshape(cs, dims)

    # move to frequency domain
    F = sqrt.(abs.(fft(fftshift(C))))
    F[1] = zero(V) # set reference level

    # perform Kriging in case of conditional simulation
    z̄, krig, searcher, dinds = nothing, nothing, nothing, nothing
    if !isnothing(data)
      ddom = domain(data)

      # retrieve process paramaters
      (; minneighbors, maxneighbors, neighborhood, distance) = process

      minneighbors, maxneighbors = fixneighlimts(ddom, minneighbors, maxneighbors)

      # determine search method
      searcher = getsearcher(ddom, maxneighbors, distance, neighborhood)

      pset = PointSet(centroid(ddom, i) for i in 1:nelements(ddom))
      kdata = georef(values(data), pset)

      # estimate conditional mean
      krig = Kriging(γ, μ)
      path = LinearPath()
      z̄ = predictvar(krig, dom, kdata, var, path, searcher, minneighbors, maxneighbors)

      # find data locations in problem domain
      lsearcher = KNearestSearch(dom, 1)
      found = [search(centroid(ddom, i), lsearcher) for i in 1:nelements(ddom)]
      dinds = unique(first.(found))
    end

    # save preprocessed inputs for variable
    var => (; F, z̄, krig, searcher, minneighbors, maxneighbors, dinds)
  end

  Dict(pairs)
end

function randsingle(rng::AbstractRNG, process::SEQ, setup::RandSetup, prep)
  # retrieve domain info
  dom = setup.domain
  grid, inds = unview(dom)
  dims = size(grid)

  # retrive variogram model and mean
  γ = process.variogram
  μ = process.mean

  varreal = map(setup.vartypes, setup.varnames) do V, var
    # unpack preprocessed parameters
    (; F, z̄, krig, searcher, minneighbors, maxneighbors, dinds) = prep[var]

    # perturbation in frequency domain
    P = F .* exp.(im .* angle.(fft(rand(rng, V, dims))))

    # move back to time domain
    Z = real(ifft(P))

    # adjust mean and variance
    σ² = Statistics.var(Z, mean=zero(V))
    Z .= √(sill(γ) / σ²) .* Z .+ μ

    # unconditional realization
    zᵤ = Z[inds]

    # perform conditioning if necessary
    z = if isnothing(krig)
      zᵤ # we are all set
    else
      # view realization at data locations
      dtab = (; var => view(zᵤ, dinds))
      ddom = view(dom, dinds)

      # solve estimation problem
      kdom = PointSet(centroid(ddom, i) for i in 1:nelements(ddom))
      kdata = georef(dtab, kdom)
      path = LinearPath()
      z̄ᵤ = predictvar(krig, kdom, kdata, var, path, searcher, minneighbors, maxneighbors)

      # add residual field
      z̄ .+ (zᵤ .- z̄ᵤ)
    end

    var => z
  end

  Dict(varreal)
end

#-----------
# UTILITIES
#-----------

function fixneighlimts(domain, min, max)
  nobs = nelements(domain)
  if max > nobs || max < 1
    @warn "Invalid maximum number of neighbors. Adjusting to $nobs..."
    max = nobs
  end

  if min > max || min < 1
    @warn "Invalid minimum number of neighbors. Adjusting to 1..."
    min = 1
  end

  min, max
end

function getsearcher(domain, maxneighbors, distance, neighborhood)
  if isnothing(neighborhood)
    # nearest neighbor search with a metric
    KNearestSearch(domain, maxneighbors; metric=distance)
  else
    # neighbor search with ball neighborhood
    KBallSearch(domain, maxneighbors, neighborhood)
  end
end

function predictvar(model, domain, data, var, path, searcher, minneighbors, maxneighbors)
  # pre-allocate memory for neighbors
  neighbors = Vector{Int}(undef, maxneighbors)

  # prediction order
  inds = traverse(domain, path)

  map(inds) do ind
    # centroid of estimation
    center = centroid(domain, ind)

    # find neighbors with data
    nneigh = search!(neighbors, center, searcher)

    # predict if enough neighbors
    if nneigh ≥ minneighbors
      # final set of neighbors
      ninds = view(neighbors, 1:nneigh)

      # view neighborhood with data
      samples = view(data, ninds)

      # fit model to data
      fmodel = fit(model, samples)

      # save prediction
      predict(fmodel, var, center)
    else # missing prediction
      missing
    end
  end
end
