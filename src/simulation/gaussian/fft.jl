# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function preprocess(::AbstractRNG, process::GaussianProcess, method::FFTSIM, domain, data)
  # retrive function and mean
  f = process.func
  μ = process.mean

  # check stationarity
  if !isstationary(f)
    throw(ArgumentError("geostatistical function must be stationary"))
  end

  dom = setup.domain
  data = setup.geotable
  grid = parent(dom)
  dims = size(grid)
  center = CartesianIndex(dims .÷ 2)
  cindex = LinearIndices(dims)[center]

  # number of threads in FFTW
  FFTW.set_num_threads(setup.threads)

  # perform Kriging in case of conditional simulation
  pred = if isnothing(data)
    nothing
  else
    (; minneighbors, maxneighbors, neighborhood, distance) = method
    pred = GeoStatsModels.fitpredict(Kriging(f, μ), data, dom; minneighbors, maxneighbors, neighborhood, distance)
  end

  pairs = map(setup.vartypes, setup.varnames) do V, var
    # compute covariances between centroid and all points
    𝒟c = [centroid(grid, cindex)]
    𝒟p = [centroid(grid, i) for i in 1:nelements(grid)]
    cs = _pairwise(f, 𝒟c, 𝒟p)
    C = reshape(cs, dims)

    # move to frequency domain
    F = sqrt.(abs.(fft(fftshift(C))))
    F[1] = zero(V) # set reference level

    # get variable prediction and data locations if necessary
    z̄, dinds = nothing, nothing
    if !isnothing(pred)
      z̄ = pred[:, var]
      # find data locations in target domain
      ddom = domain(data)
      searcher = KNearestSearch(dom, 1)
      found = [search(centroid(ddom, i), searcher) for i in 1:nelements(ddom)]
      dinds = unique(first.(found))
    end

    # save preprocessed inputs for variable
    var => (; F, z̄, dinds)
  end

  Dict(pairs)
end

function randsingle(rng::AbstractRNG, process::GaussianProcess, method::FFTSIM, domain, data)
  # retrieve domain info
  dom = setup.domain
  data = setup.geotable
  grid = parent(dom)
  inds = parentindices(dom)
  dims = size(grid)

  # function and mean
  f = process.func
  μ = process.mean

  pairs = map(setup.vartypes, setup.varnames) do V, var
    # unpack preprocessed parameters
    (; F, z̄, dinds) = prep[var]

    # perturbation in frequency domain
    P = F .* exp.(im .* angle.(fft(rand(rng, V, dims))))

    # move back to time domain
    Z = real(ifft(P))

    # adjust mean and variance
    σ² = Statistics.var(Z, mean=zero(V))
    Z .= √(sill(f) / σ²) .* Z .+ μ

    # unconditional realization
    zᵤ = Z[inds]

    # perform conditioning if necessary
    z = if isnothing(data)
      zᵤ # we are all set
    else
      # view realization at data locations
      ktab = (; var => view(zᵤ, dinds))
      kdom = view(dom, dinds)
      kdata = georef(ktab, kdom)

      # perform Kriging prediction
      (; minneighbors, maxneighbors, neighborhood, distance) = method
      pred = GeoStatsModels.fitpredict(Kriging(f, μ), kdata, dom; minneighbors, maxneighbors, neighborhood, distance)
      z̄ᵤ = pred[:, var]

      # add residual field
      z̄ .+ (zᵤ .- z̄ᵤ)
    end

    var => z
  end

  (; pairs...)
end
