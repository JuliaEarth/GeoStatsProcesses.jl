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
  center = CartesianIndex(dims .÷ 2)
  cindex = LinearIndices(dims)[center]

  # number of threads in FFTW
  FFTW.set_num_threads(setup.threads)

  # perform Kriging in case of conditional simulation
  pred = nothing
  if !isnothing(data)
    krig = Kriging(γ, μ)
    (; minneighbors, maxneighbors, neighborhood, distance) = process
    pred = fitpredict(krig, data, dom; minneighbors, maxneighbors, neighborhood, distance)
  end

  pairs = map(setup.vartypes, setup.varnames) do V, var
    # compute covariances between centroid and all points
    𝒟c = [centroid(grid, cindex)]
    𝒟p = [centroid(grid, i) for i in 1:nelements(grid)]
    cs = sill(γ) .- Variography.pairwise(γ, 𝒟c, 𝒟p)
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

function randsingle(rng::AbstractRNG, process::FFTGP, setup::RandSetup, prep)
  # retrieve domain info
  dom = setup.domain
  data = setup.geotable
  grid, inds = unview(dom)
  dims = size(grid)

  # retrive variogram model and mean
  γ = process.variogram
  μ = process.mean

  varreal = map(setup.vartypes, setup.varnames) do V, var
    # unpack preprocessed parameters
    (; F, z̄, dinds) = prep[var]

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
    z = if isnothing(data)
      zᵤ # we are all set
    else
      # view realization at data locations
      ktab = (; var => view(zᵤ, dinds))
      kdom = view(dom, dinds)
      kdata = georef(ktab, kdom)

      # perform Kriging prediction
      krig = Kriging(γ, μ)
      (; minneighbors, maxneighbors, neighborhood, distance) = process
      pred = fitpredict(krig, kdata, kdom; minneighbors, maxneighbors, neighborhood, distance)
      z̄ᵤ = pred[:, var]

      # add residual field
      z̄ .+ (zᵤ .- z̄ᵤ)
    end

    var => z
  end

  Dict(varreal)
end
