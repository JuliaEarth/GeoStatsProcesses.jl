# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function preprocess(::AbstractRNG, process::GaussianProcess, method::FFTSIM, init, domain, data)
  # function and mean
  f = process.func
  μ = process.mean

  # sanity checks
  @assert isstationary(f) "geostatistical function must be stationary"

  # simulation domain
  sdom = domain

  # initialize realization and mask
  real, mask = randinit(process, sdom, data, init)

  # multivariate simulation is not supported
  @assert length(keys(real)) == 1 "FFTSIM does not support multivariate simulation"

  # retrieve variable name
  var = first(keys(real))

  # set number of threads in FFTW
  FFTW.set_num_threads(cpucores())

  # index of grid center
  grid = parent(sdom)
  dims = size(grid)
  cind = CartesianIndex(dims .÷ 2)
  lind = LinearIndices(dims)[cind]

  # compute covariances between centroid and all points
  cdom = [centroid(grid, lind)]
  gdom = [centroid(grid, i) for i in 1:nelements(grid)]
  covs = _pairwise(f, cdom, gdom)
  C = reshape(covs, dims)

  # move to frequency domain
  F = sqrt.(abs.(fft(fftshift(C))))
  F[1] = zero(eltype(F)) # set reference level

  # perform Kriging in case of conditional simulation
  z̄ = if isnothing(data)
    nothing
  else
    # conditional mean
    (; minneighbors, maxneighbors, neighborhood, distance) = method
    krig = GeoStatsModels.fitpredict(Kriging(f, μ), data, sdom; minneighbors, maxneighbors, neighborhood, distance)
    krig[:, var]
  end

  # save data locations
  dinds = findall(mask[var])

  (; var, F, z̄, dinds)
end

function randsingle(rng::AbstractRNG, process::GaussianProcess, method::FFTSIM, domain, data, preproc)
  # unpack preprocessing results
  (; var, F, z̄, dinds) = preproc

  # function and mean
  f = process.func
  μ = process.mean

  # simulation grid
  sdom = domain
  grid = parent(sdom)
  inds = parentindices(sdom)
  dims = size(grid)

  # perturbation in frequency domain
  w = rand(rng, eltype(F), dims)
  P = F .* exp.(im .* angle.(fft(w)))

  # move back to time domain
  Z = real(ifft(P))

  # adjust mean and variance
  σ² = Statistics.var(Z, mean=zero(eltype(Z)))
  Z .= √(sill(f) / σ²) .* Z .+ μ

  # unconditional realization
  zᵤ = Z[inds]

  # perform conditioning if necessary
  z = if isnothing(data)
    zᵤ # we are all set
  else
    # view realization at data locations
    ktab = (; var => view(zᵤ, dinds))
    kdom = view(sdom, dinds)
    kdata = georef(ktab, kdom)

    # perform Kriging with unconditional values
    (; minneighbors, maxneighbors, neighborhood, distance) = method
    krigᵤ = GeoStatsModels.fitpredict(Kriging(f, μ), kdata, sdom; minneighbors, maxneighbors, neighborhood, distance)
    z̄ᵤ = krigᵤ[:, var]

    # add residual field
    z̄ .+ (zᵤ .- z̄ᵤ)
  end

  (; var => z)
end
