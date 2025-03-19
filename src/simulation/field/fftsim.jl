# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    FFTSIM(; [options])

The FFT simulation method introduced by Oliver 1995 and Ravalec 2000.

The covariance function is perturbed in the frequency domain
after a fast Fourier transform. White noise is added to the
phase of the spectrum, and a realization is produced by an
inverse Fourier transform.

Data conditioning is currently performed with Kriging, which
accepts the following neighbor search options.

## Options

* `minneighbors` - Minimum number of neighbors (default to `1`)
* `maxneighbors` - Maximum number of neighbors (default to `26`)
* `neighborhood` - Search neighborhood (default to `nothing`)
* `distance`     - Distance used to find nearest neighbors (default to `Euclidean()`)

## References

* Oliver 1995. [Moving averages for Gaussian simulation in two and three dimensions]
  (https://link.springer.com/article/10.1007/BF02091660)

* Ravalec et al. 2000. [The FFT Moving Average (FFT-MA) Generator: An Efficient
  Numerical Method for Generating and Conditioning Gaussian Simulations]
  (https://link.springer.com/article/10.1023/A:1007542406333)

* Figueiredo et al. 2019. [Revisited Formulation and Applications of FFT Moving Average]
  (https://link.springer.com/article/10.1007/s11004-019-09826-4)

### Notes

The method is limited to simulations on regular grids, and care must be
taken to make sure that the correlation length is small enough compared to
the grid size. As a general rule of thumb, avoid correlation lengths greater
than 1/3 of the grid.

Visual artifacts can appear near the boundaries of the grid if the correlation
length is large compared to the grid itself.
"""
@kwdef struct FFTSIM{N,D} <: FieldSimulationMethod
  minneighbors::Int = 1
  maxneighbors::Int = 26
  neighborhood::N = nothing
  distance::D = Euclidean()
end

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
  Z = real(ifft(P)) * unit(μ)

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
