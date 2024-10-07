# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    FFTMethod(; [paramaters])

The FFT Gaussian process method introduced by Gutjahr 1997.
The covariance function is perturbed in the frequency domain
after a fast Fourier transform. White noise is added to the
phase of the spectrum, and a realization is produced by an
inverse Fourier transform.

## Parameters
* `minneighbors` - Minimum number of neighbors (default to `1`)
* `maxneighbors` - Maximum number of neighbors (default to `10`)
* `neighborhood` - Search neighborhood (default to `nothing`)
* `distance`     - Distance used to find nearest neighbors (default to `Euclidean()`)

## References

* Gutjahr 1997. [General joint conditional simulations using a fast
  Fourier transform method](https://link.springer.com/article/10.1007/BF02769641)

* Gómez-Hernández, J. & Srivastava, R. 2021. [One Step at a Time: The Origins
  of Sequential Simulation and Beyond](https://link.springer.com/article/10.1007/s11004-021-09926-0)

### Notes

* The method is limited to simulations on Cartesian grids, and care must be
  taken to make sure that the correlation length is small enough compared to
  the grid size.

* As a general rule of thumb, avoid correlation lengths greater than 1/3 of
  the grid.

* The method is extremely fast, and can be used to generate large 3D realizations.
"""
@kwdef struct FFTMethod{N,D} <: RandMethod
  minneighbors::Int = 1
  maxneighbors::Int = 10
  neighborhood::N = nothing
  distance::D = Euclidean()
end

function randprep(::AbstractRNG, process::GaussianProcess, method::FFTMethod, setup::RandSetup)
  # retrive variogram model and mean
  γ = process.variogram
  μ = process.mean

  # check stationarity
  if !isstationary(γ)
    throw(ArgumentError("variogram model must be stationary"))
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
    pred = GeoStatsModels.fitpredict(Kriging(γ, μ), data, dom; minneighbors, maxneighbors, neighborhood, distance)
  end

  pairs = map(setup.vartypes, setup.varnames) do V, var
    # compute covariances between centroid and all points
    𝒟c = [centroid(grid, cindex)]
    𝒟p = [centroid(grid, i) for i in 1:nelements(grid)]
    cs = sill(γ) .- GeoStatsFunctions.pairwise(γ, 𝒟c, 𝒟p)
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

function randsingle(rng::AbstractRNG, process::GaussianProcess, method::FFTMethod, setup::RandSetup, prep)
  # retrieve domain info
  dom = setup.domain
  data = setup.geotable
  grid = parent(dom)
  inds = parentindices(dom)
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
      (; minneighbors, maxneighbors, neighborhood, distance) = method
      pred = GeoStatsModels.fitpredict(Kriging(γ, μ), kdata, dom; minneighbors, maxneighbors, neighborhood, distance)
      z̄ᵤ = pred[:, var]

      # add residual field
      z̄ .+ (zᵤ .- z̄ᵤ)
    end

    var => z
  end

  (; varreal...)
end
