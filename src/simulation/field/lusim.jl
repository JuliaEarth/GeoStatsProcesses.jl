# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    LUSIM()

The LU simulation method introduced by Alabert 1987 and Davis 1987.

The full covariance matrix is built to include all locations
of the data, and samples from the multivariate Gaussian are
drawn via a lower-upper (LU) matrix factorization.

## References

* Alabert 1987. [The practice of fast conditional simulations
  through the LU decomposition of the covariance matrix]
  (https://link.springer.com/article/10.1007/BF00897191)

* Davis 1987. [Production of conditional simulations via the LU
  triangular decomposition of the covariance matrix]
  (https://link.springer.com/article/10.1007/BF00898189)

* Oliver 2003. [Gaussian cosimulation: modeling of the cross-covariance]
  (https://link.springer.com/article/10.1023%2FB%3AMATG.0000002984.56637.ef)

### Notes

The method is only adequate for domains with relatively small
number of elements (e.g. 100x100 grids) where it is feasible to
factorize the full covariance.

For larger domains (e.g. 3D grids), other methods are preferred
such as [`SEQSIM`](@ref) and [`FFTSIM`](@ref).
"""
struct LUSIM <: FieldSimulationMethod end

function preprocess(::AbstractRNG, process::GaussianProcess, method::LUSIM, init, domain, data)
  # process parameters
  f = process.func
  μ = process.mean

  # sanity checks
  isvalid(f) = isstationary(f) && issymmetric(f) && isbanded(f)
  if !isvalid(f)
    throw(ArgumentError("""
      LUSIM requires a geostatistical function that is stationary, symmetric and banded.
      Covariances or composite functions of covariances satisfy these properties.
    """))
  end

  # initialize realization and mask
  real, mask = randinit(process, domain, data, init)

  # variable names
  vars = keys(real)

  # sanity checks
  @assert length(vars) == nvariates(f) "incompatible number of variables for geostatistical function"
  @assert length(vars) ∈ (1, 2) "LUSIM only supports univariate and bivariate simulation"

  # number of variables
  nvars = length(vars)

  # preprocess parameters for variable
  preproc = map(1:nvars) do j
    # current variable
    var = vars[j]

    # retrieve data and simulation indices
    dinds = findall(mask[var])
    sinds = setdiff(1:nelements(domain), dinds)

    # data for variable
    z₁ = view(real[var], dinds)

    # mean for variable
    μ₁ = μ[j]

    # centroids for data and simulation locations
    ddom = [centroid(domain, i) for i in dinds]
    sdom = [centroid(domain, i) for i in sinds]

    # marginalize function into covariance for variable
    cov = nvars > 1 ? _marginalize(f, j) : f

    # covariance between simulation locations
    C₂₂ = _pairwise(cov, sdom)

    if isempty(dinds)
      d₂ = zero(eltype(z₁))
      L₂₂ = cholesky(Symmetric(C₂₂)).L
    else
      # covariance beween data locations
      C₁₁ = _pairwise(cov, ddom)
      C₁₂ = _pairwise(cov, ddom, sdom)

      L₁₁ = cholesky(Symmetric(C₁₁)).L
      B₁₂ = L₁₁ \ C₁₂
      A₂₁ = transpose(B₁₂)

      d₂ = A₂₁ * (L₁₁ \ z₁)
      L₂₂ = cholesky(Symmetric(C₂₂ - A₂₁ * B₁₂)).L
    end

    (; var, z₁, μ₁, d₂, L₂₂, dinds, sinds)
  end

  preproc
end

function randsingle(rng::AbstractRNG, process::GaussianProcess, ::LUSIM, domain, data, preproc)
  # simulate first variable
  var₁, Y₁, w₁ = _lusim(rng, preproc[1])
  cols = [var₁ => Y₁]

  # simulate second variable
  if length(preproc) > 1
    ρ = _rho(process.func)
    var₂, Y₂, _ = _lusim(rng, preproc[2], ρ, w₁)
    push!(cols, var₂ => Y₂)
  end

  (; cols...)
end

#------------------
# HELPER FUNCTIONS
#------------------

function _marginalize(cov, j)
  cₒ, cs, covs = structures(cov)
  cnug = cₒ[j, j]
  csum = sum(c[j, j] * cov for (c, cov) in zip(cs, covs))
  iszero(cnug) ? csum : NuggetEffect(cnug) + csum
end

function _rho(cov)
  c₁₂ = cov(0)[1, 2]
  s₁, s₂ = diag(sill(cov))
  c₁₂ / √(s₁ * s₂)
end

function _lusim(rng, preprocⱼ, ρ=nothing, w₁=nothing)
  # unpack parameters
  var, z₁, μ₁, d₂, L₂₂, dinds, sinds = preprocⱼ

  # number of elements in simulation domain
  n = length(dinds) + length(sinds)

  # allocate memory for result
  y = Vector{eltype(z₁)}(undef, n)

  # conditional simulation
  w₂ = randn(rng, size(L₂₂, 2))
  if isnothing(ρ)
    y₂ = d₂ .+ L₂₂ * w₂
  else
    y₂ = d₂ .+ L₂₂ * (ρ * w₁ + √(1 - ρ^2) * w₂)
  end

  # hard data and simulated values
  y[dinds] = z₁
  y[sinds] = y₂

  # adjust mean in case of unconditional simulation
  isempty(dinds) && (y .+= μ₁)

  var, y, w₂
end
