# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function preprocess(::AbstractRNG, process::GaussianProcess, method::LUSIM, init, domain, data)
  # process parameters
  f = process.func
  μ = process.mean

  # method options
  factorization = method.factorization

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

    # retrieve data and simulation locations
    dlocs = findall(mask[var])
    slocs = setdiff(1:nelements(domain), dlocs)

    # data for variable
    z₁ = view(real[var], dlocs)

    # centroids for data and simulation locations
    ddom = [centroid(domain, i) for i in dlocs]
    sdom = [centroid(domain, i) for i in slocs]

    # marginalize function into covariance for variable
    cov = nvars > 1 ? _marginalize(f, j) : f

    # covariance between simulation locations
    C₂₂ = _pairwise(cov, sdom)

    if isempty(dlocs)
      d₂ = zero(eltype(z₁))
      L₂₂ = factorization(Symmetric(C₂₂)).L
    else
      # covariance beween data locations
      C₁₁ = _pairwise(cov, ddom)
      C₁₂ = _pairwise(cov, ddom, sdom)

      L₁₁ = factorization(Symmetric(C₁₁)).L
      B₁₂ = L₁₁ \ C₁₂
      A₂₁ = transpose(B₁₂)

      d₂ = A₂₁ * (L₁₁ \ z₁)
      L₂₂ = factorization(Symmetric(C₂₂ - A₂₁ * B₁₂)).L
    end

    (var, (; z₁, d₂, L₂₂, μ, dlocs, slocs))
  end

  preproc
end

function randsingle(rng::AbstractRNG, process::GaussianProcess, ::LUSIM, domain, data, preproc)
  # unpack preprocessing results
  var₁, params₁ = first(preproc)
  var₂, params₂ = last(preproc)

  # simulate first variable
  Y₁, w₁ = _lusim(rng, params₁)
  cols = [var₁ => Y₁]

  # simulate second variable
  if length(preproc) > 1
    ρ = _rho(process.func)
    Y₂, _ = _lusim(rng, params₂, ρ, w₁)
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

function _lusim(rng, params, ρ=nothing, w₁=nothing)
  # unpack parameters
  z₁, d₂, L₂₂, μ, dlocs, slocs = params

  # number of elements in simulation domain
  n = length(dlocs) + length(slocs)

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
  y[dlocs] = z₁
  y[slocs] = y₂

  # adjust mean in case of unconditional simulation
  isempty(dlocs) && (y .+= μ)

  y, w₂
end
