# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function preprocess(::AbstractRNG, process::GaussianProcess, method::LUSIM, domain, data)
  # process parameters
  f = process.func
  μ = process.mean

  # method options
  fact = method.factorization
  init = method.init

  # variable names
  vars = if isnothing(data)
    ntuple(i -> Symbol(:Z, i), nvariates(f))
  else
    data |> values |> Tables.columns |> Tables.columnnames
  end

  # sanity checks
  if length(vars) != nvariates(f)
    throw(ArgumentError("incompatible number of variables for geostatistical function"))
  end

  # initialize buffers for realizations and simulation mask
  buff, mask = initbuff(domain, vars, init, data=data)

  # preprocess parameters for individual variables
  pairs = map(enumerate(varnames)) do (i, var)
    # get variable specific parameters
    f = _getparam(process.func, i)
    μ = _getparam(process.mean, i)

    # check stationarity
    @assert isstationary(f) "geostatistical function must be stationary"

    # retrieve data locations and data values in domain
    dlocs = findall(mask[var])
    z₁ = view(buff[var], dlocs)

    # retrieve simulation locations
    slocs = setdiff(1:nelements(domain), dlocs)

    # create views of the domain
    𝒟d = [centroid(domain, i) for i in dlocs]
    𝒟s = [centroid(domain, i) for i in slocs]

    # covariance between simulation locations
    C₂₂ = _pairwise(f, 𝒟s)

    if isempty(dlocs)
      d₂ = zero(eltype(z₁))
      L₂₂ = fact(Symmetric(C₂₂)).L
    else
      # covariance beween data locations
      C₁₁ = _pairwise(f, 𝒟d)
      C₁₂ = _pairwise(f, 𝒟d, 𝒟s)

      L₁₁ = fact(Symmetric(C₁₁)).L
      B₁₂ = L₁₁ \ C₁₂
      A₂₁ = B₁₂'

      d₂ = A₂₁ * (L₁₁ \ z₁)
      L₂₂ = fact(Symmetric(C₂₂ - A₂₁ * B₁₂)).L
    end

    # save preprocessed parameters for variable
    var => (z₁, d₂, L₂₂, μ, dlocs, slocs)
  end

  Dict(pairs)
end

function randsingle(rng::AbstractRNG, ::GaussianProcess, method::LUSIM, domain, data)
  # list of variable names
  vars = setup.varnames

  # simulate first variable
  v₁ = first(vars)
  Y₁, w₁ = _lusim(rng, prep[v₁])
  pairs = [v₁ => Y₁]

  # simulate second variable
  if length(vars) == 2
    ρ = method.correlation
    v₂ = last(vars)
    Y₂, _ = _lusim(rng, prep[v₂], ρ, w₁)
    push!(pairs, v₂ => Y₂)
  end

  (; pairs...)
end

#-----------
# UTILITIES
#-----------

function _lusim(rng, params, ρ=nothing, w₁=nothing)
  # unpack parameters
  z₁, d₂, L₂₂, μ, dlocs, slocs = params

  # number of points in domain
  npts = length(dlocs) + length(slocs)

  # allocate memory for result
  y = Vector{eltype(z₁)}(undef, npts)

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
