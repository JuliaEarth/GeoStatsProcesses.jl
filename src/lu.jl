# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

@kwdef struct LUGP{V,T,F,C,I} <: GeoStatsProcess
  variogram::V = GaussianVariogram()
  mean::T = nothing
  factorization::F = cholesky
  correlation::C = 0.0
  init::I = NearestInit()
end

function randprep(::AbstractRNG, process::LUGP, setup::RandSetup)
  # retrieve setup paramaters
  (; domain, geotable, varnames, vartypes) = setup

  # check the number of variables
  nvars = length(varnames)
  @assert nvars ∈ (1, 2) "only 1 or 2 variables can be simulated simultaneously"

  # check process paramaters
  _checkparam(process.variogram, nvars)
  _checkparam(process.mean, nvars)

  # retrieve process parameters
  fact = process.factorization
  init = process.init

  # initialize buffers for realizations and simulation mask
  vars = Dict(zip(varnames, vartypes))
  buff, mask = initbuff(domain, vars, init, data=geotable)

  # preprocess parameters for individual variables
  pairs = map(enumerate(varnames)) do (i, var)
    # get variable specific parameters
    γ = _getparam(process.variogram, i)
    vmean = _getparam(process.mean, i)

    # check stationarity
    @assert isstationary(γ) "variogram model must be stationary"

    # retrieve data locations and data values in domain
    dlocs = findall(mask[var])
    z₁ = view(buff[var], dlocs)

    # retrieve simulation locations
    slocs = setdiff(1:nelements(domain), dlocs)

    # create views of the domain
    𝒟d = [centroid(domain, i) for i in dlocs]
    𝒟s = [centroid(domain, i) for i in slocs]

    # covariance between simulation locations
    C₂₂ = sill(γ) .- Variography.pairwise(γ, 𝒟s)

    if isempty(dlocs)
      d₂ = zero(eltype(z₁))
      L₂₂ = fact(Symmetric(C₂₂)).L
    else
      # covariance beween data locations
      C₁₁ = sill(γ) .- Variography.pairwise(γ, 𝒟d)
      C₁₂ = sill(γ) .- Variography.pairwise(γ, 𝒟d, 𝒟s)

      L₁₁ = fact(Symmetric(C₁₁)).L
      B₁₂ = L₁₁ \ C₁₂
      A₂₁ = B₁₂'

      d₂ = A₂₁ * (L₁₁ \ z₁)
      L₂₂ = fact(Symmetric(C₂₂ - A₂₁ * B₁₂)).L
    end

    if !isnothing(vmean) && !isempty(dlocs)
      @warn "mean can only be specified in unconditional simulation"
    end

    # mean for unconditional simulation
    μ = isnothing(vmean) ? zero(eltype(z₁)) : vmean

    # save preprocessed parameters for variable
    var => (z₁, d₂, L₂₂, μ, dlocs, slocs)
  end

  Dict(pairs)
end

function randsingle(rng::AbstractRNG, process::LUGP, setup::RandSetup, prep)
  # list of variable names
  vars = setup.varnames

  # simulate first variable
  v₁ = first(vars)
  Y₁, w₁ = _lusim(rng, prep[v₁])
  varreal = Dict(v₁ => Y₁)

  # simulate second variable
  if length(vars) == 2
    ρ = process.correlation
    v₂ = last(vars)
    Y₂, _ = _lusim(rng, prep[v₂], ρ, w₁)
    push!(varreal, v₂ => Y₂)
  end

  varreal
end

#-----------
# UTILITIES
#-----------

function _checkparam(param, nvars)
  if param isa Tuple
    @assert length(param) == nvars "the number of parameters must be equal to the number of variables"
  end
end

_getparam(param, i) = param
_getparam(params::Tuple, i) = params[i]

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