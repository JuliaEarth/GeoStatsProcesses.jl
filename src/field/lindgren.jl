# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    LindgrenProcess(range=1.0, sill=1.0; init=NearestInit())

Lindgren process with given `range` (correlation length)
and `sill` (total variance) as described in Lindgren 2011.
Optionally, specify the data initialization method `init`.

The process relies relies on a discretization of the Laplace-Beltrami 
operator on meshes and is adequate for highly curved domains (e.g. surfaces).

## References

* Lindgren et al. 2011. [An explicit link between Gaussian fields and
  Gaussian Markov random fields: the stochastic partial differential
  equation approach](https://rss.onlinelibrary.wiley.com/doi/10.1111/j.1467-9868.2011.00777.x)
"""
struct LindgrenProcess{ℒ<:Len,V,I} <: FieldProcess
  range::ℒ
  sill::V
  init::I
  LindgrenProcess(range::ℒ, sill::V, init::I) where {ℒ<:Len,V,I} = new{float(ℒ),float(V),I}(range, sill, init)
end

LindgrenProcess(range=1.0u"m", sill=1.0; init=NearestInit()) = LindgrenProcess(addunit(range, u"m"), sill, init)

function randprep(::AbstractRNG, process::LindgrenProcess, ::DefaultRandMethod, setup::RandSetup)
  # retrieve setup paramaters
  (; domain, geotable, varnames, vartypes) = setup

  # retrieve sill and range
  𝓁 = process.range
  σ = process.sill

  # retrieve initialization method
  init = process.init

  # sanity checks
  @assert domain isa Mesh "domain must be a `Mesh`"
  @assert 𝓁 > zero(𝓁) "range must be positive"
  @assert σ > zero(σ) "sill must be positive"

  # Laplace-Beltrami operator
  W = laplacematrix(domain)
  M = measurematrix(domain)
  Δ = inv(M) * W

  # retrieve parametric dimension
  d = paramdim(domain)

  # LHS of SPDE (κ² - Δ)Z = τW with Δ = M⁻¹W
  α = 2
  ν = α - d / 2
  κ = 1 / 𝓁
  A = κ^2 * I - Δ

  # Matérn precision matrix
  τ² = σ^2 * κ^(2ν) * (4π)^(d / 2) * gamma(α) / gamma(ν)
  Q = ustrip.(A'A / τ²)

  # factorization
  F = cholesky(Array(Q))
  L = inv(Array(F.U))

  # initialize buffers for realizations and simulation mask
  pset = PointSet(vertices(domain))
  vars = Dict(zip(varnames, vartypes))
  buff, mask = initbuff(pset, vars, init, data=geotable)

  # result of preprocessing
  pairs = map(varnames) do var
    # retrieve buffer and mask for variable
    z = buff[var]
    m = mask[var]

    # retrieve data locations and data values
    i₁ = findall(m)
    z₁ = view(z, i₁)

    # retrieve simulation locations
    i₂ = setdiff(1:nvertices(domain), i₁)

    # interpolate at simulation locations if necessary
    z̄ = if isempty(i₁)
      nothing
    else
      z[i₂] .= -Q[i₂,i₂] \ (Q[i₂,i₁] * z₁)
      z
    end

    # save preprocessed inputs for variable
    var => (; Q, L, i₁, i₂, z̄)
  end

  Dict(pairs)
end

function randsingle(rng::AbstractRNG, ::LindgrenProcess, ::DefaultRandMethod, setup::RandSetup, prep)
  # retrieve setup paramaters
  (; domain, geotable, varnames, vartypes) = setup

  # simulation at vertices
  varreal = map(varnames, vartypes) do var, V
    # unpack preprocessed parameters
    (; Q, L, i₁, i₂, z̄) = prep[var]

    # unconditional realization
    w = randn(rng, V, nvertices(domain))
    zᵤ = L * w

    # perform conditioning if necessary
    z = if isempty(i₁)
      zᵤ # we are all set
    else
      # view realization at data locations
      zᵤ₁ = view(zᵤ, i₁)

      # interpolate at simulation locations
      z̄ᵤ = similar(zᵤ)
      zᵤ₂ = -Q[i₂,i₂] \ (Q[i₂,i₁] * zᵤ₁)
      z̄ᵤ[i₁] .= zᵤ₁
      z̄ᵤ[i₂] .= zᵤ₂

      # add residual field
      z̄ .+ (zᵤ .- z̄ᵤ)
    end

    var => z
  end

  # vertex table
  vtable = (; varreal...)

  # change of support
  vdata = GeoTable(domain; vtable)
  edata = integrate(vdata, varnames...)

  # columns of element table
  cols = Tables.columns(values(edata))
  Dict(var => Tables.getcolumn(cols, var) for var in varnames)
end
