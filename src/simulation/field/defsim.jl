# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function preprocess(::AbstractRNG, process::LindgrenProcess, ::Nothing, init, domain, data)
  # process parameters
  σ² = process.sill

  # initialize realization and mask at vertices
  real, mask = initialize(process, PointSet(vertices(domain)), data, init)

  # multivariate simulation is not supported
  @assert length(keys(real)) == 1 "Lindgren's process does not support multivariate simulation"

  # retrieve variable name
  var = first(keys(real))

  # retrieve data and simulation locations
  i₁ = findall(mask[var])
  i₂ = setdiff(1:nvertices(domain), i₁)

  # Matérn precision matrix
  Q = precisionmatrix(process, domain)

  # matrix factorization
  F = cholesky(Symmetric(Q))

  # interpolate at simulation locations if necessary
  z̄ = if isempty(i₁)
    nothing
  else
    z = ustrip.(real[var])
    z₁ = view(z, i₁)
    z₂ = view(z, i₂)
    z₂ .= -Q[i₂,i₂] \ (Q[i₂,i₁] * z₁)
    z
  end

  (; var, Q, F, σ², i₁, i₂, z̄)
end

function randsingle(rng::AbstractRNG, ::LindgrenProcess, ::Nothing, domain, data, preproc)
  # unpack preprocessing results
  (; var, Q, F, σ², i₁, i₂, z̄) = preproc

  # unconditional realization at vertices
  w = randn(rng, eltype(F), size(F, 1))
  zᵤ = F \ w

  # adjust variance
  s² = Statistics.var(zᵤ, mean=zero(eltype(zᵤ)))
  zᵤ .= √(ustrip(σ²) / s²) .* zᵤ

  # perform conditioning if necessary
  z = if isempty(i₁)
    zᵤ # we are all set
  else
    # view realization at data locations
    zᵤ₁ = view(zᵤ, i₁)

    # interpolate at simulation locations
    zᵤ₂ = -Q[i₂,i₂] \ (Q[i₂,i₁] * zᵤ₁)

    # merge the above results
    z̄ᵤ = similar(zᵤ)
    z̄ᵤ[i₁] .= zᵤ₁
    z̄ᵤ[i₂] .= zᵤ₂

    # add residual field
    z̄ .+ (zᵤ .- z̄ᵤ)
  end

  # vertex table
  vtable = (; var => z * √unit(σ²))

  # change of support
  edata = average(GeoTable(domain; vtable))

  # return attribute table
  values(edata)
end
