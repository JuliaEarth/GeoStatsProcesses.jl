# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function priormean(process::LindgrenProcess, domain; kwargs...)
  # process parameters
  σ² = process.sill
  V = typeof(√σ²)

  # fill domain with mean values
  vars = (:field,)
  vals = (fill(zero(V), nelements(domain)),)
  # georeference mean values
  georef((; zip(vars, vals)...), domain)
end

function posteriormean(process::LindgrenProcess, domain, data; init, kwargs...)
  # process parameters
  σ² = process.sill

  # initialize mean and mask at vertices
  mean, mask = initialize(process, PointSet(vertices(domain)), data, init)

  # multivariate expectation is not supported
  @assert length(keys(mean)) == 1 "Lindgren's process does not support multivariate expectation"

  # retrieve variable name
  var = first(keys(mean))

  # retrieve data and expectation locations
  i₁ = findall(mask[var])
  i₂ = setdiff(1:nvertices(domain), i₁)

  # Matérn precision matrix
  Q = precisionmatrix(process, domain)

  # perform interpolation
  z = ustrip.(mean[var])
  z₁ = view(z, i₁)
  z₂ = view(z, i₂)
  z₂ .= -Q[i₂,i₂] \ (Q[i₂,i₁] * z₁)

  # vertex table
  vtable = (; var => z * √unit(σ²))

  # change of support
  vdata = GeoTable(domain; vtable)
  integrate(vdata, var)
end
