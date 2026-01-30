# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function priormean(process::LindgrenProcess, domain; kwargs...)
  # retrieve parameters
  σ² = process.sill
  V = typeof(√σ²)

  # fill domain with mean values
  vars = (:field,)
  vals = (fill(zero(V), nelements(domain)),)
  # georeference mean values
  georef((; zip(vars, vals)...), domain)
end

function posteriormean(process::LindgrenProcess, domain, data; kwargs...)
  error("not implemented")
end
