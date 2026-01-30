# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    expectedvalue(process::FieldProcess, domain::Domain; data=nothing, kwargs...)

Compute the expected value of the field `process` over the given `domain`,
conditioned on the `data` values if provided. Optionally, forward `kwargs`.
"""
function expectedvalue(process::FieldProcess, domain::Domain; data=nothing, kwargs...)
  if iscontinuous(process)
    mean(process, domain; data, kwargs...)
  else
    error("not implemented")
  end
end

"""
    mean(process::FieldProcess, domain::Domain; data=nothing, kwargs...)

Compute the mean of the continuous field `process` over the given `domain`,
conditioned on the `data` values if provided. Optionally, forward `kwargs`.
"""
mean(process::FieldProcess, domain::Domain; data=nothing, kwargs...) =
  isnothing(data) ? priormean(process, domain; kwargs...) : posteriormean(process, domain, data; kwargs...)

# ----------------
# IMPLEMENTATIONS
# ----------------

include("field/gaussian.jl")
include("field/lindgren.jl")
