# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    expectedvalue(process::FieldProcess, domain::Domain; data=nothing, init=NearestInit(), kwargs...)

Compute the expected value of the field `process` over the given `domain`,
conditioned on the `data` values if provided.

By default, the `data` is assigned to the nearest geometry of the `domain`,
but other `init`ialization methods are available. Additional `kwargs` are
forwarded to specialized implementations.
"""
function expectedvalue(process::FieldProcess, domain::Domain; kwargs...)
  if iscontinuous(process)
    mean(process, domain; kwargs...)
  else
    error("not implemented")
  end
end

"""
    mean(process::FieldProcess, domain::Domain; data=nothing, init=NearestInit(), kwargs...)

Compute the mean of the continuous field `process` over the given `domain`,
conditioned on the `data` values if provided.

By default, the `data` is assigned to the nearest geometry of the `domain`,
but other `init`ialization methods are available. Additional `kwargs` are
forwarded to specialized implementations.
"""
mean(process::FieldProcess, domain::Domain; data=nothing, init=NearestInit(), kwargs...) =
  isnothing(data) ? priormean(process, domain; kwargs...) : posteriormean(process, domain, data; init, kwargs...)

# ----------------
# IMPLEMENTATIONS
# ----------------

include("field/gaussian.jl")
include("field/lindgren.jl")
