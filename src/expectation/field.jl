# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    expectedvalue(process::FieldProcess, domain::Domain; data=nothing)

Compute the expected value of the field `process` over the given `domain`,
conditioned on the `data` values if provided.
"""
function expectedvalue(process::FieldProcess, domain::Domain; data=nothing)
  if iscontinuous(process)
    mean(process, domain; data)
  else
    error("not implemented")
  end
end

"""
    mean(process::FieldProcess, domain::Domain; data=nothing)

Compute the prior or posterior mean of the continuous field `process`
over the given `domain`. If `data` is provided, compute the posterior
mean conditioned on the `data` values; otherwise, compute the prior mean.
"""
mean(process::FieldProcess, domain::Domain; data=nothing) =
  isnothing(data) ? priormean(process, domain) : posteriormean(process, domain, data)

# ----------------
# IMPLEMENTATIONS
# ----------------

include("field/gaussian.jl")
include("field/lindgren.jl")
