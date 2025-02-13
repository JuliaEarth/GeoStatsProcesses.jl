# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    InitMethod

A method to initialize buffers in geostatistical solvers.
"""
abstract type InitMethod end

"""
    initbuff(domain, vars, method; [data])

Initialize buffers for all variables `vars` with given `method`
based on the location of elements in `domain` and (optionally) `data`.
"""
function initbuff(domain, vars, method::InitMethod; data=nothing)
  nelm = nelements(domain)
  buff = Dict(var => Vector{V}(undef, nelm) for (var, V) in pairs(vars))
  mask = Dict(var => falses(nelm) for (var, V) in pairs(vars))

  if !isnothing(data)
    ivars = keys(vars)
    dvars = setdiff(propertynames(data), [:geometry])
    preproc = preprocess(data, domain, method)
    for var in ivars ∩ dvars
      cols = Tables.columns(values(data))
      vals = Tables.getcolumn(cols, var)
      initbuff!(buff[var], mask[var], vals, method, preproc)
    end
  end

  buff, mask
end

# ----------------
# IMPLEMENTATIONS
# ----------------

include("initbuff/nearest.jl")
include("initbuff/explicit.jl")
