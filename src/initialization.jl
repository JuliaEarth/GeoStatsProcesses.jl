# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    InitMethod

A method to initialize realizations in geostatistical simulation.
"""
abstract type InitMethod end

"""
    initialize!(real, mask, domain, data, method)

Initialize `real`ization and `mask` for given `domain` and `data`
using a initialization `method`.
"""
function initialize! end

# ----------------
# IMPLEMENTATIONS
# ----------------

include("initialization/nearest.jl")
include("initialization/explicit.jl")
