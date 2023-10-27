# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    ThinningMethod

A method for thinning spatial point processes and patterns.
"""
abstract type ThinningMethod end

"""
    thin(process, method)

Thin spatial point `process` with thinning `method`.
"""
function thin end

#-----------------
# IMPLEMENTATIONS
#-----------------

include("thinning/random.jl")
