# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    ThinningMethod

A method for thinning point processes and patterns.
"""
abstract type ThinningMethod end

"""
    thin(process, method)

Thin point `process` with thinning `method`.
"""
function thin end

#-----------------
# IMPLEMENTATIONS
#-----------------

include("thinning/random.jl")
