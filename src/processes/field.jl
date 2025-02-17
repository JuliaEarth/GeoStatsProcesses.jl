# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    FieldProcess <: GeoStatsProcess

A field process (or random field) is defined over a fixed
geospatial domain. The realizations of the process are stored
in an ensemble, which is an indexable collection of geotables.
"""
abstract type FieldProcess <: GeoStatsProcess end

#-----------------
# IMPLEMENTATIONS
#-----------------

include("field/gaussian.jl")
include("field/lindgren.jl")
include("field/quilting.jl")
include("field/turing.jl")
include("field/strata.jl")
