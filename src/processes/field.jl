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

"""
    iscontinuous(process::FieldProcess)

Tells whether or not the field `process` is continuous.
"""
iscontinuous(::FieldProcess) = false

"""
    isanalytical(process::FieldProcess)

Tells whether or not the field `process` is analytical,
i.e., has closed-form expressions for conditional mean, etc.
"""
isanalytical(::FieldProcess) = false

"""
    defaultschema(process::FieldProcess)

Default schema of realizations of field `process`.
"""
defaultschema(::FieldProcess) = Tables.Schema((:field,), (Float64,))

"""
    initialize(process::FieldProcess, domain, data, init)

Initialize attribute table of realization based on the field `process`,
the geospatial `domain` and the geospatial `data` using an `init`ialization
method.
"""
function initialize(process::FieldProcess, domain, data, init)
  # retrieve appropriate schema
  schema = isnothing(data) ? defaultschema(process) : dataschema(data)

  # allocate memory for realization and simulation mask
  nelm = nelements(domain)
  buff = map(T -> Vector{T}(undef, nelm), schema.types)
  bits = map(_ -> falses(nelm), schema.types)
  real = (; zip(schema.names, buff)...)
  mask = (; zip(schema.names, bits)...)

  # initialize realization and mask with data
  isnothing(data) || initialize!(real, mask, domain, data, init)

  real, mask
end

#-----------------
# IMPLEMENTATIONS
#-----------------

include("field/gaussian.jl")
include("field/indicator.jl")
include("field/lindgren.jl")
include("field/quilting.jl")
include("field/turing.jl")
include("field/strata.jl")
