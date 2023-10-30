# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    GeoStatsProcess

Parent type of all geostatistical processes.
"""
abstract type GeoStatsProcess end

"""
    PointProcess <: GeoStatsProcess

Parent type of all point processes.
"""
abstract type PointProcess <: GeoStatsProcess end

"""
    FieldProcess <: GeoStatsProcess

Parent type of all field processes.
"""
abstract type FieldProcess <: GeoStatsProcess end

"""
    RandMethod

Parent type of all processes methods.
"""
abstract type RandMethod end

struct DefaultRandMethod <: RandMethod end

"""
    defaultmethod(process::FieldProcess, setup) -> RandMethod

Returns the default method for the corresponding field `process`.
"""
defaultmethod(::FieldProcess, setup) = DefaultRandMethod()

#-----------------
# FIELD PROCESSES
#-----------------

"""
    rand([rng], process::FieldProcess, domain, data, [method]; paramaters...)
    rand([rng], process::FieldProcess, domain, data, nreals, [method]; paramaters...)

Generate one or `nreals` realizations of the field `process` with `method`
over the `domain` with `data` and optional `paramaters`.

The `data` can be a geotable or an iterable of pairs of the form `var => T`,
where `var` is a symbol or string with the variable name and `T` is the corresponding
data type.

# Examples

```julia
julia> rand(process, domain, [:z => Float64])
julia> rand(process, domain, geotable, 3)
```
"""
Base.rand(process::FieldProcess, domain::Domain, data, method=nothing; kwargs...) =
  rand(Random.default_rng(), process, domain, data, method; kwargs...)

function Base.rand(
  rng::AbstractRNG, 
  process::FieldProcess, 
  domain::Domain, 
  data, 
  method=nothing; 
  threads=cpucores()
)
  setup = randsetup(domain, data, threads)
  rmethod = isnothing(method) ? defaultmethod(process, setup) : method
  prep = randprep(rng, process, rmethod, setup)
  real = randsingle(rng, process, rmethod, setup, prep)
  table = (; (var => real[var] for var in setup.varnames)...)
  georef(table, domain)
end

Base.rand(process::FieldProcess, domain::Domain, data, nreals::Integer, method=nothing; kwargs...) =
  rand(Random.default_rng(), process, domain, data, nreals, method; kwargs...)

function Base.rand(
  rng::AbstractRNG,
  process::FieldProcess,
  domain::Domain,
  data,
  nreals::Integer,
  method=nothing;
  pool=[myid()],
  threads=cpucores(),
  progress=true
)
  setup = randsetup(domain, data, threads)
  rmethod = isnothing(method) ? defaultmethod(process, setup) : method
  prep = randprep(rng, process, rmethod, setup)

  # pool of worker processes
  pool = CachingPool(pool)

  # simulation loop
  reals = if progress
    message = "Simulating $(join(setup.varnames, " ,", " and ")):"
    @showprogress desc = message pmap(pool, 1:nreals) do _
      randsingle(rng, process, rmethod, setup, prep)
    end
  else
    pmap(pool, 1:nreals) do _
      randsingle(rng, process, rmethod, setup, prep)
    end
  end

  # rearrange realizations
  (; varnames, vartypes) = setup
  varvects = [Vector{V}[] for V in vartypes]
  varreals = (; zip(varnames, varvects)...)
  for real in reals
    for var in varnames
      push!(varreals[var], real[var])
    end
  end

  Ensemble(domain, varreals)
end

#-----------
# UTILITIES
#-----------

struct RandSetup{D<:Domain,T}
  domain::D
  geotable::T
  varnames::Vector{Symbol}
  vartypes::Vector{DataType}
  threads::Int
end

function randsetup(domain::Domain, data, threads)
  geotable, names, types = _extract(data)
  RandSetup(domain, geotable, collect(names), collect(types), threads)
end

function _extract(geotable::AbstractGeoTable)
  table = values(geotable)
  sch = Tables.schema(table)
  geotable, sch.names, sch.types
end

function _extract(pairs)
  if !(eltype(pairs) <: Pair{Symbol,DataType})
    throw(ArgumentError("the data argument must be a geotable or an iterable of pairs"))
  end
  nothing, first.(pairs), last.(pairs)
end

#-----------------
# IMPLEMENTATIONS
#-----------------

include("field/sequential.jl")
include("field/gaussian.jl")
include("field/lindgren.jl")
include("field/quilting.jl")
include("field/turing.jl")
include("field/stratigraphy.jl")

#-----------------
# POINT PROCESSES
#-----------------

"""
    ishomogeneous(process::PointProcess)

Tells whether or not the spatial point process `process` is homogeneous.
"""
ishomogeneous(process::PointProcess) = false

"""
    rand([rng], process::PointProcess, g)
    rand([rng], process::PointProcess, g)

Generate `n` realizations of spatial point process `process`
inside geometry or domain `g`. Optionally specify the
random number generator `rng`.
"""
Base.rand(rng::AbstractRNG, p::PointProcess, g, n::Int) = [randsingle(rng, p, g) for _ in 1:n]

Base.rand(rng::AbstractRNG, p::PointProcess, g) = randsingle(rng, p, g)

Base.rand(p::PointProcess, g, n::Int) = rand(Random.default_rng(), p, g, n)

Base.rand(p::PointProcess, g) = rand(Random.default_rng(), p, g)

#-----------------
# IMPLEMENTATIONS
#-----------------

include("point/binomial.jl")
include("point/poisson.jl")
include("point/inhibition.jl")
include("point/cluster.jl")
include("point/union.jl")
