# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    GeoStatsProcess

A geostatistical process that can be simulated with the `rand` function.
"""
abstract type GeoStatsProcess end

"""
    RandMethod

A `rand` method for geostatistical processes.
"""
abstract type RandMethod end

"""
    DefaultRandMethod()

Default `rand` method used by some geostatistical processes.
"""
struct DefaultRandMethod <: RandMethod end

"""
    defaultmethod(process, setup) -> RandMethod

Returns the default method for the `process` and `setup`.
"""
defaultmethod(process, setup) = DefaultRandMethod()

#-----------------
# POINT PROCESSES
#-----------------

"""
    PointProcess <: GeoStatsProcess

Parent type of all point processes.
"""
abstract type PointProcess <: GeoStatsProcess end

"""
    ishomogeneous(process::PointProcess)

Tells whether or not the spatial point process `process` is homogeneous.
"""
ishomogeneous(process::PointProcess) = false

"""
    rand([rng], process::PointProcess, geometry, [nreals])
    rand([rng], process::PointProcess, domain, [nreals])

Generate one or `nreals` realizations of the point `process` inside
`geometry` or `domain`. Optionally specify the random number generator
`rng`.
"""
Base.rand(process::PointProcess, geomdom) = rand(Random.default_rng(), process, geomdom)

Base.rand(process::PointProcess, geomdom, nreals::Int) = rand(Random.default_rng(), process, geomdom, nreals)

Base.rand(rng::AbstractRNG, process::PointProcess, geomdom) = randsingle(rng, process, geomdom)

Base.rand(rng::AbstractRNG, process::PointProcess, geomdom, nreals::Int) = [randsingle(rng, process, geomdom) for _ in 1:nreals]

#-----------------
# IMPLEMENTATIONS
#-----------------

include("point/binomial.jl")
include("point/poisson.jl")
include("point/inhibition.jl")
include("point/cluster.jl")
include("point/union.jl")

#-----------------
# FIELD PROCESSES
#-----------------

"""
    FieldProcess <: GeoStatsProcess

Parent type of all field processes.
"""
abstract type FieldProcess <: GeoStatsProcess end

"""
    rand([rng], process::FieldProcess, domain, data, [nreals], [method]; [parameters])

Generate one or `nreals` realizations of the field `process` with `method`
over the `domain` with `data` and optional `paramaters`. Optionally, specify
the random number generator `rng` and global `parameters`.

The `data` can be a geotable, a pair, or an iterable of pairs of the form `var => T`,
where `var` is a symbol or string with the variable name and `T` is the corresponding
data type.

## Parameters

* `pool`    - Pool of worker processes (default to `[myid()]`)
* `threads` - Number of threads (default to `cpucores()`)
* `verbose` - Show progress and other information (default to `true`)

# Examples

```julia
julia> rand(process, domain, geotable, 3)
julia> rand(process, domain, :z => Float64)
julia> rand(process, domain, "z" => Float64)
julia> rand(process, domain, [:a => Float64, :b => Float64])
julia> rand(process, domain, ["a" => Float64, "b" => Float64])
```
"""
Base.rand(process::FieldProcess, domain::Domain, data, method=nothing; kwargs...) =
  rand(Random.default_rng(), process, domain, data, method; kwargs...)

Base.rand(process::FieldProcess, domain::Domain, data, nreals::Int, method=nothing; kwargs...) =
  rand(Random.default_rng(), process, domain, data, nreals, method; kwargs...)

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

function Base.rand(
  rng::AbstractRNG,
  process::FieldProcess,
  domain::Domain,
  data,
  nreals::Int,
  method=nothing;
  pool=[myid()],
  threads=cpucores(),
  verbose=true
)
  # perform preprocessing step
  setup = randsetup(domain, data, threads)
  rmethod = isnothing(method) ? defaultmethod(process, setup) : method
  prep = randprep(rng, process, rmethod, setup)

  # pool of worker processes
  pool = CachingPool(pool)

  # simulation loop
  reals = if verbose
    pname = prettyname(process)
    mname = prettyname(rmethod)
    vname = join(setup.varnames, " ,", " and ")
    message = "$pname with $mname → $vname"
    @showprogress desc = message pmap(pool, 1:nreals) do _
      randsingle(rng, process, rmethod, setup, prep)
    end
  else
    pmap(pool, 1:nreals) do _
      randsingle(rng, process, rmethod, setup, prep)
    end
  end

  varreals = (; (var => [real[var] for real in reals] for var in setup.varnames)...)
  Ensemble(domain, varreals)
end

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
  schema = Tables.schema(values(geotable))
  geotable, schema.names, schema.types
end

_extract(pair::Pair{Symbol,DataType}) = nothing, [first(pair)], [last(pair)]
_extract(pair::Pair{T,DataType}) where {T<:AbstractString} = nothing, [Symbol(first(pair))], [last(pair)]

_extract(pairs) = _extract(eltype(pairs), pairs)
_extract(::Type{Pair{Symbol,DataType}}, pairs) = nothing, first.(pairs), last.(pairs)
_extract(::Type{Pair{T,DataType}}, pairs) where {T<:AbstractString} = nothing, Symbol.(first.(pairs)), last.(pairs)
_extract(::Type, pairs) = throw(ArgumentError("the data argument must be a geotable, a pair, or an iterable of pairs"))

# -----------
# IO METHODS
# -----------

Base.summary(io::IO, process::GeoStatsProcess) = print(io, prettyname(process))

function Base.show(io::IO, process::GeoStatsProcess)
  name = prettyname(process)
  ioctx = IOContext(io, :compact => true)
  print(io, "$name(")
  printfields(ioctx, process, singleline=true)
  print(io, ")")
end

function Base.show(io::IO, ::MIME"text/plain", process::GeoStatsProcess)
  summary(io, process)
  printfields(io, process)
end

#-----------------
# IMPLEMENTATIONS
#-----------------

include("field/gaussian.jl")
include("field/lindgren.jl")
include("field/quilting.jl")
include("field/turing.jl")
include("field/strata.jl")
