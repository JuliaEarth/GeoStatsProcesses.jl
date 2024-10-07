# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

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

* `workers` - Worker processes (default to `[myid()]`)
* `threads` - Number of threads (default to `cpucores()`)
* `verbose` - Show progress and other information (default to `true`)
* `async`   - Evaluate each realization on available `workers` only when requested (default to `false`) 

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
  table = randsingle(rng, process, rmethod, setup, prep)
  georef(table, domain)
end

function Base.rand(
  rng::AbstractRNG,
  process::FieldProcess,
  domain::Domain,
  data,
  nreals::Int,
  method=nothing;
  workers=[myid()],
  threads=cpucores(),
  verbose=true,
  async=false
)
  # perform preprocessing step
  setup = randsetup(domain, data, threads)
  rmethod = isnothing(method) ? defaultmethod(process, setup) : method
  prep = randprep(rng, process, rmethod, setup)

  # generate a single realization
  realization() = randsingle(rng, process, rmethod, setup, prep)

  # pool of worker processes
  pool = CachingPool(workers)

  # number of workers
  nworkers = length(pool)

  # simulation loop
  reals = if async
    map(1:nreals) do i
      w = workers[mod1(i, nworkers)]
      @spawnat w realization()
    end
  else
    if verbose
      pname = prettyname(process)
      mname = prettyname(rmethod)
      vname = join(setup.varnames, " ,", " and ")
      message = "$pname with $mname → $vname"
      @showprogress desc = message pmap(pool, 1:nreals) do _
        realization()
      end
    else
      pmap(pool, 1:nreals) do _
        realization()
      end
    end
  end

  Ensemble(domain, setup.varnames, reals, fetch=async ? fetch : identity)
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

#-----------------
# IMPLEMENTATIONS
#-----------------

include("field/gaussian.jl")
include("field/lindgren.jl")
include("field/quilting.jl")
include("field/turing.jl")
include("field/strata.jl")
