# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    GeoStatsProcess

TODO
"""
abstract type GeoStatsProcess end

"""
    randsingle(rng::AbstractRNG, process::GeoStatsProcess, setup::RandSetup, prep)

TODO
"""
function randsingle end

"""
    randprep(rng::AbstractRNG, process::GeoStatsProcess, setup::RandSetup)

TODO
"""
function randprep end

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

Base.rand(process::GeoStatsProcess, domain::Domain, data; kwargs...) =
  rand(Random.default_rng(), process, domain, data; kwargs...)

function Base.rand(rng::AbstractRNG, process::GeoStatsProcess, domain::Domain, data; threads=cpucores())
  setup = randsetup(domain, data, threads)
  prep = randprep(rng, process, setup)
  real = randsingle(rng, process, setup, prep)
  table = (; (var => real[var] for var in setup.varnames)...)
  georef(table, domain)
end

Base.rand(process::GeoStatsProcess, domain::Domain, data, n::Integer; kwargs...) =
  rand(Random.default_rng(), process, domain, data, n; kwargs...)

function Base.rand(
  rng::AbstractRNG,
  process::GeoStatsProcess,
  domain::Domain,
  data,
  n::Integer;
  pool=[myid()],
  threads=cpucores(),
  progress=true
)
  setup = randsetup(domain, data, threads)
  prep = randprep(rng, process, setup; threads)

  # pool of worker processes
  pool = CachingPool(pool)

  # simulation loop
  reals = if progress
    message = "Simulating $(join(covars.names, " ,", " and ")):"
    @showprogress desc = message pmap(pool, 1:n) do _
      randsingle(rng, process, setup, prep)
    end
  else
    pmap(pool, 1:n) do _
      randsingle(rng, process, setup, prep)
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
