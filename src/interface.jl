"""
    GeoStatsProcess

TODO
"""
abstract type GeoStatsProcess end

"""
    randsingle(rng::AbstractRNG, process::GeoStatsProcess, domain::Domain, varnames, vartypes, prep; threads=cpucores())

TODO
"""
function randsingle end

"""
    randprep(rng::AbstractRNG, process::GeoStatsProcess, domain::Domain, varnames, vartypes; threads=cpucores())

TODO
"""
function randprep end

function Base.rand(rng::AbstractRNG, process::GeoStatsProcess, domain::Domain, data; threads=cpucores())
  varnames, vartypes = _names_and_types(data)
  prep = randprep(rng, process, domain, varnames, vartypes; threads)
  real = randsingle(rng, process, domain, varnames, vartypes, prep; threads)
  table = (; (real[var] for var in varnames)...)
  georef(table, domain)
end

function Base.rand(
  rng::AbstractRNG,
  process::GeoStatsProcess,
  domain::Domain,
  data,
  n::Integer;
  pool=[myid()],
  progress=true,
  threads=cpucores()
)
  varnames, vartypes = _names_and_types(data)
  prep = randprep(rng, process, domain, varnames, vartypes; threads)

  # pool of worker processes
  pool = CachingPool(pool)

  # simulation loop
  reals = if progress
    message = "Simulating $(join(covars.names, " ,", " and ")):"
    @showprogress desc = message pmap(pool, 1:n) do _
      randsingle(rng, process, domain, varnames, vartypes, prep; threads)
    end
  else
    pmap(pool, 1:n) do _
      randsingle(rng, process, domain, varnames, vartypes, prep; threads)
    end
  end

  # rearrange realizations
  varvects = [Vector{V}[] for V in vartypes]
  varreals = (; zip(varnames, varvects)...)
  for real in reals
    for var in varnames
      push!(varreals[var], real[var])
    end
  end

  Ensemble(domain, varreals)
end

function _names_and_types(gtb::AbstractGeoTable)
  table = values(gtb)
  sch = Tables.schema(table)
  sch.names, sch.types
end

_names_and_types(pairs) = first.(pairs), last.(pairs)
