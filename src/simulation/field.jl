# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    rand([rng], process::FieldProcess, domain, [n];
         data=nothing, method=nothing, init=NearestInit(),
         workers=workers(), async=false, showprogress=true)

Simulate `n` random realizations of the field `process` over the geospatial `domain`,
using the random number generator `rng`. If a geotable is provided as `data`, perform
conditional simulation, ensuring that all realizations match the `data` values.

The number of realizations `n` can be omitted, in which case the result is a single
geotable. If `n` is provided, the result becomes an ensemble with multiple realizations.

Some processes like the [`GaussianProcess`](@ref) can be simulated with alternative
simulation `method`s from the literature such as [`LUSIM`](@ref), [`SEQSIM`](@ref)
and [`FFTSIM`](@ref). Other processes can only be simulated with a default method.
A suitable `method` is automatically chosen based on the `process`, `domain` and `data`.

In the presence of `data`, the realizations are initialized with a `init`ialization
method. By default, the data is assigned to the nearest geometry of the simulation
domain.

Multiple `workers` created with the `Distributed` standard library can be used for
parallel simulation. The function can be called `async`hrounously, in which case it
returns immediately for online consumption of realizations in a loop. The option
`showprogress` can be used to show a progress bar with estimated time of conclusion.

## Examples

```julia
rand(process, domain) # single realization as geotable
rand(process, domain, data=data) # conditional simulation
rand(process, domain, 3) # ensemble of realizations
rand(process, domain, 10, data=data) # ensemble of realizations
rand(rng, process, domain, 3) # specify random number generator
```

### Notes

If you are not an expert in random fields, avoid setting the `method` manually.
There are good heuristics in place to maximize performance in any given setup.
"""
Base.rand(process::FieldProcess, domain::Domain; kwargs...) =
  rand(Random.default_rng(), process, domain; kwargs...)

Base.rand(process::FieldProcess, domain::Domain, nreals::Int; kwargs...) =
  rand(Random.default_rng(), process, domain, nreals; kwargs...)

function Base.rand(
  rng::AbstractRNG,
  process::FieldProcess,
  domain::Domain;
  data=nothing,
  method=nothing,
  init=NearestInit(),
  kwargs...
)
  # perform processing step
  smethod = isnothing(method) ? defaultsimulation(process, domain; data) : method
  preproc = preprocess(rng, process, smethod, init, domain, data)

  # simulate a single realization
  table = randsingle(rng, process, smethod, domain, data, preproc)

  # return geotable
  georef(table, domain)
end

function Base.rand(
  rng::AbstractRNG,
  process::FieldProcess,
  domain::Domain,
  nreals::Int;
  data=nothing,
  method=nothing,
  init=NearestInit(),
  workers=workers(),
  async=false,
  showprogress=true,
)
  # perform preprocessing step
  smethod = isnothing(method) ? defaultsimulation(process, domain; data) : method
  preproc = preprocess(rng, process, smethod, init, domain, data)

  # simulate a single realization
  realization() = randsingle(rng, process, smethod, domain, data, preproc)

  # pool of worker processes
  pool = CachingPool(workers)

  # number of workers
  nworkers = length(pool)

  if async && myid() ∈ workers
    throw(ArgumentError("the `async` option is not allowed when the master process is in the `workers`"))
  end

  # simulation loop
  reals = if async
    map(1:nreals) do i
      w = workers[mod1(i, nworkers)]
      @spawnat w realization()
    end
  else
    if showprogress
      pname = prettyname(process)
      mname = prettyname(smethod)
      desc = isnothing(smethod) ? "$pname" : "$pname ($mname)"
      @showprogress desc = desc pmap(pool, 1:nreals) do _
        realization()
      end
    else
      pmap(pool, 1:nreals) do _
        realization()
      end
    end
  end

  Ensemble(domain, reals, fetch=async ? fetch : identity)
end

# --------
# METHODS
# --------

"""
    FieldSimulationMethod

A method for simulating field processes.
"""
abstract type FieldSimulationMethod end

include("field/lusim.jl")
include("field/seqsim.jl")
include("field/fftsim.jl")
include("field/defsim.jl")

# ---------
# DEFAULTS
# ---------

"""
    defaultsimulation(process, domain; data=nothing)

Default method used for the simulation of geostatistical `process`
over given geospatial `domain` and optional `data`.
"""
defaultsimulation(process::FieldProcess, domain; data=nothing) = nothing

function defaultsimulation(process::GaussianProcess, domain; data=nothing)
  d = domain
  p = parent(d)
  b = boundingbox(p)
  f = process.func
  if p isa Grid && isstationary(f) && nvariates(f) == 1 && range(f) ≤ minimum(sides(b)) / 3 && isnothing(data)
    FFTSIM()
  elseif nelements(d) < 100 * 100 && isstationary(f) && issymmetric(f) && isbanded(f)
    LUSIM()
  else
    SEQSIM()
  end
end

function defaultsimulation(process::IndicatorProcess, dom; data=nothing)
  path = if isnothing(data)
    LinearPath()
  else
    d = domain(data)
    s = KNearestSearch(dom, 1)
    inds = map(1:nelements(d)) do i
      p = centroid(d, i)
      first(search(p, s))
    end
    SourcePath(unique(inds))
  end
  SEQSIM(; path)
end
