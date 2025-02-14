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
  smethod = isnothing(method) ? defaultsimulation(process, domain) : method
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
  smethod = isnothing(method) ? defaultsimulation(process, domain) : method
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
      message = "$pname with $mname"
      @showprogress desc = message pmap(pool, 1:nreals) do _
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

"""
    randinit(process, domain, data, init)

Initialize attribute table of realization based on the field `process`,
the geospatial `domain` and the geospatial `data` using an `init`ialization
method.
"""
function randinit(process::FieldProcess, domain, data, init)
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

function dataschema(data)
  schema = data |> values |> Tables.columns |> Tables.schema
  Tables.Schema(schema.names, map(nonmissingtype, schema.types))
end

# ----------------
# IMPLEMENTATIONS
# ----------------

include("field/gaussian.jl")
include("field/lindgren.jl")

# ---------
# DEFAULTS
# ---------

"""
    defaultschema(process)

Default schema of realizations of field `process`.
"""
defaultschema(::FieldProcess) = Tables.Schema((:Z,), (Float64,))

function defaultschema(process::GaussianProcess)
  nvars = nvariates(process.func)
  names = nvars > 1 ? ntuple(i -> Symbol(:Z, i), nvars) : (:Z,)
  types = ntuple(i -> Float64, nvars)
  Tables.Schema(names, types)
end

function defaultschema(process::QuiltingProcess)
  table = process.trainimg |> values
  table |> Tables.columns |> Tables.schema
end

"""
    defaultsimulation(process, domain)

Default method used for the simulation of geostatistical `process`
over given geospatial `domain`.
"""
defaultsimulation(::FieldProcess, domain) = DefaultSimulation()

function defaultsimulation(process::GaussianProcess, domain)
  d = domain
  p = parent(d)
  b = boundingbox(p)
  f = process.func
  if p isa Grid && range(f) ≤ minimum(sides(b)) / 3
    FFTSIM()
  elseif nelements(d) < 100 * 100
    LUSIM()
  else
    SEQSIM()
  end
end
