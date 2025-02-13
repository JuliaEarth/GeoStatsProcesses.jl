# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    rand([rng], process::FieldProcess, domain, [n];
         data=nothing, method=nothing,
         workers=workers(), async=false,
         showprogress=true)

Simulate `n` random realizations of the field `process` over the geospatial `domain`,
using the random number generator `rng`. If a geotable is provided as `data`, perform
conditional simulation, ensuring that all realizations match the `data` values.

The number of realizations `n` can be omitted, in which case the result is a single
geotable. If `n` is provided, the result becomes an ensemble with multiple realizations.

Some processes like the [`GaussianProcess`](@ref) can be simulated with alternative
simulation `method`s from the literature. Other processes can only be simulated with
a default method.

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
  kwargs...
)
  # perform processing step
  smethod = isnothing(method) ? defaultmethod(process, domain, data) : method
  preproc = preprocess(rng, process, smethod, domain, data)

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
  workers=workers(),
  async=false,
  showprogress=true,
  kwargs...
)
  # perform preprocessing step
  smethod = isnothing(method) ? defaultmethod(process, domain, data) : method
  preproc = preprocess(rng, process, smethod, domain, data)

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

# ----------------
# IMPLEMENTATIONS
# ----------------

include("field/gaussian.jl")
include("field/lindgren.jl")
