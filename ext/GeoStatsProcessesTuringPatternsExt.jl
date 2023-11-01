# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

module GeoStatsProcessesTuringPatternsExt

using Random
using Meshes

using TuringPatterns: BoxBlur, Clamp, SimplePattern, Sim
using TuringPatterns: Params, step!, scale01

using GeoStatsProcesses: RandSetup
using GeoStatsProcesses: TuringProcess, DefaultRandMethod

import GeoStatsProcesses: randprep, randsingle

const PARAMS1 =
  [Params(2, 4, 0.01), Params(5, 10, 0.02), Params(10, 20, 0.03), Params(20, 40, 0.04), Params(50, 100, 0.05)]

function randprep(::AbstractRNG, process::TuringProcess, ::DefaultRandMethod, setup::RandSetup)
  # retrieve domain of simulation
  domain = setup.domain
  topo = topology(domain)

  # assert grid topology
  @assert topo isa GridTopology "process only defined over grid topology"

  # retrieve simulation size
  sz = size(topo)

  # retrieve simulation parameters
  params = isnothing(process.params) ? PARAMS1 : process.params

  pairs = map(setup.varnames) do var
    # construct patterns from parameters
    patterns = [SimplePattern(param, sz) for param in params]
    var => patterns
  end

  Dict(pairs)
end

function randsingle(::AbstractRNG, process::TuringProcess, ::DefaultRandMethod, setup::RandSetup, prep)
  # retrieve domain size
  sz = size(topology(setup.domain))

  # retrieve process parameters
  blur = isnothing(process.blur) ? BoxBlur(sz) : process.blur(sz)
  edge = isnothing(process.edge) ? Clamp() : process.edge()
  iter = process.iter

  pairs = map(setup.vartypes, setup.varnames) do V, var
    # unpack preprocessed parameters
    patterns = prep[var]

    # perform simulation
    sim = Sim(rand(V, sz), patterns, edge, blur)
    for _ in 1:iter
      step!(sim)
    end
    real = scale01(sim.fluid)

    # flatten result
    var => vec(real)
  end

  Dict(pairs)
end

end
