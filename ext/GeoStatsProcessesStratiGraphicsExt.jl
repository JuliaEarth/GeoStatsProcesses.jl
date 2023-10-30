# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

module GeoStatsProcessesStratiGraphicsExt

using Random
using Meshes

using StratiGraphics: LandState, Strata, simulate, voxelize

using GeoStatsProcesses: SP, RandSetup, DefaultRandMethod

import GeoStatsProcesses: randprep, randsingle

function randprep(::AbstractRNG, process::SP, ::DefaultRandMethod, setup::RandSetup)
  # retrieve domain info
  domain = setup.domain

  @assert embeddim(domain) == 3 "process implemented for 3D domain only"

  pairs = map(setup.varnames) do var
    # determine initial state
    state = if isnothing(process.state)
      nx, ny, _ = size(domain)
      land = zeros(nx, ny)
      LandState(land)
    else
      process.state
    end

    # save preprocessed input
    var => state
  end

  Dict(pairs)
end

function randsingle(::AbstractRNG, process::SP, ::DefaultRandMethod, setup::RandSetup, prep)
  # retrieve domain info
  domain = setup.domain
  _, __, nz = size(domain)

  # retrieve process parameters
  (; environment, stack, nepochs, fillbase, filltop) = process

  pairs = map(setup.varnames) do var
    # get parameters for the variable
    state = prep[var]

    # simulate the environment
    record = simulate(environment, state, nepochs)

    # build stratigraphy
    strata = Strata(record, stack)

    # return voxel model
    model = voxelize(strata, nz; fillbase, filltop)

    # flatten result
    var => vec(model)
  end

  Dict(pairs)
end


end
