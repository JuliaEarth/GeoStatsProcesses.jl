# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

module GeoStatsProcessesStratiGraphicsExt

using Random
using Meshes

using StratiGraphics: LandState, Strata, simulate, voxelize

using GeoStatsProcesses: StrataProcess

import GeoStatsProcesses: preprocess, randsingle

function preprocess(::AbstractRNG, process::StrataProcess, ::Nothing, init, domain, data)
  # sanity checks
  @assert paramdim(domain) == 3 "Stratigraphy process only implemented for 3D domains"
  @assert isnothing(data) "Stratigraphy process does not support conditional simulation"

  # determine initial state
  state = if isnothing(process.state)
    nx, ny, _ = size(domain)
    land = zeros(nx, ny)
    LandState(land)
  else
    process.state
  end

  state
end

function randsingle(::AbstractRNG, process::StrataProcess, ::Nothing, domain, data, preproc)
  # unpack preprocessed results
  state = preproc

  # process parameters
  (; environment, stack, nepochs) = process

  # retrieve domain size
  _, __, nz = size(domain)

  # simulate the environment
  record = simulate(environment, state, nepochs)

  # build stratigraphy
  strata = Strata(record, stack)

  # return voxel model
  model = voxelize(strata, nz)

  # replace NaN with missing
  vals = [isnan(v) ? missing : v for v in model]

  # flatten result
  (; Z=vec(vals))
end

end
