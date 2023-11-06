# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

module GeoStatsProcessesImageQuiltingExt

using Tables
using Random
using Meshes
using GeoTables

using GeoStatsBase: initbuff
using ImageQuilting: iqsim

using GeoStatsProcesses: RandSetup
using GeoStatsProcesses: QuiltingProcess, DefaultRandMethod

import GeoStatsProcesses: randprep, randsingle

function randprep(::AbstractRNG, process::QuiltingProcess, ::DefaultRandMethod, setup::RandSetup)
  # retrieve domain info
  sdomain = setup.domain
  simsize = size(sdomain)
  lin2cart = CartesianIndices(simsize)
  Dim = embeddim(sdomain)

  # retrieve process paramaters
  TI = process.trainimg
  init = process.init
  # default overlap
  overlap = isnothing(process.overlap) ? ntuple(i -> 1 / 6, Dim) : process.overlap

  # initialize buffers for realization and simulation mask
  vars = Dict(zip(setup.varnames, setup.vartypes))
  buff, mask = initbuff(sdomain, vars, init, data=setup.geotable)

  pairs = map(setup.varnames) do var
    # training image as simple array
    trainimg = asarray(TI, var)

    # create soft data object
    soft = if !isnothing(process.soft)
      data, dataTI = process.soft
      @assert domain(data) == sdomain "incompatible soft data for target domain"
      @assert domain(dataTI) == domain(TI) "incompatible soft data for training image"
      schema = Tables.schema(values(data))
      schemaTI = Tables.schema(values(dataTI))
      vars = schema.names |> collect |> sort
      varsTI = schemaTI.names |> collect |> sort
      @assert vars == varsTI "variables for target domain and training image differ"
      [(asarray(data, var), asarray(dataTI, var)) for var in vars]
    else
      []
    end

    # create hard data object
    linds = findall(mask[var])
    cinds = [lin2cart[ind] for ind in linds]
    hvals = view(buff[var], linds)
    hard = Dict{CartesianIndex{Dim},Real}()
    for (ind, val) in zip(cinds, hvals)
      push!(hard, ind => val)
    end

    # disable inactive voxels
    if !isnothing(process.inactive)
      for icoords in process.inactive
        push!(hard, icoords => NaN)
      end
    end

    var => (trainimg, simsize, overlap, soft, hard)
  end

  Dict(pairs)
end

function randsingle(rng::AbstractRNG, process::QuiltingProcess, ::DefaultRandMethod, setup::RandSetup, prep)
  pairs = map(setup.varnames) do var
    # unpack parameters
    threads = setup.threads
    (; tilesize, path, tol) = process
    trainimg, simsize, overlap, soft, hard = prep[var]

    # run image quilting core function
    reals = iqsim(
      trainimg,
      tilesize,
      simsize;
      overlap,
      path,
      soft,
      hard,
      tol,
      threads,
      rng
    )

    # replace NaN with missing
    vals = [isnan(v) ? missing : v for v in reals[1]]
    # flatten result
    var => vec(vals)
  end

  Dict(pairs)
end

end
