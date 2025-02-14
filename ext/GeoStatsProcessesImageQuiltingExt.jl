# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

module GeoStatsProcessesImageQuiltingExt

using Tables
using Random
using Meshes
using GeoTables

using ImageQuilting: iqsim

using GeoStatsProcesses: QuiltingProcess
using GeoStatsProcesses: DefaultSimulation
using GeoStatsProcesses: randinit

import GeoStatsProcesses: preprocess, randsingle

function preprocess(::AbstractRNG, process::QuiltingProcess, ::DefaultSimulation, init, dom, data)
  # sanity checks
  @assert dom isa Grid "image quilting only defined for grids"

  # initialize realization and mask
  real, mask = randinit(process, dom, data, init)

  # multivariate simulation is not supported
  @assert length(keys(real)) == 1 "image quilting does not support multivariate simulation"

  # retrieve variable name
  var = first(keys(real))

  # training image
  TI = process.trainimg

  # training image as simple array
  trainimg = asarray(TI, var)

  # simulation size
  simsize = size(dom)

  # default overlap
  overlap = isnothing(process.overlap) ? ntuple(i -> 1 / 6, paramdim(dom)) : process.overlap

  # create soft data object
  soft = if !isnothing(process.soft)
    data, dataTI = process.soft
    @assert domain(data) == dom "incompatible soft data for target domain"
    @assert domain(dataTI) == domain(TI) "incompatible soft data for training image"
    vars = data |> values |> Tables.columns |> Tables.columnnames
    varsTI = dataTI |> values |> Tables.columns |> Tables.columnnames
    @assert Set(vars) == Set(varsTI) "auxiliary variables for target domain and training image differ"
    [(asarray(data, var), asarray(dataTI, var)) for var in vars]
  else
    []
  end

  # create hard data object
  linds = findall(mask[var])
  icart = CartesianIndices(simsize)
  cinds = [icart[ind] for ind in linds]
  hvals = view(real[var], linds)
  hard = Dict(cinds .=> hvals)

  # disable inactive voxels
  if !isnothing(process.inactive)
    inactive = Dict(process.inactive .=> NaN)
    hard = merge(hard, inactive)
  end

  (; var, trainimg, simsize, overlap, soft, hard)
end

function randsingle(rng::AbstractRNG, process::QuiltingProcess, ::DefaultSimulation, domain, data, preproc)
  # unpack preprocessing results
  (; var, trainimg, simsize, overlap, soft, hard) = preproc

  # process parameters
  (; tilesize, path, tol, nthreads) = process

  # temporary adjustments
  threads = nthreads

  # run image quilting main function
  real = iqsim(
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
  ) |> first

  # replace NaN with missing
  vals = [isnan(v) ? missing : v for v in real]

  # flatten result
  (; var => vec(vals))
end

end
