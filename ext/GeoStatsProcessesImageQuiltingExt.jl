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
using GeoStatsProcesses: randinit

import GeoStatsProcesses: preprocess, randsingle

# helper function
getarray(gtb, var) = reshape(getproperty(gtb, var), size(domain(gtb)))

function preprocess(::AbstractRNG, process::QuiltingProcess, ::Nothing, init, dom, data)
  # parent domain
  grid = parent(dom)

  # active and inactive indices
  ainds = parentindices(dom)
  iinds = setdiff(1:nelements(dom), ainds)

  # sanity checks
  @assert grid isa Grid "quilting process only defined for grids or views of grids"

  # initialize realization and mask
  real, mask = randinit(process, grid, data, init)

  # retrieve variable names
  vars = keys(real)

  # multivariate simulation is not supported
  @assert length(vars) == 1 "quilting process does not support multivariate simulation"

  # retrieve variable name
  var = first(vars)

  # training image
  TI = process.trainimg

  # training image as simple array
  trainimg = getarray(TI, var)

  # simulation size
  simsize = size(grid)

  # default overlap
  overlap = isnothing(process.overlap) ? ntuple(i -> 1 / 6, paramdim(grid)) : process.overlap

  # create soft data object
  soft = if !isnothing(process.soft)
    dat, datTI = process.soft
    @assert domain(dat) == dom "incompatible soft data for target domain"
    @assert domain(datTI) == domain(TI) "incompatible soft data for training image"
    vars₁ = dat |> values |> Tables.columns |> Tables.columnnames
    vars₂ = datTI |> values |> Tables.columns |> Tables.columnnames
    @assert Set(vars₁) == Set(vars₂) "auxiliary variables for target domain and training image differ"
    [(getarray(dat, var), getarray(datTI, var)) for var in vars₁]
  else
    []
  end

  # convert from linear to Cartesian indices
  cart = CartesianIndices(simsize)

  # hard data in dictionary form
  hinds = findall(mask[var])
  hvals = view(real[var], hinds)
  hdata = Dict(cart[hinds] .=> hvals)

  # inactive voxels in dictionary form
  idata = Dict(cart[iinds] .=> NaN)

  # create hard data object
  hard = merge(hdata, idata)

  (; var, ainds, trainimg, simsize, overlap, soft, hard)
end

function randsingle(rng::AbstractRNG, process::QuiltingProcess, ::Nothing, domain, data, preproc)
  # unpack preprocessing results
  (; var, ainds, trainimg, simsize, overlap, soft, hard) = preproc

  # process parameters
  (; tilesize, path, tol) = process

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
    rng
  ) |> first

  # replace NaN with missing
  vals = [isnan(v) ? missing : v for v in real[ainds]]

  # flatten result
  (; var => vec(vals))
end

end
