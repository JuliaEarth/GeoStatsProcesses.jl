# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

module GeoStatsProcessesTuringPatternsExt

using Random
using Meshes

using TuringPatterns: BoxBlur, Clamp, SimplePattern, Sim
using TuringPatterns: Params, step!, scale01

using GeoStatsProcesses: TuringProcess

import GeoStatsProcesses: preprocess, randsingle

const PARAMS1 =
  [Params(2, 4, 0.01), Params(5, 10, 0.02), Params(10, 20, 0.03), Params(20, 40, 0.04), Params(50, 100, 0.05)]

function preprocess(::AbstractRNG, process::TuringProcess, ::Nothing, init, domain, data)
  # sanity checks
  @assert domain isa Grid "Turing process only defined for grids"
  @assert isnothing(data) "Turing process does not support conditional simulation"

  # grid size
  dims = size(domain)

  # simulation parameters
  params = isnothing(process.params) ? PARAMS1 : process.params

  # construct patterns from parameters
  patterns = [SimplePattern(param, dims) for param in params]

  patterns
end

function randsingle(::AbstractRNG, process::TuringProcess, ::Nothing, domain, data, preproc)
  # unpack preprocessed results
  patterns = preproc

  # grid size
  dims = size(domain)

  # process parameters
  blur = isnothing(process.blur) ? BoxBlur(dims) : process.blur(dims)
  edge = isnothing(process.edge) ? Clamp() : process.edge()
  iter = process.iter

  # perform simulation
  sim = Sim(rand(Float64, dims), patterns, edge, blur)
  for _ in 1:iter
    step!(sim)
  end
  vals = scale01(sim.fluid)

  # flatten result
  (; field=vec(vals))
end

end
