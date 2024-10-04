using GeoStatsProcesses
using Meshes
using Unitful
using GeoTables
using GeoStatsFunctions
using GeoStatsImages
using LinearAlgebra
using Test, StableRNGs, CSV

import ImageQuilting
import TuringPatterns
using StratiGraphics: SmoothingProcess, Environment, ExponentialDuration

# environment settings
datadir = joinpath(@__DIR__, "data")

# list of tests
testfiles = ["point.jl", "field.jl", "ensembles.jl", "initbuff.jl"]

@testset "GeoStatsProcesses.jl" begin
  for testfile in testfiles
    include(testfile)
  end
end
