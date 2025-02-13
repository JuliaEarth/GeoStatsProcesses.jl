using GeoStatsProcesses
using Meshes
using Unitful
using GeoTables
using GeoStatsFunctions
using GeoStatsImages
using LinearAlgebra
using Distributed
using CSV
using Test, StableRNGs

import ImageQuilting
import TuringPatterns
using StratiGraphics: SmoothingProcess, Environment, ExponentialDuration

# environment settings
datadir = joinpath(@__DIR__, "data")

# list of tests
testfiles = ["initbuff.jl", "ensembles.jl", "point.jl", "field.jl"]

@testset "GeoStatsProcesses.jl" begin
  for testfile in testfiles
    include(testfile)
  end
end
