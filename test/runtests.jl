using GeoStatsProcesses
using Meshes
using Unitful
using GeoTables
using GeoStatsFunctions
using GeoStatsImages
using LinearAlgebra
using Test, StableRNGs

import ImageQuilting
import TuringPatterns
using StratiGraphics: SmoothingProcess, Environment, ExponentialDuration

# list of tests
testfiles = ["point.jl", "field.jl"]

@testset "GeoStatsProcesses.jl" begin
  for testfile in testfiles
    include(testfile)
  end
end
