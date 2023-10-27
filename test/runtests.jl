using GeoStatsProcesses
using GeoStatsImages
using Variography
using GeoTables
using Meshes
using LinearAlgebra
using Random
using Test

import ImageQuilting
import TuringPatterns
using StratiGraphics: SmoothingProcess, Environment, ExponentialDuration

@testset "GeoStatsProcesses.jl" begin
  include("processes.jl")
  include("operations.jl")
end
