using GeoStatsProcesses
using Variography
using GeoTables
using Meshes
using Random
using Test

@testset "GeoStatsProcesses.jl" begin
  @testset "FFTGP" begin
    # isotropic simulation
    Random.seed!(2019)
    dom = CartesianGrid(100, 100)
    process = FFTGP(variogram=GaussianVariogram(range=10.0))
    sims = rand(process, dom, [:z => Float64], 3)
  
    # anisotropic simulation
    Random.seed!(2019)
    dom = CartesianGrid(100, 100)
    process = FFTGP(variogram=GaussianVariogram(MetricBall((20.0, 5.0))))
    sims = rand(process, dom, [:z => Float64], 3)
  
    # simulation on view of grid
    Random.seed!(2022)
    grid = CartesianGrid(100, 100)
    vgrid = view(grid, 1:5000)
    process = FFTGP(variogram=GaussianVariogram(range=10.0))
    sim = rand(process, vgrid, [:z => Float64])
    @test domain(sim) == vgrid
    @test length(sim.geometry) == 5000
  
    # conditional simulation
    Random.seed!(2022)
    table = (; z=[1.0, -1.0, 1.0])
    coords = [(25.0, 25.0), (50.0, 75.0), (75.0, 50.0)]
    samples = georef(table, coords)
    sdomain = CartesianGrid(100, 100)
    process = FFTGP(variogram=GaussianVariogram(range=10.0))
    sim = rand(process, sdomain, samples)
  end
end
