using GeoStatsProcesses
using Variography
using GeoTables
using Meshes
using LinearAlgebra
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

  @testset "SGP" begin
    𝒮 = georef((; z=[1.0, 0.0, 1.0]), [25.0 50.0 75.0; 25.0 75.0 50.0])
    𝒟 = CartesianGrid((100, 100), (0.5, 0.5), (1.0, 1.0))
    N = 3
  
    process = SGP(variogram=SphericalVariogram(range=35.0), neighborhood=MetricBall(30.0))
  
    Random.seed!(2017)
    sims₁ = rand(process, 𝒟, 𝒮, 3)
    sims₂ = rand(process, 𝒟, [:z => Float64], 3)
  
    # basic checks
    reals = sims₁[:z]
    inds = LinearIndices(size(𝒟))
    @test all(reals[i][inds[25, 25]] == 1.0 for i in 1:N)
    @test all(reals[i][inds[50, 75]] == 0.0 for i in 1:N)
    @test all(reals[i][inds[75, 50]] == 1.0 for i in 1:N)
  end

  @testset "LUGP" begin
    𝒮 = georef((; z=[0.0, 1.0, 0.0, 1.0, 0.0]), [0.0 25.0 50.0 75.0 100.0])
    𝒟 = CartesianGrid(100)
  
    # ----------------------
    # conditional simulation
    # ----------------------
    rng = MersenneTwister(123)
    process = LUGP(variogram=SphericalVariogram(range=10.0))
    sims = rand(rng, process, 𝒟, 𝒮, 2)
  
    # ------------------------
    # unconditional simulation
    # ------------------------
    rng = MersenneTwister(123)
    process = LUGP(variogram=SphericalVariogram(range=10.0))
    sims = rand(rng, process, 𝒟, [:z => Float64], 2)
  
    # -------------
    # co-simulation
    # -------------
    𝒟 = CartesianGrid(500)  
    rng = MersenneTwister(123)
    process = LUGP(variogram=SphericalVariogram(range=10.0), correlation=0.95)
    sim = rand(rng, process, 𝒟, [:a => Float64, :b => Float64])
  
    # -----------
    # 2D example
    # -----------
    𝒟 = CartesianGrid(100, 100)  
    rng = MersenneTwister(123)
    process = LUGP(variogram=GaussianVariogram(range=10.0))
    sims = rand(rng, process, 𝒟, [:z => Float64], 3)
  
    # -------------------
    # anisotropy example
    # -------------------
    𝒟 = CartesianGrid(100, 100)
    rng = MersenneTwister(123)
    ball = MetricBall((20.0, 5.0))
    process = LUGP(variogram=GaussianVariogram(ball))
    sims = rand(rng, process, 𝒟, [:z => Float64], 3)
  
    # ---------------------
    # custom factorization
    # ---------------------
    𝒟 = CartesianGrid(100)
    rng = MersenneTwister(123)
    process1 = LUGP(variogram=SphericalVariogram(range=10.0), factorization=lu)
    process2 = LUGP(variogram=SphericalVariogram(range=10.0), factorization=cholesky)
    sim1 = rand(rng, process1, 𝒟, 𝒮, 2)
    sim2 = rand(rng, process2, 𝒟, 𝒮, 2)

    # throws
    𝒟 = CartesianGrid(100, 100)  
    process = LUGP(variogram=GaussianVariogram(range=10.0))
    # only 1 or 2 variables can be simulated simultaneously
    @test_throws AssertionError rand(process, 𝒟, [:a => Float64, :b => Float64, :c => Float64]) 
  end
end
