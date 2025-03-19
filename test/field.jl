@testset "FieldProcess" begin
  @testset "GaussianProcess" begin
    @testset "defaultsimulation" begin
      proc1 = GaussianProcess(GaussianVariogram())
      proc2 = GaussianProcess(GaussianCovariance())
      grid = CartesianGrid(100, 100)
      vgrid = view(grid, 1:1000)
      pset1 = PointSet(rand(Point, 1000))
      pset2 = PointSet(rand(Point, 10000))
      @test GeoStatsProcesses.defaultsimulation(proc1, grid) isa FFTSIM
      @test GeoStatsProcesses.defaultsimulation(proc1, vgrid) isa FFTSIM
      @test GeoStatsProcesses.defaultsimulation(proc1, pset1) isa SEQSIM
      @test GeoStatsProcesses.defaultsimulation(proc1, pset2) isa SEQSIM
      @test GeoStatsProcesses.defaultsimulation(proc2, pset1) isa LUSIM
    end

    @testset "LUSIM" begin
      method = LUSIM()

      # simulation over 1D domain
      rng = StableRNG(123)
      proc = GaussianProcess(SphericalCovariance(range=10.0))
      grid = CartesianGrid(100)
      data = georef((; Z=[0.0, 1.0, 0.0, 1.0, 0.0]), [(0.0,), (25.0,), (50.0,), (75.0,), (100.0,)])
      real = rand(rng, proc, grid; method)
      @test eltype(real.field) <: Float64
      real = rand(rng, proc, grid; method, data)
      @test eltype(real.Z) <: Float64

      # cosimulation
      rng = StableRNG(123)
      func = [1.0 0.95; 0.95 1.0] * GaussianCovariance(range=10.0)
      mean = [0.0, 0.0]
      proc = GaussianProcess(func, mean)
      real = rand(rng, proc, grid; method)
      @test eltype(real.field1) <: Float64
      @test eltype(real.field2) <: Float64

      # simulation over 2D domain
      rng = StableRNG(123)
      proc = GaussianProcess(GaussianCovariance(range=10.0))
      grid = CartesianGrid(100, 100)
      real = rand(rng, proc, grid; method)
      @test eltype(real.field) <: Float64
    end

    @testset "SEQSIM" begin
      method = SEQSIM(maxneighbors=3)

      # basic simulation
      rng = StableRNG(2017)
      proc = GaussianProcess(SphericalVariogram(range=35.0))
      grid = CartesianGrid((100, 100), (0.5, 0.5), (1.0, 1.0))
      data = georef((; Z=[1.0, 0.0, 1.0]), [(25.0, 25.0), (50.0, 75.0), (75.0, 50.0)])
      real = rand(rng, proc, grid; method)
      @test eltype(real.field) <: Float64
      real = rand(rng, proc, grid; method, data)
      @test eltype(real.Z) <: Float64
      inds = LinearIndices(size(grid))
      @test real.Z[inds[25, 25]] == 1.0

      # simulation with units
      rng = StableRNG(2017)
      mean = 0.0u"K"
      func = SphericalVariogram(range=35.0, sill=1.0u"K^2")
      proc = GaussianProcess(func, mean)
      grid = CartesianGrid(10, 10)
      data = georef((; Z=[1.0, 0.0, 1.0] * u"K"), [(25.0, 25.0), (50.0, 75.0), (75.0, 50.0)])
      real = rand(rng, proc, grid; method)
      @test unit(eltype(real.field)) == u"K"
      real = rand(rng, proc, grid; method, data)
      @test unit(eltype(real.Z)) == u"K"

      # multivariate simulation with units
      rng = StableRNG(2017)
      mean = [0.1u"ppm", 0.2u"ppm"]
      func = [1.0 0.8; 0.8 1.0] * SphericalCovariance(range=35.0, sill=1.0u"ppm^2")
      proc = GaussianProcess(func, mean)
      grid = CartesianGrid(10, 10)
      data = georef((; Cu=[0.0, 0.1, 0.0] * u"ppm", Zn=[0.1, 0.0, 0.1] * u"ppm"), [(25.0, 25.0), (50.0, 75.0), (75.0, 50.0)])
      real = rand(rng, proc, grid; method)
      @test unit(eltype(real.field1)) == u"ppm"
      @test unit(eltype(real.field2)) == u"ppm"
      real = rand(rng, proc, grid; method, data)
      @test unit(eltype(real.Cu)) == u"ppm"
      @test unit(eltype(real.Zn)) == u"ppm"
    end

    @testset "FFTSIM" begin
      method = FFTSIM()

      # basic simulation
      rng = StableRNG(2019)
      proc = GaussianProcess(GaussianVariogram(range=10.0))
      grid = CartesianGrid(100, 100)
      real = rand(rng, proc, grid; method)
      @test eltype(real.field) <: Float64

      # simulation on view of grid
      rng = StableRNG(2019)
      proc = GaussianProcess(GaussianVariogram(range=10.0))
      grid = CartesianGrid(100, 100)
      vgrid = view(grid, 1:5000)
      real = rand(rng, proc, vgrid; method)
      @test domain(real) == vgrid
      @test length(real.geometry) == 5000

      # conditional simulation
      rng = StableRNG(2019)
      proc = GaussianProcess(GaussianVariogram(range=35.0))
      grid = CartesianGrid((100, 100), (0.5, 0.5), (1.0, 1.0))
      data = georef((; Z=[1.0, 0.0, 1.0]), [(25.0, 25.0), (50.0, 75.0), (75.0, 50.0)])
      real = rand(rng, proc, grid, method=FFTSIM(maxneighbors=3), data=data)
      @test domain(real) == grid
      @test nrow(real) == 10000
    end
  end

  @testset "IndicatorProcess" begin
    # two categorical values
    rng = StableRNG(123)
    proc = IndicatorProcess(SphericalTransiogram(range=35.0))
    grid = CartesianGrid(100, 100)
    data = georef((; C=[1, 3, 1]), [(25.0, 25.0), (50.0, 75.0), (75.0, 50.0)])
    real = rand(rng, proc, grid)
    @test eltype(real.field) == Int
    @test Set(real.field) == Set([1, 2])
    real = rand(rng, proc, grid; data)
    @test eltype(real.C) == Int
    @test Set(real.C) == Set([1, 3])

    # three categorical values
    rng = StableRNG(123)
    proc = IndicatorProcess(SphericalTransiogram(range=35.0, proportions=(0.7, 0.2, 0.1)))
    grid = CartesianGrid(100, 100)
    data = georef((; C=["a", "b", "c"]), [(25.0, 25.0), (50.0, 75.0), (75.0, 50.0)])
    real = rand(rng, proc, grid)
    @test Set(real.field) == Set([1, 2, 3])
    real = rand(rng, proc, grid; data)
    @test eltype(real.C) == String
    @test Set(real.C) == Set(["a", "b", "c"])
  end

  @testset "LindgrenProcess" begin
    proc = LindgrenProcess()
    @test proc.range == 1u"m"
    @test proc.sill == 1.0

    proc = LindgrenProcess(2.0)
    @test proc.range == 2.0u"m"
    @test proc.sill == 1.0

    proc = LindgrenProcess(2.0u"km")
    @test proc.range == 2.0u"km"
    @test proc.sill == 1.0

    proc = LindgrenProcess(2.0, 2.0)
    @test proc.range == 2.0u"m"
    @test proc.sill == 2.0

    proc = LindgrenProcess(2.0u"km", 2.0)
    @test proc.range == 2.0u"km"
    @test proc.sill == 2.0

    proc = LindgrenProcess(2, 2)
    @test Unitful.numtype(proc.range) == Float64
    @test proc.sill isa Float64

    # simulation on sphere
    proc = LindgrenProcess(0.1)
    mesh = simplexify(Sphere((0, 0, 0)))
    data = georef((Z=[0.0, 1.0],), [(0, 0, -1), (0, 0, 1)])
    # unconditional realization
    rng = StableRNG(2024)
    real = rand(rng, proc, mesh)
    @test isapprox(sum(real.field) / length(real.field), 0.0, atol=1e-3)
    # conditional realization
    rng = StableRNG(2024)
    real = rand(rng, proc, mesh; data)
    @test isapprox(sum(real.Z) / length(real.Z), 0.0, atol=1e-3)
  end

  @testset "QuiltingProcess" begin
    rng = StableRNG(2017)
    data = georef((; facies=[1.0, 0.0, 1.0]), [(25.0, 25.0), (50.0, 75.0), (75.0, 50.0)])
    grid = CartesianGrid(100, 100)
    trainimg = geostatsimage("Strebelle")
    proc = QuiltingProcess(trainimg, (30, 30))

    # simulation on full grid
    real = rand(rng, proc, grid; data)
    @test size(domain(real)) == (100, 100)
    @test eltype(real.facies) == Float64

    # simulation on grid view
    vgrid = view(grid, Ball((50, 50), 20))
    real = rand(rng, proc, vgrid)
    @test domain(real) == vgrid
    @test length(real.facies) == nelements(vgrid)
  end

  @testset "TuringProcess" begin
    rng = StableRNG(2019)
    grid = CartesianGrid(200, 200)
    real = rand(rng, TuringProcess(), grid)
    @test size(domain(real)) == (200, 200)
    @test eltype(real.field) <: Float64
  end

  @testset "StrataProcess" begin
    rng = StableRNG(2019)
    proc = SmoothingProcess()
    env = Environment(rng, [proc, proc], [0.5 0.5; 0.5 0.5], ExponentialDuration(rng, 1.0))
    grid = CartesianGrid(50, 50, 20)
    real = rand(StrataProcess(env), grid)
    @test size(domain(real)) == (50, 50, 20)
    @test eltype(real.field) <: Union{Float64,Missing}
    @test any(ismissing, real.field)
    @test all(!isnan, skipmissing(real.field))
  end

  @testset "parallel simulation" begin
    addprocs(2)

    @everywhere using GeoStatsProcesses, StableRNGs

    rng = StableRNG(2019)
    grid = CartesianGrid(100, 100)
    proc = GaussianProcess(GaussianVariogram())

    # synchronous
    real = rand(rng, proc, grid, 3, showprogress=false)
    @test length(real) == 3
    @test domain(real[1]) == grid
    @test eltype(real[1].field) <: Float64

    # asynchronous
    real = rand(rng, proc, grid, 3, async=true)
    @test length(real) == 3
    @test domain(real[1]) == grid
    @test eltype(real[1].field) <: Float64

    # async option is not allowed when the master is in the workers
    @test_throws ArgumentError rand(rng, proc, grid, 3, workers=[myid()], async=true)

    rmprocs(workers()...)
  end
end
