@testset "FieldProcess" begin
  @testset "async" begin
    addprocs(2)

    @everywhere using GeoStatsProcesses, StableRNGs

    rng = StableRNG(2019)
    grid = CartesianGrid(100, 100)
    proc = GaussianProcess()
    real = rand(rng, proc, grid, 3, async=true)
    @test length(real) == 3
    @test domain(real[1]) == grid
    @test eltype(real[1].Z) <: Float64
    @test domain(real[2]) == grid
    @test eltype(real[2].Z) <: Float64
    @test domain(real[3]) == grid
    @test eltype(real[3].Z) <: Float64
    # error: the `async` option is not allowed when the master process is in the `workers`
    @test_throws ArgumentError rand(rng, proc, grid, 3, workers=[myid()], async=true)

    rmprocs(workers()...)
  end

  @testset "GaussianProcess" begin
    @testset "defaultsimulation" begin
      proc = GaussianProcess()
      grid = CartesianGrid(100, 100)
      vgrid = view(grid, 1:1000)
      pset1 = PointSet(rand(Point, 1000))
      pset2 = PointSet(rand(Point, 10000))
      @test GeoStatsProcesses.defaultsimulation(proc, grid, nothing) isa FFTSIM
      @test GeoStatsProcesses.defaultsimulation(proc, vgrid, nothing) isa FFTSIM
      @test GeoStatsProcesses.defaultsimulation(proc, pset1, nothing) isa LUSIM
      @test GeoStatsProcesses.defaultsimulation(proc, pset2, nothing) isa SEQSIM
    end

    @testset "LUSIM" begin
      data = georef((; Z=[0.0, 1.0, 0.0, 1.0, 0.0]), [(0.0,), (25.0,), (50.0,), (75.0,), (100.0,)])
      grid = CartesianGrid(100)

      method = LUSIM()

      # unconditional simulation
      rng = StableRNG(123)
      proc = GaussianProcess(SphericalCovariance(range=10.0))
      real = rand(rng, proc, grid; method)
      @test eltype(real.Z) <: Float64

      # conditional simulation
      rng = StableRNG(123)
      proc = GaussianProcess(SphericalCovariance(range=10.0))
      real = rand(rng, proc, grid; method, data)
      @test eltype(real.Z) <: Float64

      # cosimulation
      rng = StableRNG(123)
      grid = CartesianGrid(500)
      proc = GaussianProcess([1.0 0.95; 0.95 1.0] * GaussianCovariance(range=10.0))
      real = rand(rng, proc, grid; method)
      @test eltype(real.Z1) <: Float64
      @test eltype(real.Z2) <: Float64

      # custom factorization
      rng = StableRNG(123)
      grid = CartesianGrid(100)
      proc = GaussianProcess(SphericalCovariance(range=10.0))
      real = rand(rng, proc, grid; method=LUSIM(factorization=lu), data=data)
      @test eltype(real.Z) <: Float64

      # 2D example
      rng = StableRNG(123)
      grid = CartesianGrid(100, 100)
      proc = GaussianProcess(GaussianCovariance(range=10.0))
      real = rand(rng, proc, grid; method)
    end

    @testset "SEQSIM" begin
      proc = GaussianProcess(SphericalVariogram(range=35.0))
      grid = CartesianGrid((100, 100), (0.5, 0.5), (1.0, 1.0))
      data = georef((; Z=[1.0, 0.0, 1.0]), [(25.0, 25.0), (50.0, 75.0), (75.0, 50.0)])
      nreal = 3

      method = SEQSIM(neighborhood=MetricBall(35.0), maxneighbors=3)

      # unconditional simulation
      rng = StableRNG(2017)
      real = rand(rng, proc, grid, nreal; method)
      @test eltype(real[1].Z) <: Float64

      # conditional simulation
      real = rand(rng, proc, grid, nreal; method, data)
      @test eltype(real[1].Z) <: Float64
      inds = LinearIndices(size(grid))
      @test all(real[i].Z[inds[25, 25]] == 1.0 for i in 1:nreal)
      @test all(real[i].Z[inds[50, 75]] == 0.0 for i in 1:nreal)
      @test all(real[i].Z[inds[75, 50]] == 1.0 for i in 1:nreal)
    end

    @testset "FFTSIM" begin
      # basic simulation
      rng = StableRNG(2019)
      proc = GaussianProcess(GaussianVariogram(range=10.0))
      grid = CartesianGrid(100, 100)
      real = rand(rng, proc, grid, 3, method=FFTSIM())
      @test eltype(real[1].Z) <: Float64

      # simulation on view of grid
      rng = StableRNG(2022)
      proc = GaussianProcess(GaussianVariogram(range=10.0))
      grid = CartesianGrid(100, 100)
      vgrid = view(grid, 1:5000)
      real = rand(rng, proc, vgrid, method=FFTSIM())
      @test domain(real) == vgrid
      @test length(real.geometry) == 5000

      # conditional simulation
      rng = StableRNG(2022)
      proc = GaussianProcess(GaussianVariogram(range=35.0))
      grid = CartesianGrid((100, 100), (0.5, 0.5), (1.0, 1.0))
      data = georef((; Z=[1.0, 0.0, 1.0]), [(25.0, 25.0), (50.0, 75.0), (75.0, 50.0)])
      real = rand(rng, proc, grid, method=FFTSIM(maxneighbors=3), data=data)
    end

    @testset "show" begin
      proc = GaussianProcess()
      @test sprint(show, proc) == "GaussianProcess(func: GaussianVariogram(range: 1.0 m, sill: 1.0, nugget: 0.0), mean: 0.0)"
      @test sprint(show, MIME("text/plain"), proc) == """
      GaussianProcess
      ├─ func: GaussianVariogram(range: 1.0 m, sill: 1.0, nugget: 0.0)
      └─ mean: 0.0"""
    end
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
    p = LindgrenProcess(0.1)
    s = Sphere((0, 0, 0))
    m = simplexify(s)
    # unconditional realization
    rng = StableRNG(2024)
    r = rand(rng, p, m, 3)
    for i in 1:3
      @test isapprox(sum(r[i].Z) / length(r[i].Z), 0.0, atol=1e-3)
    end
    # conditional realization
    rng = StableRNG(2024)
    d = georef((Z=[0.0, 1.0],), [(0, 0, -1), (0, 0, 1)])
    r = rand(rng, p, m, 3, data=d)
    for i in 1:3
      @test isapprox(sum(r[i].Z) / length(r[i].Z), 0.0, atol=1e-3)
    end

    proc = LindgrenProcess()
    @test sprint(show, proc) == "LindgrenProcess(range: 1.0 m, sill: 1.0)"
    @test sprint(show, MIME("text/plain"), proc) == """
    LindgrenProcess
    ├─ range: 1.0 m
    └─ sill: 1.0"""
  end

  @testset "QuiltingProcess" begin
    data = georef((; facies=[1.0, 0.0, 1.0]), [(25.0, 25.0), (50.0, 75.0), (75.0, 50.0)])
    grid = CartesianGrid(100, 100)

    rng = StableRNG(2017)
    trainimg = geostatsimage("Strebelle")
    inactive = [CartesianIndex(i, j) for i in 1:30 for j in 1:30]
    proc = QuiltingProcess(trainimg, (30, 30); inactive)

    real = rand(rng, proc, grid, 3; data)
    @test length(real) == 3
    @test size(domain(real[1])) == (100, 100)
    @test eltype(real[1].facies) <: Union{Float64,Missing}
    @test any(ismissing, real[1].facies)
    @test all(!isnan, skipmissing(real[1].facies))

    proc = QuiltingProcess(trainimg, (30, 30))
    @test sprint(show, proc) == "QuiltingProcess(trainimg: 62500×2 GeoTable over 250×250 CartesianGrid, tilesize: (30, 30), overlap: nothing, path: :raster, inactive: nothing, soft: nothing, tol: 0.1, nthreads: 8)"
    @test sprint(show, MIME("text/plain"), proc) == """
    QuiltingProcess
    ├─ trainimg: 62500×2 GeoTable over 250×250 CartesianGrid
    ├─ tilesize: (30, 30)
    ├─ overlap: nothing
    ├─ path: :raster
    ├─ inactive: nothing
    ├─ soft: nothing
    ├─ tol: 0.1
    └─ nthreads: 8"""
  end

  @testset "TuringProcess" begin
    rng = StableRNG(2019)
    grid = CartesianGrid(200, 200)
    real = rand(rng, TuringProcess(), grid, 3)
    @test length(real) == 3
    @test size(domain(real[1])) == (200, 200)
    @test eltype(real[1].Z) <: Float64

    proc = TuringProcess()
    @test sprint(show, proc) == "TuringProcess(params: nothing, blur: nothing, edge: nothing, iter: 100)"
    @test sprint(show, MIME("text/plain"), proc) == """
    TuringProcess
    ├─ params: nothing
    ├─ blur: nothing
    ├─ edge: nothing
    └─ iter: 100"""
  end

  @testset "StrataProcess" begin
    rng = StableRNG(2019)
    proc = SmoothingProcess()
    env = Environment(rng, [proc, proc], [0.5 0.5; 0.5 0.5], ExponentialDuration(rng, 1.0))
    grid = CartesianGrid(50, 50, 20)
    real = rand(StrataProcess(env), grid, 3)
    @test length(real) == 3
    @test size(domain(real[1])) == (50, 50, 20)
    @test eltype(real[1].Z) <: Union{Float64,Missing}
    @test any(ismissing, real[1].Z)
    @test all(!isnan, skipmissing(real[1].Z))
  end
end
