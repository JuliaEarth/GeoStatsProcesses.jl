@testset "FieldProcess" begin
  @testset "data argument" begin
    grid = CartesianGrid(10, 10)
    gtb = georef((; z=rand(100)), grid)
    # geotable
    setup = GeoStatsProcesses.randsetup(grid, gtb, 1)
    @test setup.geotable == gtb
    @test setup.varnames == [:z]
    @test setup.vartypes == [Float64]
    # pair
    setup = GeoStatsProcesses.randsetup(grid, :z => Float64, 1)
    @test isnothing(setup.geotable)
    @test setup.varnames == [:z]
    @test setup.vartypes == [Float64]
    setup = GeoStatsProcesses.randsetup(grid, "z" => Float64, 1)
    @test isnothing(setup.geotable)
    @test setup.varnames == [:z]
    @test setup.vartypes == [Float64]
    # iterable of pairs
    setup = GeoStatsProcesses.randsetup(grid, [:a => Float64, :b => Int], 1)
    @test isnothing(setup.geotable)
    @test setup.varnames == [:a, :b]
    @test setup.vartypes == [Float64, Int]
    setup = GeoStatsProcesses.randsetup(grid, ["a" => Float64, "b" => Int], 1)
    @test isnothing(setup.geotable)
    @test setup.varnames == [:a, :b]
    @test setup.vartypes == [Float64, Int]
    # error: invalid iterator
    @test_throws ArgumentError GeoStatsProcesses.randsetup(grid, [:a, :b], 1)
  end

  @testset "async" begin
    rng = StableRNG(2019)
    grid = CartesianGrid(100, 100)
    process = GaussianProcess()
    sims = rand(rng, process, grid, :z => Float64, 3, async=true)
    @test length(sims) == 3
    @test domain(sims[1]) == grid
    @test eltype(sims[1].z) <: Float64
    @test domain(sims[2]) == grid
    @test eltype(sims[2].z) <: Float64
    @test domain(sims[3]) == grid
    @test eltype(sims[3].z) <: Float64
    # error: passing the master worker when using the `async` option is not allowed
    @test_throws ArgumentError rand(rng, process, grid, :z => Float64, 3, workers=[myid()], async=true)
  end

  @testset "GaussianProcess" begin
    @testset "defaultmethod" begin
      grid = CartesianGrid(100, 100)
      vgrid = view(grid, 1:1000)
      pset1 = PointSet(rand(Point, 1000))
      pset2 = PointSet(rand(Point, 10000))

      process = GaussianProcess()
      setup = GeoStatsProcesses.randsetup(grid, :z => Float64, 1)
      @test GeoStatsProcesses.defaultmethod(process, setup) isa FFTMethod
      setup = GeoStatsProcesses.randsetup(vgrid, :z => Float64, 1)
      @test GeoStatsProcesses.defaultmethod(process, setup) isa FFTMethod
      setup = GeoStatsProcesses.randsetup(pset1, :z => Float64, 1)
      @test GeoStatsProcesses.defaultmethod(process, setup) isa LUMethod
      setup = GeoStatsProcesses.randsetup(pset2, :z => Float64, 1)
      @test GeoStatsProcesses.defaultmethod(process, setup) isa SEQMethod
    end

    @testset "FFTMethod" begin
      # isotropic simulation
      rng = StableRNG(2019)
      dom = CartesianGrid(100, 100)
      process = GaussianProcess(GaussianVariogram(range=10.0))
      method = FFTMethod()
      sims = rand(rng, process, dom, :z => Float64, 3, method)
      @test eltype(sims[1].z) <: Float64

      # anisotropic simulation
      rng = StableRNG(2019)
      dom = CartesianGrid(100, 100)
      process = GaussianProcess(GaussianVariogram(MetricBall((20.0, 5.0))))
      method = FFTMethod()
      sims = rand(rng, process, dom, :z => Float64, 3, method)

      # simulation on view of grid
      rng = StableRNG(2022)
      grid = CartesianGrid(100, 100)
      vgrid = view(grid, 1:5000)
      process = GaussianProcess(GaussianVariogram(range=10.0))
      method = FFTMethod()
      sim = rand(rng, process, vgrid, :z => Float64, method)
      @test domain(sim) == vgrid
      @test length(sim.geometry) == 5000

      # conditional simulation
      rng = StableRNG(2022)
      table = (; z=[1.0, -1.0, 1.0])
      coords = [(25.0, 25.0), (50.0, 75.0), (75.0, 50.0)]
      samples = georef(table, coords)
      sdomain = CartesianGrid(100, 100)
      process = GaussianProcess(GaussianVariogram(range=10.0))
      method = FFTMethod(maxneighbors=3)
      sim = rand(rng, process, sdomain, samples, method)
    end

    @testset "SEQMethod" begin
      𝒮 = georef((; z=[1.0, 0.0, 1.0]), [(25.0, 25.0), (50.0, 75.0), (75.0, 50.0)])
      𝒟 = CartesianGrid((100, 100), (0.5, 0.5), (1.0, 1.0))
      N = 3

      process = GaussianProcess(SphericalVariogram(range=35.0))
      method = SEQMethod(neighborhood=MetricBall(30.0), maxneighbors=3)

      rng = StableRNG(2017)
      sims₁ = rand(rng, process, 𝒟, 𝒮, 3)
      sims₂ = rand(rng, process, 𝒟, :z => Float64, 3)
      @test eltype(sims₁[1].z) <: Float64
      @test eltype(sims₂[1].z) <: Float64

      # basic checks
      inds = LinearIndices(size(𝒟))
      @test all(sims₁[i].z[inds[25, 25]] == 1.0 for i in 1:N)
      @test all(sims₁[i].z[inds[50, 75]] == 0.0 for i in 1:N)
      @test all(sims₁[i].z[inds[75, 50]] == 1.0 for i in 1:N)
    end

    @testset "LUMethod" begin
      𝒮 = georef((; z=[0.0, 1.0, 0.0, 1.0, 0.0]), [(0.0,), (25.0,), (50.0,), (75.0,), (100.0,)])
      𝒟 = CartesianGrid(100)

      # ----------------------
      # conditional simulation
      # ----------------------
      rng = StableRNG(123)
      process = GaussianProcess(SphericalVariogram(range=10.0))
      method = LUMethod()
      sims = rand(rng, process, 𝒟, 𝒮, 2, method)
      @test eltype(sims[1].z) <: Float64

      # ------------------------
      # unconditional simulation
      # ------------------------
      rng = StableRNG(123)
      process = GaussianProcess(SphericalVariogram(range=10.0))
      method = LUMethod()
      sims = rand(rng, process, 𝒟, :z => Float64, 2, method)

      # -------------
      # co-simulation
      # -------------
      rng = StableRNG(123)
      𝒟 = CartesianGrid(500)
      process = GaussianProcess((SphericalVariogram(range=10.0), GaussianVariogram(range=10.0)))
      method = LUMethod(correlation=0.95)
      sim = rand(rng, process, 𝒟, [:a => Float64, :b => Float64], method)

      # -----------
      # 2D example
      # -----------
      rng = StableRNG(123)
      𝒟 = CartesianGrid(100, 100)
      process = GaussianProcess(GaussianVariogram(range=10.0))
      method = LUMethod()
      sims = rand(rng, process, 𝒟, :z => Float64, 3, method)

      # -------------------
      # anisotropy example
      # -------------------
      rng = StableRNG(123)
      𝒟 = CartesianGrid(100, 100)
      ball = MetricBall((20.0, 5.0))
      process = GaussianProcess(GaussianVariogram(ball))
      method = LUMethod()
      sims = rand(rng, process, 𝒟, :z => Float64, 3, method)

      # ---------------------
      # custom factorization
      # ---------------------
      rng = StableRNG(123)
      𝒟 = CartesianGrid(100)
      process = GaussianProcess(SphericalVariogram(range=10.0))
      method1 = LUMethod(factorization=lu)
      method2 = LUMethod(factorization=cholesky)
      sim1 = rand(rng, process, 𝒟, 𝒮, 2, method1)
      sim2 = rand(rng, process, 𝒟, 𝒮, 2, method2)

      # throws
      𝒟 = CartesianGrid(100, 100)
      process = GaussianProcess(GaussianVariogram(range=10.0))
      method = LUMethod()
      # only 1 or 2 variables can be simulated simultaneously
      @test_throws AssertionError rand(process, 𝒟, [:a => Float64, :b => Float64, :c => Float64], method)
      process = GaussianProcess((GaussianVariogram(range=10.0),))
      # the number of parameters must be equal to the number of variables
      @test_throws AssertionError rand(process, 𝒟, [:a => Float64, :b => Float64], method)
    end

    @testset "show" begin
      process = GaussianProcess()
      @test sprint(show, process) == "GaussianProcess(variogram: GaussianVariogram(sill: 1.0, nugget: 0.0, range: 1.0 m, distance: Euclidean), mean: 0.0)"
      @test sprint(show, MIME("text/plain"), process) == """
      GaussianProcess
      ├─ variogram: GaussianVariogram(sill: 1.0, nugget: 0.0, range: 1.0 m, distance: Euclidean)
      └─ mean: 0.0"""
    end
  end

  @testset "LindgrenProcess" begin
    process = LindgrenProcess()
    @test process.range == 1u"m"
    @test process.sill == 1.0

    process = LindgrenProcess(2.0)
    @test process.range == 2.0u"m"
    @test process.sill == 1.0

    process = LindgrenProcess(2.0u"km")
    @test process.range == 2.0u"km"
    @test process.sill == 1.0

    process = LindgrenProcess(2.0, 2.0)
    @test process.range == 2.0u"m"
    @test process.sill == 2.0

    process = LindgrenProcess(2.0u"km", 2.0)
    @test process.range == 2.0u"km"
    @test process.sill == 2.0

    process = LindgrenProcess(2, 2)
    @test Unitful.numtype(process.range) == Float64
    @test process.sill isa Float64

    # simulation on sphere
    p = LindgrenProcess(0.1)
    s = Sphere((0,0,0))
    m = simplexify(s)
    # unconditional realization
    rng = StableRNG(2024)
    r = rand(rng, p, m, :z => Float64, 3)
    for i in 1:3
      @test isapprox(sum(r[i].z) / length(r[i].z), 0.0, atol=1e-3)
    end
    # conditional realization
    rng = StableRNG(2024)
    d = georef((z=[0.,1.],), [(0,0,-1), (0,0,1)])
    r = rand(rng, p, m, d, 3)
    for i in 1:3
      @test isapprox(sum(r[i].z) / length(r[i].z), 0.0, atol=1e-3)
    end

    process = LindgrenProcess()
    @test sprint(show, process) == "LindgrenProcess(range: 1.0 m, sill: 1.0, init: NearestInit())"
    @test sprint(show, MIME("text/plain"), process) == """
    LindgrenProcess
    ├─ range: 1.0 m
    ├─ sill: 1.0
    └─ init: NearestInit()"""
  end

  @testset "QuiltingProcess" begin
    sdata = georef((; facies=[1.0, 0.0, 1.0]), [(25.0, 25.0), (50.0, 75.0), (75.0, 50.0)])
    sdomain = CartesianGrid(100, 100)

    rng = StableRNG(2017)
    trainimg = geostatsimage("Strebelle")
    inactive = [CartesianIndex(i, j) for i in 1:30 for j in 1:30]
    process = QuiltingProcess(trainimg, (30, 30); inactive)

    sims = rand(rng, process, sdomain, sdata, 3)
    @test length(sims) == 3
    @test size(domain(sims[1])) == (100, 100)
    @test eltype(sims[1].facies) <: Union{Float64,Missing}
    @test any(ismissing, sims[1].facies)
    @test all(!isnan, skipmissing(sims[1].facies))

    process = QuiltingProcess(trainimg, (30, 30))
    @test sprint(show, process) == "QuiltingProcess(trainimg: 62500×2 GeoTable over 250×250 CartesianGrid, tilesize: (30, 30), overlap: nothing, path: :raster, inactive: nothing, soft: nothing, tol: 0.1, init: NearestInit())"
    @test sprint(show, MIME("text/plain"), process) == """
    QuiltingProcess
    ├─ trainimg: 62500×2 GeoTable over 250×250 CartesianGrid
    ├─ tilesize: (30, 30)
    ├─ overlap: nothing
    ├─ path: :raster
    ├─ inactive: nothing
    ├─ soft: nothing
    ├─ tol: 0.1
    └─ init: NearestInit()"""
  end

  @testset "TuringProcess" begin
    rng = StableRNG(2019)
    sdomain = CartesianGrid(200, 200)
    sims = rand(rng, TuringProcess(), sdomain, :z => Float64, 3)
    @test length(sims) == 3
    @test size(domain(sims[1])) == (200, 200)
    @test eltype(sims[1].z) <: Float64

    process = TuringProcess()
    @test sprint(show, process) == "TuringProcess(params: nothing, blur: nothing, edge: nothing, iter: 100)"
    @test sprint(show, MIME("text/plain"), process) == """
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
    sdomain = CartesianGrid(50, 50, 20)
    sims = rand(StrataProcess(env), sdomain, :z => Float64, 3)
    @test length(sims) == 3
    @test size(domain(sims[1])) == (50, 50, 20)
    @test eltype(sims[1].z) <: Union{Float64,Missing}
    @test any(ismissing, sims[1].z)
    @test all(!isnan, skipmissing(sims[1].z))

    rng = StableRNG(2019)
    proc = SmoothingProcess()
    env = Environment(rng, [proc, proc], [0.5 0.5; 0.5 0.5], ExponentialDuration(rng, 1.0))
    process = StrataProcess(env)
  end
end
