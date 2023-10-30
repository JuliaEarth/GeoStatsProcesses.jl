@testset "Processes" begin
  @testset "FieldProcess" begin
    @testset "GaussianProcess" begin
      @testset "FFTMethod" begin
        # isotropic simulation
        Random.seed!(2019)
        dom = CartesianGrid(100, 100)
        process = GaussianProcess(variogram=GaussianVariogram(range=10.0))
        method = FFTMethod()
        sims = rand(process, dom, [:z => Float64], 3, method)
    
        # anisotropic simulation
        Random.seed!(2019)
        dom = CartesianGrid(100, 100)
        process = GaussianProcess(variogram=GaussianVariogram(MetricBall((20.0, 5.0))))
        method = FFTMethod()
        sims = rand(process, dom, [:z => Float64], 3, method)
    
        # simulation on view of grid
        Random.seed!(2022)
        grid = CartesianGrid(100, 100)
        vgrid = view(grid, 1:5000)
        process = GaussianProcess(variogram=GaussianVariogram(range=10.0))
        method = FFTMethod()
        sim = rand(process, vgrid, [:z => Float64], method)
        @test domain(sim) == vgrid
        @test length(sim.geometry) == 5000
    
        # conditional simulation
        Random.seed!(2022)
        table = (; z=[1.0, -1.0, 1.0])
        coords = [(25.0, 25.0), (50.0, 75.0), (75.0, 50.0)]
        samples = georef(table, coords)
        sdomain = CartesianGrid(100, 100)
        process = GaussianProcess(variogram=GaussianVariogram(range=10.0))
        method = FFTMethod(maxneighbors=3)
        sim = rand(process, sdomain, samples, method)
      end
    
      @testset "SEQMethod" begin
        𝒮 = georef((; z=[1.0, 0.0, 1.0]), [25.0 50.0 75.0; 25.0 75.0 50.0])
        𝒟 = CartesianGrid((100, 100), (0.5, 0.5), (1.0, 1.0))
        N = 3
    
        process = GaussianProcess(variogram=SphericalVariogram(range=35.0))
        method = SEQMethod(neighborhood=MetricBall(30.0), maxneighbors=3)
    
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
    
      @testset "LUMethod" begin
        𝒮 = georef((; z=[0.0, 1.0, 0.0, 1.0, 0.0]), [0.0 25.0 50.0 75.0 100.0])
        𝒟 = CartesianGrid(100)
    
        # ----------------------
        # conditional simulation
        # ----------------------
        rng = MersenneTwister(123)
        process = GaussianProcess(variogram=SphericalVariogram(range=10.0))
        method = LUMethod()
        sims = rand(rng, process, 𝒟, 𝒮, 2, method)
    
        # ------------------------
        # unconditional simulation
        # ------------------------
        rng = MersenneTwister(123)
        process = GaussianProcess(variogram=SphericalVariogram(range=10.0))
        method = LUMethod()
        sims = rand(rng, process, 𝒟, [:z => Float64], 2, method)
    
        # -------------
        # co-simulation
        # -------------
        𝒟 = CartesianGrid(500)
        rng = MersenneTwister(123)
        process = GaussianProcess(variogram=(SphericalVariogram(range=10.0), GaussianVariogram(range=10.0)))
        method = LUMethod(correlation=0.95)
        sim = rand(rng, process, 𝒟, [:a => Float64, :b => Float64], method)
    
        # -----------
        # 2D example
        # -----------
        𝒟 = CartesianGrid(100, 100)
        rng = MersenneTwister(123)
        process = GaussianProcess(variogram=GaussianVariogram(range=10.0))
        method = LUMethod()
        sims = rand(rng, process, 𝒟, [:z => Float64], 3, method)
    
        # -------------------
        # anisotropy example
        # -------------------
        𝒟 = CartesianGrid(100, 100)
        rng = MersenneTwister(123)
        ball = MetricBall((20.0, 5.0))
        process = GaussianProcess(variogram=GaussianVariogram(ball))
        method = LUMethod()
        sims = rand(rng, process, 𝒟, [:z => Float64], 3, method)
    
        # ---------------------
        # custom factorization
        # ---------------------
        𝒟 = CartesianGrid(100)
        rng = MersenneTwister(123)
        process = GaussianProcess(variogram=SphericalVariogram(range=10.0))
        method1 = LUMethod(factorization=lu)
        method2 = LUMethod(factorization=cholesky)
        sim1 = rand(rng, process, 𝒟, 𝒮, 2, method1)
        sim2 = rand(rng, process, 𝒟, 𝒮, 2, method2)
    
        # throws
        𝒟 = CartesianGrid(100, 100)
        process = GaussianProcess(variogram=GaussianVariogram(range=10.0))
        method = LUMethod()
        # only 1 or 2 variables can be simulated simultaneously
        @test_throws AssertionError rand(process, 𝒟, [:a => Float64, :b => Float64, :c => Float64], method)
        process = GaussianProcess(variogram=(GaussianVariogram(range=10.0),))
        # the number of parameters must be equal to the number of variables
        @test_throws AssertionError rand(process, 𝒟, [:a => Float64, :b => Float64], method)
      end
    end
  
    @testset "QuiltingProcess" begin
      sdata = georef((; facies=[1.0, 0.0, 1.0]), [25.0 50.0 75.0; 25.0 75.0 50.0])
      sdomain = CartesianGrid(100, 100)
    
      rng = MersenneTwister(2017)
      trainimg = geostatsimage("Strebelle")
      inactive = [CartesianIndex(i, j) for i in 1:30 for j in 1:30]
      process = QuiltingProcess(trainimg=trainimg, tilesize=(30, 30), inactive=inactive)
    
      sims = rand(rng, process, sdomain, sdata, 3)
      @test length(sims) == 3
      @test size(domain(sims[1])) == (100, 100)
    end
  
    @testset "TuringProcess" begin
      Random.seed!(2019)
      sdomain = CartesianGrid(200, 200)
      sims = rand(TuringProcess(), sdomain, [:z => Float64], 3)
      @test length(sims) == 3
      @test size(domain(sims[1])) == (200, 200)
    end
  
    @testset "StrataProcess" begin
      rng = MersenneTwister(2019)
      proc = SmoothingProcess()
      env = Environment(rng, [proc, proc], [0.5 0.5; 0.5 0.5], ExponentialDuration(rng, 1.0))
      sdomain = CartesianGrid(50, 50, 20)
      sims = rand(StrataProcess(environment=env), sdomain, [:z => Float64], 3)
      @test length(sims) == 3
      @test size(domain(sims[1])) == (50, 50, 20)
    end
  end

  @testset "PointProcess" begin
    # geometries and domains
    seg = Segment((0.0, 0.0), (11.3, 11.3))
    tri = Triangle((0.0, 0.0), (5.65, 0.0), (5.65, 5.65))
    quad = Quadrangle((0.0, 0.0), (0.0, 4.0), (4.0, 4.0), (4.0, 0.0))
    box = Box((0.0, 0.0), (4.0, 4.0))
    ball = Ball((1.0, 1.0), 2.25)
    outer = [(0, -4), (4, -1), (4, 1.5), (0, 3)]
    hole1 = [(0.2, -0.2), (1.4, -0.2), (1.4, 0.6), (0.2, 0.6)]
    hole2 = [(2, -0.2), (3, -0.2), (3, 0.4), (2, 0.4)]
    poly = PolyArea([outer, hole1, hole2])
    grid = CartesianGrid((0, 0), (4, 4), dims=(10, 10))
    points = Point2[(0, 0), (4.5, 0), (0, 4.2), (4, 4.3), (1.5, 1.5)]
    connec = connect.([(1, 2, 5), (2, 4, 5), (4, 3, 5), (3, 1, 5)], Triangle)
    mesh = SimpleMesh(points, connec)
    geoms = [seg, tri, quad, box, ball, poly, grid, mesh]
  
    # point processes
    λ(s) = sum(coordinates(s) .^ 2)
    binom = BinomialProcess(100)
    poisson1 = PoissonProcess(100.0)
    poisson2 = PoissonProcess(λ)
    inhibit = InhibitionProcess(0.1)
    procs = [binom, poisson1, poisson2, inhibit]
  
    @testset "Basic" begin
      for p in procs, g in geoms
        pp = rand(p, g)
        @test all(∈(g), pp)
      end
    end
  
    @testset "Binomial" begin
      p = BinomialProcess(10)
      for g in geoms
        pp = rand(p, g)
        @test nelements(pp) == 10
      end
    end
  
    @testset "Poisson" begin
      # inhomogeneous with piecewise constant intensity
      for g in [grid, mesh]
        p = PoissonProcess(λ.(centroid.(g)))
        pp = rand(p, g)
        @test all(∈(g), pp)
      end
  
      # empty pointsets
      for g in geoms
        @test isnothing(rand(PoissonProcess(0.0), seg))
      end
  
      pp = PointSet(rand(Point2, 10))
      @test isnothing(rand(PoissonProcess(100.0), pp))
    end
  
    @testset "Inhibition" begin end
  
    @testset "Cluster" begin
      binom = BinomialProcess(100)
      poisson = PoissonProcess(100.0)
      procs = [binom, poisson]
  
      ofun1 = parent -> rand(BinomialProcess(10), Ball(parent, 0.2))
      ofun2 = parent -> rand(PoissonProcess(100), Ball(parent, 0.2))
      ofun3 = parent -> rand(PoissonProcess(x -> 100 * sum((x - parent) .^ 2)), Ball(parent, 0.5))
      ofun4 = parent -> PointSet(sample(Sphere(parent, 0.1), RegularSampling(10)))
      ofuns = [ofun1, ofun2, ofun3, ofun4]
  
      box = Box((0.0, 0.0), (4.0, 4.0))
      ball = Ball((1.0, 1.0), 2.25)
      tri = Triangle((0.0, 0.0), (5.65, 0.0), (5.65, 5.65))
      grid = CartesianGrid((0, 0), (4, 4), dims=(10, 10))
      geoms = [box, ball, tri, grid]
  
      for p in procs, ofun in ofuns, g in geoms
        cp = ClusterProcess(p, ofun)
        pp = rand(cp, g)
        if !isnothing(pp)
          @test all(∈(g), pp)
        end
      end
    end
  
    @testset "Union" begin
      b = Box((0.0, 0.0), (100.0, 100.0))
      p₁ = BinomialProcess(50)
      p₂ = BinomialProcess(50)
      p = p₁ ∪ p₂ # 100 points
  
      s = rand(p, b, 2)
      @test length(s) == 2
      @test s[1] isa PointSet
      @test s[2] isa PointSet
      @test nelements.(s) == [100, 100]
    end
  end  
end