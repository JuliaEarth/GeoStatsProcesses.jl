@testset "PointProcess" begin
  # geometries and domains
  seg = Segment((0.0, 0.0), (11.3, 11.3))
  tri = Triangle((0.0, 0.0), (5.65, 0.0), (5.65, 5.65))
  quad = Quadrangle((0.0, 0.0), (0.0, 4.0), (4.0, 4.0), (4.0, 0.0))
  box = Box((0.0, 0.0), (4.0, 4.0))
  ball = Ball((1.0, 1.0), 2.25)
  outer = [(0.0, -4.0), (4.0, -1.0), (4.0, 1.5), (0.0, 3.0)]
  hole1 = [(0.2, 0.6), (1.4, 0.6), (1.4, -0.2), (0.2, -0.2)]
  hole2 = [(2.0, 0.4), (3.0, 0.4), (3.0, -0.2), (2.0, -0.2)]
  poly = PolyArea([outer, hole1, hole2])
  grid = CartesianGrid((0, 0), (4, 4), dims=(10, 10))
  points = [(0, 0), (4.5, 0), (0, 4.2), (4, 4.3), (1.5, 1.5)]
  connec = connect.([(1, 2, 5), (2, 4, 5), (4, 3, 5), (3, 1, 5)], Triangle)
  mesh = SimpleMesh(points, connec)
  geoms = [seg, tri, quad, box, ball, poly, grid, mesh]

  # point processes
  λ(s) = sum(to(s) .^ 2)
  binom = BinomialProcess(100)
  poisson1 = PoissonProcess(100.0)
  poisson2 = PoissonProcess(λ)
  inhibit = InhibitionProcess(0.1)
  procs = [binom, poisson1, poisson2, inhibit]

  @testset "Basic" begin
    rng = StableRNG(2024)

    for p in procs
      pp = rand(rng, p, box)
      @test all(∈(box), pp)
    end

    for g in geoms
      pp = rand(rng, binom, g)
      @test all(∈(g), pp)
    end
  end

  @testset "Binomial" begin
    rng = StableRNG(2024)

    p = BinomialProcess(10)
    for g in geoms
      pp = rand(rng, p, g)
      @test nelements(pp) == 10
    end

    p = BinomialProcess(10)
    @test sprint(show, p) == "BinomialProcess(n: 10)"
    @test sprint(show, MIME("text/plain"), p) == """
    BinomialProcess
    └─ n: 10"""
  end

  @testset "Poisson" begin
    rng = StableRNG(2024)

    # inhomogeneous with piecewise constant intensity
    for g in [grid, mesh]
      p = PoissonProcess(λ.(centroid.(g)))
      pp = rand(rng, p, g)
      @test all(∈(g), pp)
    end

    # empty pointsets
    for g in geoms
      @test isnothing(rand(rng, PoissonProcess(0.0), seg))
    end

    pp = PointSet(rand(rng, Point, 10))
    @test isnothing(rand(rng, PoissonProcess(100.0), pp))

    p = PoissonProcess(100.0)
    @test sprint(show, p) == "PoissonProcess(λ: 100.0)"
    @test sprint(show, MIME("text/plain"), p) == """
    PoissonProcess
    └─ λ: 100.0"""
  end

  @testset "Inhibition" begin
    p = InhibitionProcess(0.1)
    @test sprint(show, p) == "InhibitionProcess(δ: 0.1)"
    @test sprint(show, MIME("text/plain"), p) == """
    InhibitionProcess
    └─ δ: 0.1"""
  end

  @testset "Cluster" begin
    rng = StableRNG(2024)

    binom = BinomialProcess(100)
    poisson = PoissonProcess(100.0)
    procs = [binom, poisson]

    ofun1 = parent -> rand(rng, BinomialProcess(10), Ball(parent, 0.2))
    ofun2 = parent -> rand(rng, PoissonProcess(100), Ball(parent, 0.2))
    ofun3 = parent -> rand(rng, PoissonProcess(x -> 100 * sum((x - parent) .^ 2)), Ball(parent, 0.5))
    ofun4 = parent -> PointSet(sample(Sphere(parent, 0.1), RegularSampling(10)))
    ofuns = [ofun1, ofun2, ofun3, ofun4]

    box = Box((0.0, 0.0), (4.0, 4.0))
    ball = Ball((1.0, 1.0), 2.25)
    tri = Triangle((0.0, 0.0), (5.65, 0.0), (5.65, 5.65))
    grid = CartesianGrid((0, 0), (4, 4), dims=(10, 10))
    geoms = [box, ball, tri, grid]

    for p in procs, ofun in ofuns
      cp = ClusterProcess(p, ofun)
      pp = rand(rng, cp, box)
      if !isnothing(pp)
        @test all(∈(box), pp)
      end
    end

    for g in geoms
      cp = ClusterProcess(binom, ofun1)
      pp = rand(rng, cp, g)
      if !isnothing(pp)
        @test all(∈(g), pp)
      end
    end

    p = ClusterProcess(binom, identity)
    @test sprint(show, p) == "ClusterProcess(proc: BinomialProcess(n: 100), ofun: identity)"
    @test sprint(show, MIME("text/plain"), p) == """
    ClusterProcess
    ├─ proc: BinomialProcess(n: 100)
    └─ ofun: identity"""
  end

  @testset "Union" begin
    rng = StableRNG(2024)

    b = Box((0.0, 0.0), (100.0, 100.0))
    p₁ = BinomialProcess(50)
    p₂ = BinomialProcess(50)
    p = p₁ ∪ p₂ # 100 points

    s = rand(rng, p, b, 2)
    @test length(s) == 2
    @test s[1] isa PointSet
    @test s[2] isa PointSet
    @test nelements.(s) == [100, 100]

    p = BinomialProcess(50) ∪ BinomialProcess(100)
    @test sprint(show, p) == "UnionProcess(p₁: BinomialProcess(n: 50), p₂: BinomialProcess(n: 100))"
    @test sprint(show, MIME("text/plain"), p) == """
    UnionProcess
    ├─ p₁: BinomialProcess(n: 50)
    └─ p₂: BinomialProcess(n: 100)"""
  end

  @testset "Thinning" begin
    rng = StableRNG(2024)

    p = PoissonProcess(10)
    q = Quadrangle((0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0))
    pp = rand(rng, p, q)

    tp = thin(pp, RandomThinning(0.3))
    @test length(tp) ≤ length(pp)
    xs = to.(tp)
    @test all(0u"m" .≤ first.(xs) .≤ 4.0u"m")
    @test all(0u"m" .≤ last.(xs) .≤ 4.0u"m")

    tp = thin(pp, RandomThinning(s -> λ(s) / λ(Point(4.0, 4.0))))
    @test length(tp) ≤ length(pp)
    xs = to.(tp)
    @test all(0u"m" .≤ first.(xs) .≤ 4.0u"m")
    @test all(0u"m" .≤ last.(xs) .≤ 4.0u"m")
  end
end
