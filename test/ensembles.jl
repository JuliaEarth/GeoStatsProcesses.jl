@testset "Ensembles" begin
  d = PointSet(rand(Point, 100))
  r = [(; z=1:100) for _ in 1:10]
  s = Ensemble(d, r)
  for i in 1:10
    @test s[i] == georef((z=1:100,), d)
  end

  @test sprint(show, s) == "3D Ensemble"
  @test sprint(show, MIME"text/plain"(), s) ==
        "3D Ensemble\n  domain:    100 PointSet\n  variables: z\n  N° reals:  10"

  d = CartesianGrid(10, 10)
  r = [(; z=1:100) for _ in 1:10]
  s = Ensemble(d, r)
  for i in 1:10
    @test s[i] == georef((z=1:100,), d)
  end

  @test sprint(show, s) == "2D Ensemble"
  @test sprint(show, MIME"text/plain"(), s) ==
        "2D Ensemble\n  domain:    10×10 CartesianGrid\n  variables: z\n  N° reals:  10"

  grid = CartesianGrid(3, 3)
  reals = [(; z=i * ones(nelements(grid))) for i in 1:3]
  ensemble = Ensemble(grid, reals)

  # mean
  mean2D = mean(ensemble)
  @test mean2D.z == 2.0 * ones(nrow(mean2D))

  # variance
  var2D = var(ensemble)
  @test var2D.z == 1.0 * ones(nrow(var2D))

  # cdf
  cdf2D = cdf(ensemble, 1)
  @test cdf2D.z == 1 / 3 * ones(nrow(cdf2D))
  cdf2D = cdf(ensemble, 2)
  @test cdf2D.z == 2 / 3 * ones(nrow(cdf2D))
  cdf2D = cdf(ensemble, 3)
  @test cdf2D.z == 3 / 3 * ones(nrow(cdf2D))

  # quantile (scalar)
  p = 0.5
  quant2D = quantile(ensemble, p)
  @test quant2D.z == 2.0 * ones(nrow(quant2D))

  # quantile (vector)
  ps = [0.0, 0.5, 1.0]
  quants2D = quantile(ensemble, ps)
  @test quants2D[2].z == quant2D.z
end
