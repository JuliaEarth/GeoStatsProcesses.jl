@testset "Initialization" begin
  data1D = georef(CSV.File(joinpath(datadir, "data1D.tsv")), (:x,))
  data2D = georef(CSV.File(joinpath(datadir, "data2D.tsv")), (:x, :y))

  @testset "NearestInit" begin
    proc = GaussianProcess(GaussianVariogram())
    init = NearestInit()

    grid = CartesianGrid((100,), (-0.5,), (1.0,))
    real, mask = GeoStatsProcesses.randinit(proc, grid, data1D, init)
    for (i, j) in (100 => 11, 81 => 9, 11 => 2, 21 => 3, 91 => 10, 51 => 6, 61 => 7, 71 => 8, 31 => 4, 41 => 5, 1 => 1)
      @test real[:value][i] == data1D.value[j]
      @test mask[:value][i] == true
    end

    grid = CartesianGrid((100, 100), (-0.5, -0.5), (1.0, 1.0))
    real, mask = GeoStatsProcesses.randinit(proc, grid, data2D, init)
    for (i, j) in (5076 => 3, 2526 => 1, 7551 => 2)
      @test real[:value][i] == data2D.value[j]
      @test mask[:value][i] == true
    end

    pset = PointSet([(25, 25), (50, 75), (75, 50)])
    real, mask = GeoStatsProcesses.randinit(proc, pset, data2D, init)
    for (i, j) in (2 => 2, 3 => 3, 1 => 1)
      @test real[:value][i] == data2D.value[j]
      @test mask[:value][i] == true
    end
  end

  @testset "ExplicitInit" begin
    proc = GaussianProcess(GaussianVariogram())
    grid = CartesianGrid(10, 10, 10)
    data = georef((z=rand(10),), rand(Point, 10))

    # copy data to first locations in domain
    init = ExplicitInit(1:10)
    real, mask = GeoStatsProcesses.randinit(proc, grid, data, init)
    for (i, j) in (i => i for i in 1:10)
      @test real[:z][i] == data.z[j]
      @test mask[:z][i] == true
    end

    # copy data to last locations in domain
    init = ExplicitInit(991:1000)
    real, mask = GeoStatsProcesses.randinit(proc, grid, data, init)
    for (i, j) in (j => i for (i, j) in enumerate(991:1000))
      @test real[:z][i] == data.z[j]
      @test mask[:z][i] == true
    end

    # copy first 3 data points to last 3 domain locations
    init = ExplicitInit(1:3, 998:1000)
    real, mask = GeoStatsProcesses.randinit(proc, grid, data, init)
    for (i, j) in (j => i for (i, j) in enumerate(998:1000))
      @test real[:z][i] == data.z[j]
      @test mask[:z][i] == true
    end
  end
end
