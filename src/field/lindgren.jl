# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    LindgrenProcess([paramaters])

The SPDE process introduced by Lindgren 2011. It relies on a
discretization of the Laplace-Beltrami operator on meshes and
is adequate for highly curved domains (e.g. surfaces).

## Parameters

* `range` - Range or correlation length (default to `1.0`)
* `sill`  - Sill or total variance (default to `1.0`)

### References

* Lindgren et al. 2011. [An explicit link between Gaussian fields and
  Gaussian Markov random fields: the stochastic partial differential
  equation approach](https://rss.onlinelibrary.wiley.com/doi/10.1111/j.1467-9868.2011.00777.x)
"""
@kwdef struct LindgrenProcess <: FieldProcess
  range::Float64 = 1.0
  sill::Float64 = 1.0
end

function randprep(::AbstractRNG, process::LindgrenProcess, ::DefaultRandMethod, setup::RandSetup)
  isnothing(setup.geotable) || @error "conditional simulation is not implemented"

  # retrieve sill and range
  𝓁 = process.range
  σ = process.sill

  @assert 𝓁 > 0 "range must be positive"
  @assert σ > 0 "sill must be positive"

  # retrieve domain info
  𝒟 = setup.domain
  d = paramdim(𝒟)

  # Beltrami-Laplace discretization
  B = laplacematrix(𝒟)
  M = measurematrix(𝒟)
  Δ = inv(M) * B

  # result of preprocessing
  pairs = map(setup.varnames) do var
    # LHS of SPDE (κ² - Δ)Z = τW
    α = 2one(σ + 𝓁)
    ν = α - d / 2
    κ = 1 / 𝓁
    A = κ^2 * I - Δ

    # covariance structure
    τ² = σ^2 * κ^(2ν) * (4π)^(d / 2) * gamma(α) / gamma(ν)
    Q = A'A / τ²

    # factorization
    F = cholesky(Array(Q))
    L = inv(Array(F.U))

    # save preprocessed inputs for variable
    var => (; L)
  end

  Dict(pairs)
end

function randsingle(rng::AbstractRNG, ::LindgrenProcess, ::DefaultRandMethod, setup::RandSetup, prep)
  # retrieve problem info
  𝒟 = setup.domain
  n = nvertices(𝒟)

  # simulation at vertices
  varreal = map(setup.vartypes, setup.varnames) do V, var
    # unpack preprocessed parameters
    L = prep[var].L

    # perform simulation
    w = randn(rng, V, n)
    z = L * w

    var => z
  end

  # vertex table
  vtable = (; varreal...)

  # change of support
  vdata = georef(vtable, 𝒟)
  edata = integrate(vdata, setup.varnames...)

  # columns of element table
  cols = Tables.columns(values(edata))
  Dict(var => Tables.getcolumn(cols, var) for var in setup.varnames)
end
