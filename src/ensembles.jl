# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Ensemble

An ensemble of geostatistical realizations from a geostatistical process.
"""
struct Ensemble{D,V,T,R,F}
  domain::D
  varnames::V
  vartypes::T
  reals::R
  nreals::Int
  fetch::F
end

Ensemble(domain, varnames, vartypes, reals; fetch=identity) = Ensemble(domain, varnames, vartypes, reals, length(reals), fetch)

# -------------
# ITERATOR API
# -------------

Base.iterate(e::Ensemble, state=1) = state > e.nreals ? nothing : (e[state], state + 1)
Base.length(e::Ensemble) = e.nreals

# --------------
# INDEXABLE API
# --------------

function Base.getindex(e::Ensemble, ind::Int)
  domain = e.domain
  real = e.fetch(e.reals[ind])
  table = (; (var => real[var] for var in e.varnames)...)
  georef(table, domain)
end
Base.getindex(e::Ensemble, inds::AbstractVector{Int}) = [getindex(e, ind) for ind in inds]
Base.firstindex(e::Ensemble) = 1
Base.lastindex(e::Ensemble) = length(e)

# -----------
# STATISTICS
# -----------

function mean(e::Ensemble)
  reals = map(e.fetch, e.reals)
  calc(var) = mean(real[var] for real in reals)
  mtable = (; (var => calc(var) for var in e.varnames)...)
  georef(mtable, e.domain)
end

function var(e::Ensemble)
  reals = map(e.fetch, e.reals)
  calc(v) = var([real[v] for real in reals])
  vtable = (; (v => calc(v) for v in e.varnames)...)
  georef(vtable, e.domain)
end

function cdf(e::Ensemble, x::Number)
  reals = map(e.fetch, e.reals)
  function calc(var)
    map(1:nelements(e.domain)) do ind
      slice = [real[var][ind] for real in reals]
      _cdf(slice, x)
    end
  end
  ctable = (; (var => calc(var) for var in e.varnames)...)
  georef(ctable, e.domain)
end

function quantile(e::Ensemble, p::Number)
  reals = map(e.fetch, e.reals)
  function calc(var)
    map(1:nelements(e.domain)) do ind
      slice = [real[var][ind] for real in reals]
      quantile(slice, p)
    end
  end
  qtable = (; (var => calc(var) for var in e.varnames)...)
  georef(qtable, e.domain)
end

quantile(e::Ensemble, ps::AbstractVector) = [quantile(e, p) for p in ps]

# -----------------
# HELPER FUNCTIONS
# -----------------

_cdf(xs, x) = count(≤(x), xs) / length(xs)

# -----------
# IO METHODS
# -----------

function Base.show(io::IO, e::Ensemble)
  N = embeddim(e.domain)
  print(io, "$(N)D Ensemble")
end

function Base.show(io::IO, ::MIME"text/plain", e::Ensemble)
  vars = ["$n ($t)" for (n, t) in zip(e.varnames, e.vartypes)]
  println(io, e)
  println(io, "  domain:    ", e.domain)
  println(io, "  variables: ", join(vars, ", ", " and "))
  print(io, "  N° reals:  ", e.nreals)
end
