# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Ensemble

An ensemble of geostatistical realizations from a geostatistical process.
"""
struct Ensemble{D,V,R,F}
  domain::D
  vars::V
  reals::R
  fetch::F
end

Ensemble(domain, vars, reals; fetch=identity) = Ensemble(domain, vars, reals, fetch)

# -------------
# ITERATOR API
# -------------

Base.iterate(e::Ensemble, state=1) = state > length(e) ? nothing : (e[state], state + 1)
Base.length(e::Ensemble) = length(e.reals)

# --------------
# INDEXABLE API
# --------------

function Base.getindex(e::Ensemble, ind::Int)
  domain = e.domain
  real = e.fetch(e.reals[ind])
  table = (; (var => real[var] for var in e.vars)...)
  georef(table, domain)
end
Base.getindex(e::Ensemble, inds::AbstractVector{Int}) = [getindex(e, ind) for ind in inds]
Base.firstindex(e::Ensemble) = 1
Base.lastindex(e::Ensemble) = length(e)

# -----------
# STATISTICS
# -----------

mean(e::Ensemble) = ereduce(mean, e)

var(e::Ensemble) = ereduce(var, e)

cdf(e::Ensemble, x::Number) = ereduce(vals -> count(≤(x), vals) / length(vals), e)

quantile(e::Ensemble, p::Number) = ereduce(vals -> quantile(vals, p), e)

quantile(e::Ensemble, ps::AbstractVector) = [quantile(e, p) for p in ps]

# -----------------
# HELPER FUNCTIONS
# -----------------

function ereduce(f, e)
  function reducevar(var)
    map(1:nelements(e.domain)) do i
      rowvals = (e.fetch(real)[var][i] for real in e.reals)
      f(rowvals)
    end
  end
  table = (; (var => reducevar(var) for var in e.vars)...)
  georef(table, e.domain)
end

# -----------
# IO METHODS
# -----------

function Base.show(io::IO, e::Ensemble)
  N = embeddim(e.domain)
  print(io, "$(N)D Ensemble")
end

function Base.show(io::IO, ::MIME"text/plain", e::Ensemble)
  println(io, e)
  println(io, "  domain:    ", e.domain)
  println(io, "  variables: ", join(e.vars, ", ", " and "))
  print(io, "  N° reals:  ", length(e))
end
