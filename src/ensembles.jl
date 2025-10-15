# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Ensemble

An ensemble of geostatistical realizations from a geostatistical process.
"""
struct Ensemble{D,R,F}
  domain::D
  reals::R
  fetch::F
end

Ensemble(domain, reals; fetch=identity) = Ensemble(domain, reals, fetch)

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
  table = e.fetch(e.reals[ind])
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

ccdf(e::Ensemble, x::Number) = ereduce(vals -> count(>(x), vals) / length(vals), e)

quantile(e::Ensemble, p::Number) = ereduce(vals -> quantile(vals, p), e)

quantile(e::Ensemble, ps::AbstractVector) = [quantile(e, p) for p in ps]

# -----------
# IO METHODS
# -----------

function Base.show(io::IO, e::Ensemble)
  dim = embeddim(e.domain)
  print(io, "$(dim)D Ensemble")
end

function Base.show(io::IO, ::MIME"text/plain", e::Ensemble)
  println(io, e)
  println(io, "  domain:    ", e.domain)
  println(io, "  variables: ", join(evars(e), ", ", " and "))
  print(io, "  N° reals:  ", length(e))
end

# -----------------
# HELPER FUNCTIONS
# -----------------

evars(e) = e.reals |> first |> Tables.columns |> Tables.columnnames

function ereduce(f, e)
  function reducevar(var)
    map(1:nelements(e.domain)) do i
      vals = (e.fetch(real)[var][i] for real in e.reals)
      f(vals)
    end
  end
  table = (; (var => reducevar(var) for var in evars(e))...)
  georef(table, e.domain)
end
