# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    ExplicitInit([orig], dest)

A method to initialize realizations using explicit lists of
indices `orig` and `dest` in the data and domain, respectively.
If `orig` is ommitted, use all the indices of the `data`.
"""
struct ExplicitInit{V1,V2} <: InitMethod
  orig::V1
  dest::V2
end

ExplicitInit(dest) = ExplicitInit(nothing, dest)

function initialize!(real, mask, dom, data, method::ExplicitInit)
  # auxiliary variables
  ddom = domain(data)
  cols = Tables.columns(values(data))

  # retrieve origin and destination indices
  orig = isnothing(method.orig) ? (1:nelements(ddom)) : method.orig
  dest = method.dest

  # sanity checks
  @assert dest isa AbstractVector "destination indices must be provided as a vector"
  @assert length(orig) == length(dest) "invalid explicit initialization"

  # loop over pairs of locations
  @inbounds for ind in eachindex(orig, dest)
    # extract specific pair
    i, j = orig[ind], dest[ind]

    # update realization and mask if not missing
    for var in keys(real)
      valᵢ = Tables.getcolumn(cols, var)[i]
      if !ismissing(valᵢ)
        real[var][j] = valᵢ
        mask[var][j] = true
      end
    end
  end
end
