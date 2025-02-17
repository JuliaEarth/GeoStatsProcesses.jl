# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    NearestInit()

A method to initialize realizations using the nearest sample in the domain.
"""
struct NearestInit <: InitMethod end

function initialize!(real, mask, dom, data, ::NearestInit)
  # initialize nearest search
  s = KNearestSearch(dom, 1)

  # auxiliary variables
  ddom = domain(data)
  cols = Tables.columns(values(data))

  # loop over data locations
  @inbounds for i in 1:nelements(ddom)
    # find nearest location in simulation domain
    j = search(centroid(ddom, i), s)[1]

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
