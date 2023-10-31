# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    TuringProcess([paramaters])

Turing process as described in Turing 1952.

## Parameters

* `params` - basic parameters (default to `nothing`)
* `blur`   - blur algorithm (default to `nothing`)
* `edge`   - edge condition (default to `nothing`)
* `iter`   - number of iterations (default to `100`)

### References

Turing 1952. *The chemical basis of morphogenesis.*
"""
struct TuringProcess{P,B,E,I} <: FieldProcess
  params::P
  blur::B
  edge::E
  iter::I

  function TuringProcess(params::P=nothing, blur::B=nothing, edge::E=nothing, iter::I=100) where {P,B,E,I}
    new{P,B,E,I}(params, blur, edge, iter)
  end
end

TuringProcess(; params=nothing, blur=nothing, edge=nothing, iter=100) = TuringProcess(params, blur, edge, iter)
