# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    TuringProcess(; [options])

Turing process as described in Turing 1952.

## Options

* `params` - basic parameters (default to `nothing`)
* `blur`   - blur algorithm (default to `nothing`)
* `edge`   - edge condition (default to `nothing`)
* `iter`   - number of iterations (default to `100`)

## References

* Turing 1952. [The chemical basis of morphogenesis](https://royalsocietypublishing.org/doi/10.1098/rstb.1952.0012)
"""
struct TuringProcess{P,B,E,I} <: FieldProcess
  params::P
  blur::B
  edge::E
  iter::I
end

TuringProcess(; params=nothing, blur=nothing, edge=nothing, iter=100) = TuringProcess(params, blur, edge, iter)
