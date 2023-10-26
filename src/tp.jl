# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    TP([paramaters])

Turing pattern simulation.

## Parameters

* `params` - basic parameters (default to `nothing`)
* `blur`   - blur algorithm (default to `nothing`)
* `edge`   - edge condition (default to `nothing`)
* `iter`   - number of iterations (default to `100`)

### References

Turing 1952. *The chemical basis of morphogenesis.*
"""
@kwdef struct TP{P,B,E,I} <: FieldProcess
  params::P = nothing
  blur::B = nothing
  edge::E = nothing
  iter::I = 100
end
