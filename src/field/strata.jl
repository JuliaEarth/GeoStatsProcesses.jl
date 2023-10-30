# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    StrataProcess([paramaters])

Stratigraphy simulation with Markov-Poisson sampling.

## Parameters

* `environment` - Geological environment
* `state`       - Initial geological state
* `stack`       - Stacking scheme (:erosional or :depositional)
* `nepochs`     - Number of epochs (default to 10)
* `fillbase`    - Fill value for the bottom layer (default to `NaN`)
* `filltop`     - Fill value for the top layer (default to `NaN`)

### References

Hoffimann 2018. *Morphodynamic analysis and statistical
synthesis of geormorphic data.*
"""
@kwdef struct StrataProcess{E,S,ST,N,FB,FT} <: FieldProcess
  environment::E
  state::S = nothing
  stack::ST = :erosional
  nepochs::N = 10
  fillbase::FB = NaN
  filltop::FT = NaN
end
