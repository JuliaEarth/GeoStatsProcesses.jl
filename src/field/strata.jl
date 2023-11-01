# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    StrataProcess(environment; [paramaters])

Strata process with given geological `environment`
as described in Hoffimann 2018.

## Parameters

* `state`       - Initial geological state
* `stack`       - Stacking scheme (:erosional or :depositional)
* `nepochs`     - Number of epochs (default to 10)
* `fillbase`    - Fill value for the bottom layer (default to `NaN`)
* `filltop`     - Fill value for the top layer (default to `NaN`)

### References

Hoffimann 2018. *Morphodynamic analysis and statistical
synthesis of geormorphic data.*
"""
struct StrataProcess{E,S,ST,N,FB,FT} <: FieldProcess
  environment::E
  state::S
  stack::ST
  nepochs::N
  fillbase::FB
  filltop::FT
end

StrataProcess(environment; state=nothing, stack=:erosional, nepochs=10, fillbase=NaN, filltop=NaN) =
  StrataProcess(environment, state, stack, nepochs, fillbase, filltop)
