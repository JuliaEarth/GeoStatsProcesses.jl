# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    StrataProcess(environment; [options])

Strata process with given geological `environment`
as described in Hoffimann 2018.

## Options

* `state`       - Initial geological state
* `stack`       - Stacking scheme (:erosional or :depositional)
* `nepochs`     - Number of epochs (default to 10)

## References

* Hoffimann 2018. [Morphodynamic analysis and statistical synthesis of
  geormorphic data](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2019JF005245)
"""
struct StrataProcess{E,S,ST,N} <: FieldProcess
  environment::E
  state::S
  stack::ST
  nepochs::N
end

StrataProcess(environment; state=nothing, stack=:erosional, nepochs=10) =
  StrataProcess(environment, state, stack, nepochs)
