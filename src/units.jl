# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

const Len{T} = Quantity{T,u"𝐋"}

addunit(x::Number, u) = x * u
addunit(x::Quantity, _) = x
