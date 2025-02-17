# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    InhibitionProcess(δ)

An inhibition point process with minimum distance `δ`.
"""
struct InhibitionProcess{D<:Real} <: PointProcess
  δ::D
end
