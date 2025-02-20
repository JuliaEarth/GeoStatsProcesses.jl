# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    QuiltingProcess(trainimg, tilesize; [options])

Image quilting process with training image `trainimg` and tile size `tilesize`
as described in Hoffimann et al. 2017.

## Options

* `overlap`  - Overlap size (default to (1/6, 1/6, ..., 1/6))
* `path`     - Process path (`:raster` (default), `:dilation`, or `:random`)
* `soft`     - A pair `(data,dataTI)` of geospatial data objects (default to `nothing`)
* `tol`      - Initial relaxation tolerance in (0,1] (default to `0.1`)

## References

* Hoffimann et al 2017. [Stochastic simulation by image quilting
  of process-based geological models](https://www.sciencedirect.com/science/article/abs/pii/S0098300417301139)

* Hoffimann et al 2015. [Geostatistical modeling of evolving landscapes
  by means of image quilting](https://www.researchgate.net/publication/295902985_Geostatistical_Modeling_of_Evolving_Landscapes_by_Means_of_Image_Quilting)
"""
struct QuiltingProcess{TR,TS,O,P,S,T} <: FieldProcess
  trainimg::TR
  tilesize::TS
  overlap::O
  path::P
  soft::S
  tol::T
end

QuiltingProcess(
  trainimg,
  tilesize;
  overlap=nothing,
  path=:raster,
  soft=nothing,
  tol=0.1
) = QuiltingProcess(trainimg, tilesize, overlap, path, soft, tol)
