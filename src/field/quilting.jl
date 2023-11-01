# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    QuiltingProcess(trainimg, tilesize; [paramaters])

Image quilting process as described in Hoffimann et al. 2017.

## Parameters

### Required

* `trainimg` - Training image from which to extract tiles
* `tilesize` - Tuple with tile size for each dimension

### Optional

* `overlap`  - Overlap size (default to (1/6, 1/6, ..., 1/6))
* `path`     - Simulation path (`:raster` (default), `:dilation`, or `:random`)
* `inactive` - Vector of inactive voxels (i.e. `CartesianIndex`) in the grid
* `soft`     - A pair `(data,dataTI)` of geospatial data objects (default to `nothing`)
* `tol`      - Initial relaxation tolerance in (0,1] (default to `0.1`)
* `init`     - Data initialization method (default to `NearestInit()`)

## References

* Hoffimann et al 2017. *Stochastic simulation by image quilting of process-based geological models.*
* Hoffimann et al 2015. *Geostatistical modeling of evolving landscapes by means of image quilting.*
"""
struct QuiltingProcess{TR,TS,O,P,IN,S,T,I} <: FieldProcess
  trainimg::TR
  tilesize::TS
  overlap::O
  path::P
  inactive::IN
  soft::S
  tol::T
  init::I
end

QuiltingProcess(
  trainimg,
  tilesize;
  overlap=nothing,
  path=:raster,
  inactive=nothing,
  soft=nothing,
  tol=0.1,
  init=NearestInit()
) = QuiltingProcess(trainimg, tilesize, overlap, path, inactive, soft, tol, init)
