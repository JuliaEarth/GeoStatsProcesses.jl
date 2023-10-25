# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

module GeoStatsProcesses

using Meshes
using GeoTables
using Variography
using GeoStatsModels

using CpuId
using Tables
using Distances
using Distributions
using ProgressMeter

using Random
using Distributed
using LinearAlgebra

using GeoStatsBase: Ensemble, integrate
using GeoStatsBase: NearestInit, initbuff
using Bessels: gamma

include("interface.jl")
include("spde.jl")
include("seq.jl")
include("fft.jl")

export SPDEGP, SEQ, SGS

end
