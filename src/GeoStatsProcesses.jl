# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

module GeoStatsProcesses

using Meshes
using GeoTables
using Variography
using GeoStatsModels

using FFTW
using CpuId
using Tables
using Distances
using Distributions
using ProgressMeter

using Random
using Statistics
using Distributed
using LinearAlgebra

using GeoStatsBase: Ensemble, integrate
using GeoStatsBase: NearestInit, initbuff
using GeoStatsModels: fitpredict
using Bessels: gamma

include("interface.jl")
include("spde.jl")
include("seq.jl")
include("fft.jl")

export SPDEGP, SEQ, SGP, FFTGP

end
