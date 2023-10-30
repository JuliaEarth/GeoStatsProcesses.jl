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

include("processes.jl")
include("operations.jl")

export
  # field processes
  SPDEGP, 
  SEQ, 
  GaussianProcess, 
  IQP, 
  TP, 
  SP,

  # random methods
  LUMethod,
  FFTMethod,
  SEQMethod,

  # point processes
  BinomialProcess, 
  ClusterProcess, 
  InhibitionProcess, 
  PoissonProcess, 
  UnionProcess,

  # operations
  RandomThinning, 
  thin

end
