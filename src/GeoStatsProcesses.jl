# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

module GeoStatsProcesses

using Meshes
using GeoTables
using GeoStatsFunctions
using GeoStatsModels

using FFTW
using CpuId
using Tables
using Distances
using Distributions
using ProgressMeter
using Bessels: gamma

using Random
using Statistics
using Distributed
using LinearAlgebra

using GeoStatsBase: Ensemble, integrate
using GeoStatsBase: NearestInit, initbuff

include("processes.jl")
include("operations.jl")

export
  # point processes
  BinomialProcess,
  ClusterProcess,
  InhibitionProcess,
  PoissonProcess,
  UnionProcess,
  ishomogeneous,

  # point operations
  RandomThinning,
  thin,

  # field processes
  GaussianProcess,
  LindgrenProcess,
  QuiltingProcess,
  TuringProcess,
  StrataProcess,

  # field methods
  LUMethod,
  FFTMethod,
  SEQMethod

end
