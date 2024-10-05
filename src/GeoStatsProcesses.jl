# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

module GeoStatsProcesses

using Meshes
using Unitful
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

import Distributions: mean, var, cdf, quantile
import Base: ==

include("units.jl")
include("ioutils.jl")
include("ensembles.jl")
include("initbuff.jl")
include("processes.jl")
include("operations.jl")

export
  # ensembles
  Ensemble,
  mean,
  var,
  cdf,
  quantile,

  # initialization
  InitMethod,
  NearestInit,
  ExplicitInit,
  initbuff,

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
