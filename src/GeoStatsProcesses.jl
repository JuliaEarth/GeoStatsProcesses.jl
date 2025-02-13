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

include("utils.jl")
include("initbuff.jl")
include("ensembles.jl")
include("processes.jl")
include("operations.jl")
include("simulation.jl")

export
  # initialization
  NearestInit,
  ExplicitInit,

  # ensembles
  Ensemble,
  mean,
  var,
  cdf,
  quantile,

  # point processes
  BinomialProcess,
  ClusterProcess,
  InhibitionProcess,
  PoissonProcess,
  UnionProcess,
  ishomogeneous,

  # field processes
  GaussianProcess,
  LindgrenProcess,
  QuiltingProcess,
  TuringProcess,
  StrataProcess,

  # point operations
  RandomThinning,
  thin,

  # field simulation
  LUSIM,
  SEQSIM,
  FFTSIM

end
