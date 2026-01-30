# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

module GeoStatsProcesses

using Meshes
using Unitful
using GeoTables
using GeoStatsBase
using GeoStatsFunctions
using GeoStatsModels

using FFTW
using Tables
using Distances
using Distributions
using TableTransforms
using ProgressMeter
using Bessels: gamma
using CpuId: cpucores

using Random
using Statistics
using Distributed
using LinearAlgebra

import Distributions: mean, var
import Distributions: cdf, ccdf
import Distributions: quantile
import Base: ==

include("utils.jl")
include("ensembles.jl")
include("processes.jl")
include("operations.jl")
include("initialization.jl")
include("simulation.jl")
include("expectation.jl")

export
  # ensembles
  Ensemble,
  mean,
  var,
  cdf,
  ccdf,
  quantile,

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
  IndicatorProcess,
  LindgrenProcess,
  QuiltingProcess,
  TuringProcess,
  StrataProcess,
  iscontinuous,

  # field simulation
  LUSIM,
  SEQSIM,
  FFTSIM,

  # field expectation
  expectedvalue,
  mean,

  # initialization
  NearestInit,
  ExplicitInit

end
