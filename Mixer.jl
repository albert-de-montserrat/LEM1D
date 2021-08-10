using BenchmarkTools
using LinearAlgebra
using SparseArrays
using Interpolations
using DelimitedFiles
using NumericalIntegration
using LoopVectorization
using StaticArrays
using TimerOutputs
using HDF5
using GLMakie
# using CSV,DataFrames
include("src/ProcessParameters.jl")
include("src/read_txt_files.jl")
include("src/solvers.jl")
include("src/Knickpoints.jl")
include("solver_dispatch.jl")
include("src/others.jl")
include("src/FEM.jl")
include("src/river_functions.jl")

# TIME CONSTANTS ==========================================================
const yr = 60.0 * 60.0 * 24.0 * 365.0;