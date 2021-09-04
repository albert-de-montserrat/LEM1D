module LEM1D

using LinearAlgebra
using SparseArrays
using Interpolations
using DelimitedFiles
using NumericalIntegration
using LoopVectorization
using StaticArrays
using TimerOutputs
using HDF5
using Polyester

# using CSV,DataFrames
include("ProcessParameters.jl")
export FluxLimiter
export KnickPointsParameters
export KnickPointsArrays
export TerraceParameters
export HillSlopeParameters
export SubmarineParameters
export LimiterFunction
export minmod
export superbee
export vanAlbada
export vanAlbada2
export vanLeer
export ospre
export hquick
export koren

include("read_txt_files.jl")
export SeaLevel
export Bintanja
export Spratt450
export Spratt800
export SeaLevelOpts
export sea_level_starting_time
export get_sea_level
export sea_level_corrections
export read_uplift

include("Knickpoints.jl")
export fluxlimiter
export updater!
export TVD!

include("solvers.jl")
export Model 
export TerraceOnly
export TerraceKnickpoint
export TerraceDiffusion
export TerraceKnickpointDiffusion
export solve
export ModelOpts

include("others.jl")
export uplift!

include("FEM.jl")
export RiverProfile
export TerraceProfile
export Profile
export mesher 
export sparsitymatrix
export ScratchTerraceFEM

include("river_functions.jl")
export find_shore_id

# TIME CONSTANTS ==========================================================
const yr = 60.0 * 60.0 * 24.0 * 365.0
const km = 1e3
export km
export yr

end # module
