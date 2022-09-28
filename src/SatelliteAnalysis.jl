module SatelliteAnalysis

using Dates
using DataFrames
using LinearAlgebra
using PrettyTables
using ReferenceFrameRotations
using StaticArrays
using SatelliteToolbox
using Statistics
using TerminalPager

################################################################################
#                                   Includes
################################################################################

include("./beta_angle.jl")
include("./eclipse_time.jl")
include("./lighting_condition.jl")
include("./misc/find_crossing.jl")

include("./ground_facilities/is_ground_facility_visible.jl")
include("./ground_facilities/ground_facility_accesses_and_gaps.jl")

end # module
