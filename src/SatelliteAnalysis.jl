module SatelliteAnalysis

using Dates
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

include("./ground_stations/is_ground_station_visible.jl")
include("./ground_stations/ground_station_accesses_and_gaps.jl")

end # module
