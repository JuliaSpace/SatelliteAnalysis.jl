module SatelliteAnalysis

using LinearAlgebra
using ReferenceFrameRotations
using StaticArrays
using SatelliteToolbox

################################################################################
#                                   Includes
################################################################################

include("./beta_angle.jl")
include("./lighting_condition.jl")
include("./misc/find_crossing.jl")

end # module
