module SatelliteAnalysis

using Dates
using Downloads
using LinearAlgebra
using Logging
using PrettyNumbers
using Reexport
using Scratch
using StaticArrays
using Statistics

@reexport using DataFrames
@reexport using ReferenceFrameRotations
@reexport using SatelliteToolbox

############################################################################################
#                                         Includes                                         #
############################################################################################

include("./beta_angle.jl")
include("./eclipse_time.jl")
include("./frozen_orbits.jl")
include("./ground_repeating_orbits.jl")
include("./ground_track.jl")
include("./lighting_condition.jl")
include("./sun_synchronous_orbits.jl")

include("./misc/find_crossing.jl")

include("./ground_facilities/is_ground_facility_visible.jl")
include("./ground_facilities/ground_facility_accesses_and_gaps.jl")
include("./ground_facilities/ground_facility_visibility_circle.jl")

include("./plotting/fetch_country_polygons.jl")
include("./plotting/ground_facilities.jl")
include("./plotting/ground_track.jl")

end # module
