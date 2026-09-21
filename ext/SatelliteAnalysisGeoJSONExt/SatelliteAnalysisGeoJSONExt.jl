module SatelliteAnalysisGeoJSONExt

using GeoJSON
using Makie
using SatelliteAnalysis

############################################################################################
#                                        Constants                                         #
############################################################################################

include("./constants.jl")

############################################################################################
#                                         Includes                                         #
############################################################################################

include("./ground_facilities.jl")
include("./ground_track.jl")
include("./world_map.jl")

end # module SatelliteAnalysisGeoJSONExt
