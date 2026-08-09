module SatelliteAnalysisGeoMakieExt

using SatelliteAnalysis
using GeoMakie

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

end # module SatelliteAnalysisGeoMakieExt
