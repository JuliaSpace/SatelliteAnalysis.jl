module SatelliteAnalysisMakieExt

using Makie
using Makie: Colorant
using SatelliteAnalysis
import SatelliteAnalysis: makie_theme, makie_palette

############################################################################################
#                                        Constants                                         #
############################################################################################

include("./constants.jl")

############################################################################################
#                                         Includes                                         #
############################################################################################

include("./theme.jl")

end # module SatelliteAnalysisMakieExt
