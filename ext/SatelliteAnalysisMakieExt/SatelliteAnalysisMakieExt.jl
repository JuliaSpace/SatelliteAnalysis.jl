module SatelliteAnalysisMakieExt

using Dates
using Makie
using Makie: Colorant
using SatelliteAnalysis
import SatelliteAnalysis: makie_theme, makie_palette, plot_decay_analysis

############################################################################################
#                                        Constants                                         #
############################################################################################

include("./constants.jl")

############################################################################################
#                                         Includes                                         #
############################################################################################

include("./decay_analysis.jl")
include("./theme.jl")

end # module SatelliteAnalysisMakieExt
