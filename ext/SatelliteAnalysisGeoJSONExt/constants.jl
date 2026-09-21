## Description #############################################################################
#
# Definition of constants.
#
############################################################################################

# Colors used to draw the country polygons in the world map, keyed by the theme variant.
# They are colors of the SatelliteAnalysis.jl Makie theme, which are defined in the main
# package because one extension cannot access the constants of another extension.
const _COUNTRY_FILL_COLOR = Dict(
    :light => SatelliteAnalysis._THEME_SURFACE_HEX,
    :dark  => SatelliteAnalysis._THEME_SEPARATOR_DARK_HEX,
)

const _COUNTRY_STROKE_COLOR = Dict(
    :light => SatelliteAnalysis._THEME_BORDER_LIGHT_HEX,
    :dark  => SatelliteAnalysis._THEME_BORDER_DARK_HEX,
)
