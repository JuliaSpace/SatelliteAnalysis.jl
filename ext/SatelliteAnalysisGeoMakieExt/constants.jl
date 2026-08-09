## Description #############################################################################
#
# Definition of constants.
#
############################################################################################

# Colors used to draw the country polygons in the world map, keyed by the theme variant.
# The dark values mirror `SEPARATOR_DARK` and `BORDER_DARK` from the extension
# SatelliteAnalysisMakieExt (see its `constants.jl`). They are duplicated here as plain hex
# strings because one extension cannot access another extension's constants at precompile
# time.
const _COUNTRY_FILL_COLOR   = Dict(:light => "#FFFFFF", :dark => "#162940")
const _COUNTRY_STROKE_COLOR = Dict(:light => "#000000", :dark => "#1E3A5F")
