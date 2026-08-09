## Description #############################################################################
#
# Tests related to the Makie themes provided by the extension SatelliteAnalysisMakieExt.
#
############################################################################################

# == File: ./ext/SatelliteAnalysisMakieExt/theme.jl ========================================

SAMakieExt = Base.get_extension(SatelliteAnalysis, :SatelliteAnalysisMakieExt)

# -- Function makie_palette ----------------------------------------------------------------

@testset "Function makie_palette" begin
    @test SatelliteAnalysis.makie_palette(6) == SAMakieExt.CATEGORICAL_LIGHT
    @test SatelliteAnalysis.makie_palette(6; dark = true) == SAMakieExt.CATEGORICAL_DARK
    @test SatelliteAnalysis.makie_palette(3) == SAMakieExt.CATEGORICAL_LIGHT[1:3]
    @test isempty(SatelliteAnalysis.makie_palette(0))
end

@testset "Function makie_palette [ERRORS]" begin
    @test_throws ArgumentError SatelliteAnalysis.makie_palette(7)
    @test_throws ArgumentError SatelliteAnalysis.makie_palette(-1)
end

# -- Function makie_theme ------------------------------------------------------------------

@testset "Function makie_theme" begin
    theme = SatelliteAnalysis.makie_theme()
    @test theme isa Makie.Theme
    @test to_value(theme.backgroundcolor) == SAMakieExt.SURFACE

    dark_theme = SatelliteAnalysis.makie_theme(:dark)
    @test dark_theme isa Makie.Theme
    @test to_value(dark_theme.backgroundcolor) == SAMakieExt.NAVY_PRIMARY

    # The Symbol and Val call forms must build the same theme.
    @test to_value(SatelliteAnalysis.makie_theme(Val(:dark)).backgroundcolor) ==
        to_value(dark_theme.backgroundcolor)
    @test to_value(SatelliteAnalysis.makie_theme(Val(:light)).backgroundcolor) ==
        to_value(theme.backgroundcolor)

    # The keyword `fontscale` must scale the font sizes.
    scaled_theme = SatelliteAnalysis.makie_theme(; fontscale = 2)
    @test to_value(scaled_theme.Axis.xticklabelsize) ==
        2 * to_value(theme.Axis.xticklabelsize)

    # The keyword `mono_ticklabels` must switch the tick label font to the bundled IBM Plex
    # Mono.
    mono_theme = SatelliteAnalysis.makie_theme(; mono_ticklabels = true)
    @test to_value(theme.Axis.xticklabelfont) == :regular
    @test endswith(to_value(mono_theme.Axis.xticklabelfont), "IBMPlexMono-Regular.ttf")
end

@testset "Function makie_theme [ERRORS]" begin
    @test_throws ArgumentError SatelliteAnalysis.makie_theme(:blue)
    @test_throws ArgumentError SatelliteAnalysis.makie_theme(Val(:blue))
end

# NOTE: The fallback errors thrown when Makie.jl is not loaded cannot be tested here because
# Makie.jl is loaded for this test set and package extensions cannot be unloaded.
