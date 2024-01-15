## Description #############################################################################
#
# Tests related to the plotting extension.
#
############################################################################################

# == File: ./src/plotting/fetch_country_polygons.jl ========================================

# -- Function: fetch_country_polygons ------------------------------------------------------

@testset "Function fetch_country_polygons" begin
    f1 = @test_logs(
        (:info, "Downloading the file 'countries.geojson' from 'https://datahub.io/core/geo-countries/r/countries.geojson'..."),
        fetch_country_polygons()
    )

    f2 = @test_logs(
        min_level = Logging.Warn,
        fetch_country_polygons()
    )

    @test f1 == f2
end

# == File: ./src/plotting/ground_track.jl ==================================================

# -- Function: plot_ground_track -----------------------------------------------------------

@testset "Function plot_ground_track" begin
    @test_logs(
        (:error, "The package GeoMakie.jl must be loaded to use this functionality."),
        plot_ground_track(1)
    )
end
