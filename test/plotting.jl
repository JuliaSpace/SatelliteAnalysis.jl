## Description #############################################################################
#
# Tests related to the plotting extension.
#
############################################################################################

# == File: ./src/plotting/fetch_country_polygons.jl ========================================

# -- Function: fetch_country_polygons ------------------------------------------------------

@testset "Function fetch_country_polygons" begin
    f1 = @test_logs(
        (:info, "Downloading the file 'countries.geojson' from 'https://pkgstore.datahub.io/core/geo-countries/countries/archive/23f420f929e0e09c39d916b8aaa166fb/countries.geojson'..."),
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

@testset "Function plot_ground_track [EXT]" begin
    using GeoMakie

    jd₀ = date_to_jd(2021, 1, 1)
    2.4592155e6

    orb = KeplerianElements(
        jd₀,
        7130.982e3,
        0.001111,
        98.405 |> deg2rad,
        ltdn_to_raan(10.5, jd₀),
        π / 2,
        0
    )

    orbp = Propagators.init(Val(:J2), orb)

    gt = ground_track(orbp; track_types = :descending, duration = 5 * 86400);

    fig, ax = plot_ground_track(gt; size = (2000, 1000))

    @test fig isa Figure
    @test ax isa Axis
end

# == File: ./src/plotting/ground_facilities.jl =============================================

# -- Function: plot_ground_facility_visibility_circles -------------------------------------

@testset "Function plot_ground_facility_visibility_circles" begin
    @test_logs(
        (:error, "The package GeoMakie.jl must be loaded to use this functionality."),
        plot_ground_facility_visibility_circles(1)
    )
end

@testset "Function plot_ground_facility_visibility_circles [EXT]" begin
    using GeoMakie

    gfv1 = ground_facility_visibility_circle((0, 0, 0), EARTH_EQUATORIAL_RADIUS + 700e3);
    gfv2 = ground_facility_visibility_circle((-40 |> deg2rad, -60 |> deg2rad, 0), EARTH_EQUATORIAL_RADIUS + 700e3);

    fig, ax = plot_ground_facility_visibility_circles(
        [gfv1, gfv2];
        ground_facility_names = ["GF 1", "GF 2"]
    )

    @test fig isa Figure
    @test ax isa Axis

    # == Errors ============================================================================

    @test_throws ArgumentError plot_ground_facility_visibility_circles(
        [gfv1, gfv2];
        ground_facility_names = ["GF 1", "GF 2", "GF 3"]
    )
end

