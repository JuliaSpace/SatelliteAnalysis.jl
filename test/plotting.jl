## Description #############################################################################
#
# Tests related to the plotting extension.
#
############################################################################################

# == File: ./src/plotting/decay_analysis.jl ================================================

# -- Function: plot_decay_analysis ---------------------------------------------------------

@testset "Function plot_decay_analysis" begin
    @test_throws(
        "Wrong input or the package Makie.jl is not loaded.", plot_decay_analysis(1)
    )
end

@testset "Function plot_decay_analysis [EXT]" begin
    # Hand-built decay analysis result, mimicking the output of `decay_analysis` without
    # requiring the numerical integration.
    time             = collect(range(0, 0.1; length = 20))
    date             = julian2datetime.(date_to_jd(2024, 1, 1) .+ 365.25 .* time)
    f107             = collect(range(140.0, 180.0; length = 20))
    perigee_altitude = collect(range(300.0, 120.0; length = 20))
    apogee_altitude  = perigee_altitude .+ 10

    df = DataFrame(;
        date             = date,
        time             = time,
        f107             = f107,
        apogee_altitude  = apogee_altitude,
        perigee_altitude = perigee_altitude,
    )

    metadata!(df, "Atmospheric Model",   "NRLMSISE-00";        style = :note)
    metadata!(df, "Drag Coefficient",    2.2;                  style = :note)
    metadata!(df, "F10.7 Source",        "Constant (140 sfu)"; style = :note)
    metadata!(df, "Satellite Mass",      100.0;                style = :note)
    metadata!(df, "Satellite Mean Area", 1.0;                  style = :note)
    metadata!(df, "SRP Coefficient",     1.25;                 style = :note)
    metadata!(df, "Terminate Altitude",  120e3;                style = :note)

    colmetadata!(df, :time,             "Unit", :y;  style = :note)
    colmetadata!(df, :apogee_altitude,  "Unit", :km; style = :note)
    colmetadata!(df, :perigee_altitude, "Unit", :km; style = :note)

    fig, ax = plot_decay_analysis(df)

    @test fig isa Figure
    @test ax isa Axis

    fig, ax = plot_decay_analysis(df; theme = :dark)

    @test fig isa Figure
    @test ax isa Axis

    # The ballistic coefficient and assumptions cards can be turned off.
    fig, ax = plot_decay_analysis(df; show_assumptions = false)

    @test fig isa Figure
    @test ax isa Axis

    # The in-plot reentry callout can be turned off.
    fig, ax = plot_decay_analysis(df; show_reentry_callout = false)

    @test fig isa Figure
    @test ax isa Axis

    # Optional decorations: mission name, reentry date, and F10.7 twin y-axis.
    fig, ax = plot_decay_analysis(
        df;
        mission_name      = "Amazonia-1",
        show_f107         = true,
        show_reentry_date = true
    )

    @test fig isa Figure
    @test ax isa Axis

    # A `DataFrame` without the metadata must still be plottable. In this case, the
    # information panel is omitted.
    df_no_metadata = DataFrame(;
        date             = date,
        time             = time,
        apogee_altitude  = apogee_altitude,
        perigee_altitude = perigee_altitude,
    )

    fig, ax = plot_decay_analysis(df_no_metadata)

    @test fig isa Figure
    @test ax isa Axis

    # The keywords must override the missing metadata.
    fig, ax = plot_decay_analysis(
        df_no_metadata;
        satellite_mass      = 42.0,
        satellite_mean_area = 0.5,
        terminate_altitude  = 120e3
    )

    @test fig isa Figure
    @test ax isa Axis

    # Partial assumption metadata must render only the resolvable lines.
    df_partial = copy(df_no_metadata)
    metadata!(df_partial, "Drag Coefficient", 2.0; style = :note)

    fig, ax = plot_decay_analysis(df_partial)

    @test fig isa Figure
    @test ax isa Axis

    # Analysis without a reentry.
    df_no_reentry = DataFrame(;
        date             = date,
        time             = time,
        apogee_altitude  = apogee_altitude,
        perigee_altitude = collect(range(300.0, 250.0; length = 20)),
    )

    fig, ax = plot_decay_analysis(df_no_reentry; terminate_altitude = 120e3)

    @test fig isa Figure
    @test ax isa Axis

    # == Errors ============================================================================

    @test_throws ArgumentError plot_decay_analysis(df; theme = :blue)
    @test_throws ArgumentError plot_decay_analysis(DataFrame(; a = [1]))
    @test_throws ArgumentError plot_decay_analysis(empty!(copy(df)))

    # The keyword `show_f107` requires the column `f107`.
    @test_throws ArgumentError plot_decay_analysis(df_no_metadata; show_f107 = true)
end

# == File: ./src/plotting/fetch_country_polygons.jl ========================================

# -- Function: fetch_country_polygons ------------------------------------------------------

@testset "Function fetch_country_polygons" begin
    SatelliteAnalysis.Scratch.clear_scratchspaces!(SatelliteAnalysis)

    f1 = @test_logs(
        (
            :info,
            "Downloading the file 'countries.geojson' from 'https://pkgstore.datahub.io/core/geo-countries/countries/archive/23f420f929e0e09c39d916b8aaa166fb/countries.geojson'...",
        ),
        fetch_country_polygons()
    )

    f2 = @test_logs(min_level = Logging.Warn, fetch_country_polygons())

    @test f1 == f2
end

# == File: ./src/plotting/ground_track.jl ==================================================

# -- Function: plot_ground_track -----------------------------------------------------------

@testset "Function plot_ground_track" begin
    @test_throws(
        "Wrong input or the package GeoMakie.jl is not loaded.", plot_ground_track(1)
    )

    @test_throws(
        "Wrong input or the package GeoMakie.jl is not loaded.", plot_ground_track!(1)
    )
end

@testset "Function plot_ground_track [EXT]" begin
    using GeoMakie

    jd₀ = date_to_jd(2021, 1, 1)
    2.4592155e6

    orb = KeplerianElements(
        jd₀, 7130.982e3, 0.001111, 98.405 |> deg2rad, ltdn_to_raan(10.5, jd₀), π / 2, 0
    )

    orbp = Propagators.init(Val(:J2), orb)

    gt = ground_track(orbp; track_types = :descending, duration = 5 * 86400);

    fig, ax = plot_ground_track(gt; size = (2000, 1000))

    @test fig isa Figure
    @test ax isa Axis

    fig, ax = plot_ground_track(gt; theme = :dark)

    @test fig isa Figure
    @test ax isa Axis
end

# == File: ./src/plotting/ground_facilities.jl =============================================

# -- Function: plot_ground_facility_visibility_circles -------------------------------------

@testset "Function plot_ground_facility_visibility_circles" begin
    @test_throws(
        "Wrong input or the package GeoMakie.jl is not loaded.",
        plot_ground_facility_visibility_circles(1)
    )

    @test_throws(
        "Wrong input or the package GeoMakie.jl is not loaded.",
        plot_ground_facility_visibility_circles!(1)
    )
end

@testset "Function plot_ground_facility_visibility_circles [EXT]" begin
    using GeoMakie

    gfv1 = ground_facility_visibility_circle((0, 0, 0), EARTH_EQUATORIAL_RADIUS + 700e3);
    gfv2 = ground_facility_visibility_circle(
        (-40 |> deg2rad, -60 |> deg2rad, 0), EARTH_EQUATORIAL_RADIUS + 700e3
    );

    fig, ax = plot_ground_facility_visibility_circles(
        [gfv1, gfv2]; ground_facility_names = ["GF 1", "GF 2"]
    )

    @test fig isa Figure
    @test ax isa Axis

    # == Errors ============================================================================

    @test_throws ArgumentError plot_ground_facility_visibility_circles(
        [gfv1, gfv2]; ground_facility_names = ["GF 1", "GF 2", "GF 3"]
    )
end

# == File: ./src/plotting/world_map.jl =====================================================

@testset "Function plot_world_map" begin
    @test_throws("Wrong input or the package GeoMakie.jl is not loaded.", plot_world_map(1))
end

@testset "Function plot_world_map [EXT]" begin
    using GeoMakie

    fig, ax = plot_world_map()

    @test fig isa Figure
    @test ax isa Axis

    fig, ax = plot_world_map(; theme = :dark)

    @test fig isa Figure
    @test ax isa Axis

    # == Errors ============================================================================

    @test_throws ArgumentError plot_world_map(; theme = :blue)
end
