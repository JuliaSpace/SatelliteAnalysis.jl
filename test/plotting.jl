## Description #############################################################################
#
# Tests related to the plotting extension.
#
############################################################################################

# == File: ./src/plotting/decay_analysis.jl ================================================

# -- Function: plot_decay_analysis ---------------------------------------------------------

@testset "Function plot_decay_analysis" begin
    @test_throws(
        "The function `plot_decay_analysis` is provided by a package extension.", plot_decay_analysis(1)
    )
end

@testset "Function plot_decay_analysis [EXT]" begin
    # Hand-built decay analysis result, mimicking the output of `decay_analysis` without
    # requiring the numerical integration.
    time             = collect(range(0, 0.1; length = 20))
    date             = julian2datetime.(date_to_jd(2024, 1, 1) .+ 365.25 .* time)
    space_indices    = [
        (f107 = f, f107_avg = f - 10.0, ap = 9.0) for f in range(140.0, 180.0; length = 20)
    ]
    perigee_altitude = collect(range(300.0, 120.0; length = 20))
    apogee_altitude  = perigee_altitude .+ 10

    df = DataFrame(;
        date             = date,
        time             = time,
        space_indices    = space_indices,
        apogee_altitude  = apogee_altitude,
        perigee_altitude = perigee_altitude,
    )

    metadata!(df, "Atmospheric Model",    "NRLMSISE-00";  style = :note)
    metadata!(df, "Drag Coefficient",     2.2;            style = :note)
    metadata!(df, "Satellite Mass",       100.0;          style = :note)
    metadata!(df, "Satellite Mean Area",  1.0;            style = :note)
    metadata!(df, "Space Indices Source", "User function"; style = :note)
    metadata!(df, "SRP Coefficient",      1.25;           style = :note)
    metadata!(df, "Terminate Altitude",   120e3;          style = :note)

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

    # Custom title and subtitle, and disabled subtitle.
    fig, ax = plot_decay_analysis(df; title = "My Decay", subtitle = "Worst case")

    @test fig isa Figure
    @test ax.title[] == "My Decay"
    @test ax.subtitle[] == "Worst case"

    fig, ax = plot_decay_analysis(df; subtitle = nothing)

    @test fig isa Figure
    @test ax.subtitle[] == ""

    # The automatic subtitle must show the analysis timespan only when the dates are
    # enabled.
    fig, ax = plot_decay_analysis(df)

    @test ax.subtitle[] == ""

    fig, ax = plot_decay_analysis(df; show_dates = true)

    @test ax.subtitle[] == "2024-01-01 → 2024-02-06 UTC"

    # The axis labels must use human-readable unit names.
    @test ax.xlabel[] == "Time [years]"
    @test ax.ylabel[] == "Altitude [km]"

    # The card number formatter must group the digits of large integers.
    ext = Base.get_extension(SatelliteAnalysis, :SatelliteAnalysisMakieExt)

    @test ext._format_number(12340)   == "12 340"
    @test ext._format_number(1234567) == "1 235 000"
    @test ext._format_number(-12340)  == "-12 340"
    @test ext._format_number(1234)    == "1234"
    @test ext._format_number(3.14159) == "3.142"

    # Theme and layout keywords.
    fig, ax = plot_decay_analysis(
        df;
        fontscale       = 1.4,
        mono_ticklabels = true,
        panel_width     = 300,
        xlims           = (0.0, 0.2),
        ylims           = (100.0, 350.0)
    )

    @test fig isa Figure
    @test ax isa Axis

    # Optional decorations: mission name, absolute dates, and F10.7 twin y-axis with the
    # default getters.
    fig, ax = plot_decay_analysis(
        df;
        mission_name = "Amazonia-1",
        show_dates   = true,
        show_f107    = true
    )

    @test fig isa Figure
    @test ax isa Axis

    # Custom getters must support space indices with different field names.
    df_custom = copy(df)
    df_custom.space_indices = [
        (daily = f, mean81 = f - 10.0) for f in range(140.0, 180.0; length = 20)
    ]

    fig, ax = plot_decay_analysis(
        df_custom;
        show_f107       = true,
        f107_getter     = si -> si.daily,
        f107_avg_getter = si -> si.mean81
    )

    @test fig isa Figure
    @test ax isa Axis

    # Passing `nothing` to a getter must omit the related curve.
    fig, ax = plot_decay_analysis(df; show_f107 = true, f107_avg_getter = nothing)

    @test fig isa Figure
    @test ax isa Axis

    fig, ax = plot_decay_analysis(df; show_f107 = true, f107_getter = nothing)

    @test fig isa Figure

    # Non-finite F10.7 values must be neglected when aligning the ticks of the twin axis.
    # Notice that the function was never returning in this case.
    num_rows = size(df, 1)
    k = Ref(0)

    fig, ax = plot_decay_analysis(
        df;
        show_f107       = true,
        f107_getter     = si -> (k[] += 1) <= div(num_rows, 2) ? NaN : si.f107,
        f107_avg_getter = nothing
    )

    @test fig isa Figure

    fig, ax = plot_decay_analysis(
        df; show_f107 = true, f107_getter = si -> NaN, f107_avg_getter = nothing
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

    # The keywords with text must accept any string, and the panel width any real number.
    fig, ax = plot_decay_analysis(
        df;
        mission_name = SubString("Amazonia-1 ", 1, 10),
        panel_width  = 300.0,
        subtitle     = SubString("Subtitle ", 1, 8),
        title        = SubString("Title ", 1, 5)
    )

    @test ax.title[] == "Title"

    # The keyword `theme` also accepts a Makie theme, which is applied as it is, and
    # `nothing`, which keeps the current Makie theme.
    fig, ax = plot_decay_analysis(df; theme = Theme(; Axis = (; titlesize = 31,)))
    @test ax.titlesize[] == 31

    fig, ax = with_theme(Theme(; Axis = (; titlesize = 33,))) do
        plot_decay_analysis(df; theme = nothing)
    end
    @test ax.titlesize[] == 33

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

    # The keyword `show_f107` requires the column `space_indices`.
    @test_throws ArgumentError plot_decay_analysis(df_no_metadata; show_f107 = true)

    # At least one getter must be provided when `show_f107` is `true`.
    @test_throws ArgumentError plot_decay_analysis(
        df;
        show_f107       = true,
        f107_getter     = nothing,
        f107_avg_getter = nothing
    )

    # A getter that does not match the space indices must raise a clear error.
    @test_throws ArgumentError plot_decay_analysis(
        df;
        show_f107   = true,
        f107_getter = si -> si.not_a_field
    )

    # The keyword `subtitle` only accepts the symbol `:auto`.
    @test_throws ArgumentError plot_decay_analysis(df; subtitle = :date)
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
        "The function `plot_ground_track` is provided by a package extension.",
        plot_ground_track(1)
    )

    @test_throws(
        "The function `plot_ground_track!` is provided by a package extension.",
        plot_ground_track!(1)
    )
end

@testset "Function plot_ground_track [EXT]" begin
    using GeoJSON

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

    # == In-Place Version ==================================================================

    # The in-place version must return the plot, and it must pass the keywords to `lines!`.
    fig = Figure()
    ax  = Axis(fig[1, 1])

    plt = plot_ground_track!(ax, gt; color = :red, label = "Ground Track", linewidth = 5)

    @test plt isa Lines
    @test plt.linewidth[] == 5
    @test plt.label[] == "Ground Track"
    @test length(ax.scene.plots) == 1

    # The line width must not be fixed. Hence, it is obtained from the current theme.
    plt = with_theme(Theme(; Lines = (; linewidth = 7,))) do
        fig = Figure()
        ax  = Axis(fig[1, 1])
        plot_ground_track!(ax, gt)
    end

    @test plt.linewidth[] == 7
end

# == File: ./src/plotting/ground_facilities.jl =============================================

# -- Function: plot_ground_facility_visibility_circles -------------------------------------

@testset "Function plot_ground_facility_visibility_circles" begin
    @test_throws(
        "The function `plot_ground_facility_visibility_circles` is provided by a package",
        plot_ground_facility_visibility_circles(1)
    )

    @test_throws(
        "The function `plot_ground_facility_visibility_circles!` is provided by a package",
        plot_ground_facility_visibility_circles!(1)
    )

    # The fallback methods must also accept keywords. Otherwise, the user would see a
    # `MethodError` instead of the message describing the valid call.
    @test_throws(
        "check the arguments: the valid call is",
        plot_ground_facility_visibility_circles!(1; ground_facility_names = ["A"])
    )

    @test_throws(
        "check the arguments: the valid call is", plot_ground_track!(1; color = :red)
    )
end

@testset "Function plot_ground_facility_visibility_circles [EXT]" begin
    using GeoJSON

    gfv1 = ground_facility_visibility_circle((0, 0, 0), EARTH_EQUATORIAL_RADIUS + 700e3);
    gfv2 = ground_facility_visibility_circle(
        (-40 |> deg2rad, -60 |> deg2rad, 0), EARTH_EQUATORIAL_RADIUS + 700e3
    );

    fig, ax = plot_ground_facility_visibility_circles(
        [gfv1, gfv2]; ground_facility_names = ["GF 1", "GF 2"]
    )

    @test fig isa Figure

    fig, ax = plot_ground_facility_visibility_circles(
        [gfv1, gfv2]; ground_facility_names = split("GF1 GF2")
    )

    @test fig isa Figure
    @test ax isa Axis

    # == In-Place Version ==================================================================

    # The in-place version must return the plots of the visibility circles, and it must
    # pass the keywords to `lines!`.
    fig = Figure()
    ax  = Axis(fig[1, 1])

    plts = plot_ground_facility_visibility_circles!(
        ax, [gfv1, gfv2]; ground_facility_names = ["GF 1", "GF 2"], linewidth = 5
    )

    @test plts isa Vector{Lines}
    @test length(plts) == 2
    @test all(p -> p.linewidth[] == 5, plts)

    # Two circles, two markers, and two labels.
    @test length(ax.scene.plots) == 6

    # == Ground Facility Positions =========================================================

    # If the user provides the ground facility positions, the markers must be exactly at
    # those positions. Otherwise, they are estimated using the visibility circles.
    gfs = [(0, 0, 0), (-40 |> deg2rad, -60 |> deg2rad, 0)]

    fig = Figure()
    ax  = Axis(fig[1, 1])

    plot_ground_facility_visibility_circles!(ax, [gfv1, gfv2]; ground_facilities = gfs)

    markers = filter(p -> p isa Scatter, ax.scene.plots)

    @test length(markers) == 2
    @test markers[2][1][][1] ≈ Point2f(-60, -40)

    fig = Figure()
    ax  = Axis(fig[1, 1])

    plot_ground_facility_visibility_circles!(ax, [gfv1, gfv2])

    markers = filter(p -> p isa Scatter, ax.scene.plots)

    @test markers[2][1][][1] ≈ Point2f(-60, -40) atol = 0.5
    @test markers[2][1][][1] != Point2f(-60, -40)

    fig, ax = plot_ground_facility_visibility_circles([gfv1, gfv2]; ground_facilities = gfs)

    @test fig isa Figure

    # == Errors ============================================================================

    @test_throws ArgumentError plot_ground_facility_visibility_circles(
        [gfv1, gfv2]; ground_facility_names = ["GF 1", "GF 2", "GF 3"]
    )

    @test_throws ArgumentError plot_ground_facility_visibility_circles(
        [gfv1, gfv2]; ground_facilities = [(0, 0, 0)]
    )

    @test_throws ArgumentError plot_ground_facility_visibility_circles!(
        ax, [gfv1, gfv2]; ground_facility_names = ["GF 1"]
    )
end

# == File: ./src/plotting/world_map.jl =====================================================

@testset "Function plot_world_map" begin
    @test_throws(
        "The function `plot_world_map` is provided by a package extension.",
        plot_world_map(1)
    )

    @test_throws(
        "The function `plot_world_map!` is provided by a package extension.",
        plot_world_map!(1; theme = :dark)
    )
end

@testset "Function plot_world_map [EXT]" begin
    using GeoJSON

    fig, ax = plot_world_map()

    @test fig isa Figure
    @test ax isa Axis

    fig, ax = plot_world_map(; theme = :dark)

    @test fig isa Figure
    @test ax isa Axis

    # == Errors ============================================================================

    @test_throws ArgumentError plot_world_map(; theme = :blue)

    # == In-Place Version ==================================================================

    # The in-place version must draw only the country polygons, returning the plot.
    fig = Figure()
    ax  = Axis(fig[1, 1]; title = "My Map")

    plt = plot_world_map!(ax; theme = :dark, strokewidth = 3)

    @test plt isa Poly
    @test plt.strokewidth[] == 3
    @test ax.title[] == "My Map"
    @test length(ax.scene.plots) == 1

    @test_throws ArgumentError plot_world_map!(ax; theme = :blue)

    # The keyword `theme` also accepts a Makie theme, which is applied as it is, and
    # `nothing`, which keeps the current Makie theme.
    custom_theme = Theme(; Axis = (; titlesize = 31,))

    fig, ax = plot_world_map(; theme = custom_theme)
    @test ax.titlesize[] == 31

    fig, ax = with_theme(Theme(; Axis = (; titlesize = 33,))) do
        plot_world_map(; theme = nothing)
    end
    @test ax.titlesize[] == 33

    # The default theme must not use the current Makie theme.
    fig, ax = with_theme(Theme(; Axis = (; titlesize = 33,))) do
        plot_world_map()
    end
    @test ax.titlesize[] != 33
end
