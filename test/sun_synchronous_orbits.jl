## Description #############################################################################
#
# Tests related to the Sun-synchronous orbit functions.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# == File: ./src/sun_synchronous_orbits.jl =================================================

# -- Function: design_sun_sync_ground_repeating_orbits -------------------------------------

@testset "Designing Sun-Synchronous, Ground-Repeating Orbits" begin

    # == Default ===========================================================================

    df = design_sun_sync_ground_repeating_orbit(1, 1)

    @test size(df) == (5, 7)

    @test df[begin, :semi_major_axis] ≈ 6382.409 (atol = 1e-3)
    @test df[begin, :altitude] ≈ 4.272 (atol = 1e-3)
    @test df[begin, :inclination] ≈ 95.6949 (atol = 1e-4)
    @test df[begin, :period] ≈ 84.706 (atol = 1e-3)
    @test df[begin, :revs_per_day] == "17"
    @test df[begin, :adjacent_gt_distance] ≈ 2327.845 (atol = 1e-3)
    @test df[begin, :adjacent_gt_angle] ≈ 169.1250 (atol = 1e-4)

    @test df[end, :semi_major_axis] ≈ 7635.252 (atol = 1e-3)
    @test df[end, :altitude] ≈ 1257.115 (atol = 1e-3)
    @test df[end, :inclination] ≈ 100.7057 (atol = 1e-4)
    @test df[end, :period] ≈ 110.769 (atol = 1e-3)
    @test df[end, :revs_per_day] == "13"
    @test df[end, :adjacent_gt_distance] ≈ 2976.209 (atol = 1e-3)
    @test df[end, :adjacent_gt_angle] ≈ 91.7639 (atol = 1e-4)

    # == Metadata ==========================================================================

    @test metadata(df, "Description") ==
        "Sun-synchronous, ground-repeating orbits with repetition between 1 and 1 days."

    @test colmetadata(df, :semi_major_axis, "Unit") == :km
    @test colmetadata(df, :altitude, "Unit") == :km
    @test colmetadata(df, :inclination, "Unit") == :deg
    @test colmetadata(df, :period, "Unit") == :min
    @test colmetadata(df, :adjacent_gt_distance, "Unit") == :km
    @test colmetadata(df, :adjacent_gt_angle, "Unit") == :deg

    # All the metadata must use the style `:note` to propagate through DataFrame
    # transformations.
    @test all(k -> metadata(df, k; style = true)[2] == :note, metadatakeys(df))
    @test all(
        colmetadata(df, col, k; style = true)[2] == :note for
        (col, keys) in colmetadatakeys(df) for k in keys
    )

    df_units = design_sun_sync_ground_repeating_orbit(
        1, 1; angle_unit = :rad, distance_unit = :m, time_unit = :s
    )

    @test colmetadata(df_units, :altitude, "Unit") == :m
    @test colmetadata(df_units, :adjacent_gt_angle, "Unit") == :rad
    @test colmetadata(df_units, :period, "Unit") == :s

    # == Altitude filter ===================================================================

    df = design_sun_sync_ground_repeating_orbit(
        5, 5; minimum_altitude = 750e3, maximum_altitude = 760e3
    )

    @test size(df) == (1, 7)
    @test df[begin, :semi_major_axis] ≈ 7130.984 (atol = 1e-3)
    @test df[begin, :altitude] ≈ 752.847 (atol = 1e-3)
    @test df[begin, :inclination] ≈ 98.4106 (atol = 1e-4)
    @test df[begin, :period] ≈ 100.000 (atol = 1e-3)
    @test df[begin, :revs_per_day] == "14 + ²/₅"
    @test df[begin, :adjacent_gt_distance] ≈ 543.811 (atol = 1e-3)
    @test df[begin, :adjacent_gt_angle] ≈ 39.425 (atol = 1e-3)

    corrected = design_sun_sync_ground_repeating_orbit(
        5,
        5;
        minimum_altitude = 750e3,
        maximum_altitude = 760e3,
        pretty_revs_per_day = false,
    )
    @test corrected[begin, :revs_per_day] == (14, 2 // 5)

    # == Revolutions per day ===============================================================

    df = design_sun_sync_ground_repeating_orbit(
        5,
        5;
        minimum_altitude = 750e3,
        maximum_altitude = 760e3,
        pretty_revs_per_day = false,
    )

    @test size(df) == (1, 7)
    @test df[begin, :semi_major_axis] ≈ 7130.984 (atol = 1e-3)
    @test df[begin, :altitude] ≈ 752.847 (atol = 1e-3)
    @test df[begin, :inclination] ≈ 98.4106 (atol = 1e-4)
    @test df[begin, :period] ≈ 100.000 (atol = 1e-3)
    @test df[begin, :revs_per_day] == (14, 2 // 5)
    @test df[begin, :adjacent_gt_distance] ≈ 543.811 (atol = 1e-3)
    @test df[begin, :adjacent_gt_angle] ≈ 39.425 (atol = 1e-3)

    # == Test units ========================================================================

    # -- Angle -----------------------------------------------------------------------------

    df = design_sun_sync_ground_repeating_orbit(
        5, 5; minimum_altitude = 750e3, maximum_altitude = 760e3, angle_unit = :rad
    )

    @test df[begin, :inclination] ≈ 98.4106 |> deg2rad (atol = 2e-6)
    @test df[begin, :adjacent_gt_angle] ≈ 39.425 |> deg2rad (atol = 8e-5)

    # -- Distance --------------------------------------------------------------------------

    df = design_sun_sync_ground_repeating_orbit(
        5, 5; minimum_altitude = 750e3, maximum_altitude = 760e3, distance_unit = :m
    )

    @test df[begin, :semi_major_axis] ≈ 7130.984e3 (atol = 1)
    @test df[begin, :altitude] ≈ 752.847e3 (atol = 1)
    @test df[begin, :adjacent_gt_distance] ≈ 543.811e3 (atol = 1)

    # -- Time ------------------------------------------------------------------------------

    df = design_sun_sync_ground_repeating_orbit(
        5, 5; minimum_altitude = 750e3, maximum_altitude = 760e3, time_unit = :s
    )

    @test df[begin, :period] ≈ 6000.000 (atol = 1e-3)

    df = design_sun_sync_ground_repeating_orbit(
        5, 5; minimum_altitude = 750e3, maximum_altitude = 760e3, time_unit = :h
    )

    @test df[begin, :period] ≈ 1.66667 (atol = 1e-5)

    # == Errors ============================================================================

    @test_throws ArgumentError design_sun_sync_ground_repeating_orbit(0, 10)
    @test_throws ArgumentError design_sun_sync_ground_repeating_orbit(3, -3)
    @test_throws ArgumentError design_sun_sync_ground_repeating_orbit(3, 2)
    @test_throws ArgumentError design_sun_sync_ground_repeating_orbit(
        1, 5; eccentricity = 1
    )

    @test_throws ArgumentError design_sun_sync_ground_repeating_orbit(
        1, 5; angle_unit = :unknown
    )

    @test_throws ArgumentError design_sun_sync_ground_repeating_orbit(
        1, 5; distance_unit = :unknown
    )

    @test_throws ArgumentError design_sun_sync_ground_repeating_orbit(
        1, 5; time_unit = :unknown
    )

    # == Revolutions per Day Without a Sun-Synchronous Orbit ===============================

    # If there is no Sun-synchronous orbit for a number of revolutions per day, it must be
    # skipped without throwing exceptions or printing warnings.
    df = @test_logs design_sun_sync_ground_repeating_orbit(1, 1; int_rev_per_day = (5, 14))

    @test size(df) == (1, 7)
    @test df[begin, :revs_per_day] == "14"

    df = @test_logs design_sun_sync_ground_repeating_orbit(1, 1; int_rev_per_day = (5,))

    @test size(df) == (0, 7)
end

# -- Function: sun_sync_orbit_from_angular_velocity ----------------------------------------

@testset "Function sun_sync_orbit_from_angular_velocity" begin
    # == Float64 ===========================================================================

    a, i, c = sun_sync_orbit_from_angular_velocity(0.06 |> deg2rad)
    @test a ≈ 7.130983932846816e6
    @test i ≈ 1.7175898375139984
    @test c == true
    @test eltype(a) == Float64
    @test eltype(i) == Float64

    a, i, c = sun_sync_orbit_from_angular_velocity(0.06 |> deg2rad, 0.05)
    @test a ≈ 7.130953853502454e6
    @test i ≈ 1.7168497436578767
    @test c == true
    @test eltype(a) == Float64
    @test eltype(i) == Float64

    # == Float32 ===========================================================================

    a, i, c = sun_sync_orbit_from_angular_velocity(0.06f0 |> deg2rad)
    @test a ≈ 7.130984f6
    @test i ≈ 1.7175899f0
    @test c == true
    @test eltype(a) == Float32
    @test eltype(i) == Float32

    a, i, c = sun_sync_orbit_from_angular_velocity(0.06f0 |> deg2rad, 0.05f0)
    @test a ≈ 7.1309545f6
    @test i ≈ 1.7168498f0
    @test c == true
    @test eltype(a) == Float32
    @test eltype(i) == Float32

    # == Test When Algorithm Did Not Converge ==============================================

    a, i, c = (@test_logs(
        (:warn,), sun_sync_orbit_from_angular_velocity(0.06 |> deg2rad; max_iterations = 3)
    ))
    @test c == false

    # A loose tolerance must lead to a less accurate result without warnings.
    a_ref, i_ref, ~ = sun_sync_orbit_from_angular_velocity(0.06 |> deg2rad)
    a, i, c = @test_logs sun_sync_orbit_from_angular_velocity(
        0.06 |> deg2rad; tolerance = (1, 1)
    )
    @test c == true
    @test a ≈ a_ref rtol = 1e-2
    @test i ≈ i_ref rtol = 1e-2
    @test !isapprox(a, a_ref; rtol = 1e-12)

    # The keyword `no_warnings` must suppress all the warnings.
    a, i, c = @test_logs sun_sync_orbit_from_angular_velocity(
        0.06 |> deg2rad; max_iterations = 3, no_warnings = true
    )
    @test c == false

    a, i, c = @test_logs sun_sync_orbit_from_angular_velocity(
        0.2 |> deg2rad; no_warnings = true
    )
    @test c == true

    # == Test When the Orbit Is Not Valid ==================================================

    a, i, c = (@test_logs(
        (:warn, "The orbit is not valid because the perigee is inside the Earth."),
        sun_sync_orbit_from_angular_velocity(0.2 |> deg2rad)
    ))
    @test a ≈ 3.1843358902077232e6
    @test i ≈ 1.5795229678241955
    @test c == true
end

@testset "Function sun_sync_orbit_from_angular_velocity" begin
    @test_throws ArgumentError sun_sync_orbit_from_angular_velocity(-0.01)
    @test_throws ArgumentError sun_sync_orbit_from_angular_velocity(+0.06 |> deg2rad, -0.01)
    @test_throws ArgumentError sun_sync_orbit_from_angular_velocity(+0.06 |> deg2rad, -1.1)
    @test_throws ArgumentError sun_sync_orbit_from_angular_velocity(+0.004 |> deg2rad, -1.1)
end

# -- Function: _sun_sync_orbit__residues_and_jacobian --------------------------------------

@testset "Function _sun_sync_orbit__residues_and_jacobian" begin
    # The analytical Jacobian must match the one obtained using central finite differences.
    # We use constants with the same order of magnitude as those in a LEO design, but with a
    # large `k₂` and `k₆` to make the second-order terms relevant.
    k = (-5.0, 0.2, 0.004, 0.0045, 3.5, 1e-3)
    Ω̇_d = 0.9856
    ω_d = 3.6

    fun(x, c) =
        SatelliteAnalysis._sun_sync_orbit__residues_and_jacobian(x, c, Ω̇_d, ω_d, k...)

    for (x, c) in ((0.95, -0.15), (0.9, 0.3), (0.8, -0.6), (1.0, 0.05))
        ~, J = fun(x, c)

        Δ = 1e-6
        ∂f_∂x = (fun(x + Δ, c)[1] - fun(x - Δ, c)[1]) / 2Δ
        ∂f_∂c = (fun(x, c + Δ)[1] - fun(x, c - Δ)[1]) / 2Δ

        @test J[1, 1] ≈ ∂f_∂x[1] rtol = 1e-7
        @test J[2, 1] ≈ ∂f_∂x[2] rtol = 1e-7
        @test J[1, 2] ≈ ∂f_∂c[1] rtol = 1e-7
        @test J[2, 2] ≈ ∂f_∂c[2] rtol = 1e-7
    end
end

# -- Function: sun_sync_orbit_inclination --------------------------------------------------

@testset "Function sun_sync_orbit_inclination" begin
    # == Float64 ===========================================================================

    i, c = sun_sync_orbit_inclination(7130.982e3)
    @test i isa Float64
    @test i ≈ 1.7175896973066611
    @test c == true

    i, c = sun_sync_orbit_inclination(7130.982e3, 0.05)
    @test i isa Float64
    @test i ≈ 1.716851774960272
    @test c == true

    # == Float32 ===========================================================================

    i, c = sun_sync_orbit_inclination(7130.982f3)
    @test i isa Float32
    @test i ≈ 1.7175897f0
    @test c == true

    i, c = sun_sync_orbit_inclination(7130.982f3, 0.05f0)
    @test i isa Float32
    @test i ≈ 1.7168517f0
    @test c == true

    # == Test When Algorithm Did Not Converge ==============================================

    a, c = (@test_logs(
        (:warn,), sun_sync_orbit_inclination(7130.982e3, 0.05; max_iterations = 2)
    ))
    @test c == false

    # A loose tolerance must lead to a less accurate result without warnings.
    i_ref, ~ = sun_sync_orbit_inclination(7130.982e3, 0.05)
    i, c = @test_logs sun_sync_orbit_inclination(7130.982e3, 0.05; tolerance = 1e-2)
    @test c == true
    @test i ≈ i_ref rtol = 1e-3
    @test !isapprox(i, i_ref; rtol = 1e-12)

    # The keyword `no_warnings` must suppress the warning.
    a, c = @test_logs sun_sync_orbit_inclination(
        7130.982e3, 0.05; max_iterations = 2, no_warnings = true
    )
    @test c == false
end

@testset "Function sun_sync_orbit_inclination [ERRORS]" begin
    @test_throws ArgumentError sun_sync_orbit_inclination(7130.982e3, 1.1)
    @test_throws ArgumentError sun_sync_orbit_inclination(7130.982e3, -0.1)
    @test_throws ArgumentError sun_sync_orbit_inclination(7130.982e3, 0.4)

    # == Test When the Orbit Is Not Valid ==================================================

    @test_throws ArgumentError sun_sync_orbit_inclination(15_000e3)
end

# -- Function: sun_sync_orbit_semi_major_axis ----------------------------------------------

@testset "Function sun_sync_orbit_semi_major_axis" begin
    # == Float64 ===========================================================================

    a, c = sun_sync_orbit_semi_major_axis(98.410 |> deg2rad)
    @test a isa Float64
    @test a ≈ 7.130827866508739e6
    @test c == true

    a, c = sun_sync_orbit_semi_major_axis(98.410 |> deg2rad, 0.05)
    @test a isa Float64
    @test a ≈ 7.141033706031902e6
    @test c == true

    # == Float32 ===========================================================================

    a, c = sun_sync_orbit_semi_major_axis(98.410f0 |> deg2rad)
    @test a isa Float32
    @test a ≈ 7.1308275f6
    @test c == true

    a, c = sun_sync_orbit_semi_major_axis(98.410f0 |> deg2rad, 0.05f0)
    @test a isa Float32
    @test a ≈ 7.1410335f6
    @test c == true

    # == Test When Algorithm Did Not Converge ==============================================

    a, c = (@test_logs(
        (:warn,),
        sun_sync_orbit_semi_major_axis(98.410 |> deg2rad, 0.05; max_iterations = 2)
    ))
    @test c == false

    # A loose tolerance must lead to a less accurate result without warnings.
    a_ref, ~ = sun_sync_orbit_semi_major_axis(98.410 |> deg2rad, 0.05)
    a, c = @test_logs sun_sync_orbit_semi_major_axis(
        98.410 |> deg2rad, 0.05; tolerance = 1e-2
    )
    @test c == true
    @test a ≈ a_ref rtol = 1e-3
    @test !isapprox(a, a_ref; rtol = 1e-12)

    # The keyword `no_warnings` must suppress all the warnings.
    a, c = @test_logs sun_sync_orbit_semi_major_axis(
        98.410 |> deg2rad, 0.05; max_iterations = 2, no_warnings = true
    )
    @test c == false

    a, c = @test_logs sun_sync_orbit_semi_major_axis(90.01 |> deg2rad; no_warnings = true)
    @test c == true

    # == Test When the Orbit Is Not Valid ==================================================

    a, c = (@test_logs(
        (:warn, "The orbit is not valid because the perigee is inside the Earth."),
        sun_sync_orbit_semi_major_axis(95.5 |> deg2rad)
    ))
    @test a ≈ 6.319396335751179e6
    @test c == true
end

@testset "Function sun_sync_orbit_semi_major_axis [ERRORS]" begin
    @test_throws ArgumentError sun_sync_orbit_semi_major_axis(0.5, 1.1)
    @test_throws ArgumentError sun_sync_orbit_semi_major_axis(0.5, -0.01)

    # There is no Sun-synchronous orbit if the inclination is not higher than 90°.
    @test_throws ArgumentError sun_sync_orbit_semi_major_axis(π / 2)
    @test_throws ArgumentError sun_sync_orbit_semi_major_axis(80 |> deg2rad)
end
