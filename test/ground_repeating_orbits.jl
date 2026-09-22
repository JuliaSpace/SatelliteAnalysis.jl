## Description #############################################################################
#
# Tests related to the ground repeating orbits.
#
############################################################################################

# == File: ./src/ground_repeating_orbits.jl ================================================

# -- Functions: ground_repeating_orbit_adjacent_track_angle and _distance ------------------

@testset "Adjacent Track Angle and Distance" begin
    # Sun-synchronous orbit of the Amazonia-1 mission, which repeats the ground track after
    # 5 days with 14 + 2/5 revolutions per day.
    a = 7130.982e3
    e = 0.001111
    i = 98.410 |> deg2rad

    # == Default Keywords ==================================================================

    d = ground_repeating_orbit_adjacent_track_distance(a, e, i, 5)
    γ = ground_repeating_orbit_adjacent_track_angle(a, e, i, 5)

    @test d ≈ 543812.020224 atol = 1e-3
    @test γ ≈ 0.688106559027 atol = 1e-9

    # With a cycle of 1 day, the adjacent tracks are much farther apart.
    d = ground_repeating_orbit_adjacent_track_distance(a, e, i, 1)
    γ = ground_repeating_orbit_adjacent_track_angle(a, e, i, 1)

    @test d ≈ 2718135.854448 atol = 1e-3
    @test γ ≈ 1.967754895599 atol = 1e-9

    # == Perturbation Models ===============================================================

    d_J0 = ground_repeating_orbit_adjacent_track_distance(a, e, i, 5; perturbation = :J0)
    d_J2 = ground_repeating_orbit_adjacent_track_distance(a, e, i, 5; perturbation = :J2)
    d_J4 = ground_repeating_orbit_adjacent_track_distance(a, e, i, 5; perturbation = :J4)

    @test d_J0 ≈ 544640.853781 atol = 1e-3
    @test d_J2 ≈ 543812.020224 atol = 1e-3
    @test d_J4 ≈ 543813.291644 atol = 1e-3

    γ_J0 = ground_repeating_orbit_adjacent_track_angle(a, e, i, 5; perturbation = :J0)
    γ_J4 = ground_repeating_orbit_adjacent_track_angle(a, e, i, 5; perturbation = :J4)

    @test γ_J0 ≈ 0.689058920431 atol = 1e-9
    @test γ_J4 ≈ 0.688108020241 atol = 1e-9

    # == Constants =========================================================================

    # The keywords with the constants must change the result. Notice that the keyword `we`
    # was ignored in the versions before v0.4.0.
    d_we = ground_repeating_orbit_adjacent_track_distance(
        a, e, i, 5; we = 1.1 * EARTH_ANGULAR_SPEED
    )
    γ_we = ground_repeating_orbit_adjacent_track_angle(
        a, e, i, 5; we = 1.1 * EARTH_ANGULAR_SPEED
    )

    @test d_we ≈ 597450.791833 atol = 1e-3
    @test γ_we ≈ 0.748913296075 atol = 1e-9

    d_J2 = ground_repeating_orbit_adjacent_track_distance(
        a, e, i, 5; J2 = 1.1 * EGM_2008_J2
    )
    γ_J2 = ground_repeating_orbit_adjacent_track_angle(
        a, e, i, 5; J2 = 1.1 * EGM_2008_J2
    )

    @test d_J2 ≈ 543729.103408 atol = 1e-3
    @test γ_J2 ≈ 0.688011262719 atol = 1e-9

    # == Type Promotion ====================================================================

    d_32 = ground_repeating_orbit_adjacent_track_distance(
        Float32(a), Float32(e), Float32(i), 5
    )
    γ_32 = ground_repeating_orbit_adjacent_track_angle(
        Float32(a), Float32(e), Float32(i), 5
    )

    @test d_32 isa Float32
    @test γ_32 isa Float32
    @test d_32 ≈ 543812.02 rtol = 1e-5
    @test γ_32 ≈ 0.6881066 rtol = 1e-5

    @test ground_repeating_orbit_adjacent_track_distance(7130982, 0, i, 5) isa Float64

    # == Consistency With the Sun-Synchronous Orbit Design =================================

    df = design_sun_sync_ground_repeating_orbit(
        5,
        5;
        minimum_altitude = 750e3,
        maximum_altitude = 760e3,
        distance_unit    = :m,
        angle_unit       = :rad,
        eccentricity     = e,
    )

    @test size(df, 1) == 1

    a_df = df[begin, :semi_major_axis]
    i_df = df[begin, :inclination]

    @test df[begin, :adjacent_gt_distance] ≈
        ground_repeating_orbit_adjacent_track_distance(a_df, e, i_df, 5)
    @test df[begin, :adjacent_gt_angle] ≈
        ground_repeating_orbit_adjacent_track_angle(a_df, e, i_df, 5)
end
