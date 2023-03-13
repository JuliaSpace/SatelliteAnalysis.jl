# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==============================================================================
#
#   Tests related to the Sun syncrhonous orbit functions.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# File: ./src/sun_synchronous_orbits.jl
# ==============================================================================

# Function: design_sun_sync_ground_repeating_orbits
# ------------------------------------------------------------------------------

@testset "Designing Sun-syncrhonous, ground-repeating orbits" begin

    # Default
    # ==========================================================================

    df = design_sun_sync_ground_repeating_orbit(1, 1)

    @test size(df) == (5, 7)

    @test df[begin, :semi_major_axis]       ≈ 6382.406  (atol = 1e-3)
    @test df[begin, :altitude]              ≈    4.269  (atol = 1e-3)
    @test df[begin, :inclination]           ≈   95.6904 (atol = 1e-4)
    @test df[begin, :period]                ≈   84.706  (atol = 1e-3)
    @test df[begin, :rev_per_days]         == "17"
    @test df[begin, :adjacent_gt_distance]  ≈ 2345.182  (atol = 1e-3)
    @test df[begin, :adjacent_gt_angle]     ≈ 84.5253   (atol = 1e-4)

    @test df[end,   :semi_major_axis]       ≈ 7635.250  (atol = 1e-3)
    @test df[end,   :altitude]              ≈ 1257.113  (atol = 1e-3)
    @test df[end,   :inclination]           ≈  100.7002 (atol = 1e-4)
    @test df[end,   :period]                ≈  110.769  (atol = 1e-3)
    @test df[end,   :rev_per_days]         == "13"
    @test df[end,   :adjacent_gt_distance]  ≈ 3024.626  (atol = 1e-3)
    @test df[end,   :adjacent_gt_angle]     ≈   46.2227 (atol = 1e-4)

    # Altitude filter
    # ==========================================================================

    df = design_sun_sync_ground_repeating_orbit(
        5,
        5;
        minimum_altitude = 750e3,
        maximum_altitude = 760e3
    )

    @test size(df) == (1, 7)
    @test df[begin, :semi_major_axis]       ≈ 7130.982  (atol = 1e-3)
    @test df[begin, :altitude]              ≈  752.845  (atol = 1e-3)
    @test df[begin, :inclination]           ≈   98.4055 (atol = 1e-4)
    @test df[begin, :period]                ≈  100.000  (atol = 1e-3)
    @test df[begin, :rev_per_days]         == "14 + ²/₅"
    @test df[begin, :adjacent_gt_distance]  ≈  550.604  (atol = 1e-3)
    @test df[begin, :adjacent_gt_angle]     ≈   19.936  (atol = 1e-4)

    # Revolution per days
    # ==========================================================================

    df = design_sun_sync_ground_repeating_orbit(
        5,
        5;
        minimum_altitude = 750e3,
        maximum_altitude = 760e3,
        pretify_rev_per_days = false
    )

    @test size(df) == (1, 7)
    @test df[begin, :semi_major_axis]       ≈ 7130.982     (atol = 1e-3)
    @test df[begin, :altitude]              ≈  752.845     (atol = 1e-3)
    @test df[begin, :inclination]           ≈   98.4055    (atol = 1e-4)
    @test df[begin, :period]                ≈  100.000     (atol = 1e-3)
    @test df[begin, :rev_per_days]         == (14, 2 // 5)
    @test df[begin, :adjacent_gt_distance]  ≈  550.604     (atol = 1e-3)
    @test df[begin, :adjacent_gt_angle]     ≈   19.936     (atol = 1e-4)

    # Test units
    # ==========================================================================

    # Angle
    # --------------------------------------------------------------------------

    df = design_sun_sync_ground_repeating_orbit(
        5,
        5;
        minimum_altitude = 750e3,
        maximum_altitude = 760e3,
        angle_unit = :rad,
    )

    @test df[begin, :inclination]       ≈ 98.4055 |> deg2rad (atol = 2e-6)
    @test df[begin, :adjacent_gt_angle] ≈ 19.936  |> deg2rad (atol = 2e-6)

    # Distance
    # --------------------------------------------------------------------------

    df = design_sun_sync_ground_repeating_orbit(
        5,
        5;
        minimum_altitude = 750e3,
        maximum_altitude = 760e3,
        distance_unit = :m,
    )

    @test df[begin, :semi_major_axis]      ≈ 7130.982e3 (atol = 1)
    @test df[begin, :altitude]             ≈ 752.845e3  (atol = 1)
    @test df[begin, :adjacent_gt_distance] ≈ 550.604e3  (atol = 1)

    # Time
    # --------------------------------------------------------------------------

    df = design_sun_sync_ground_repeating_orbit(
        5,
        5;
        minimum_altitude = 750e3,
        maximum_altitude = 760e3,
        time_unit = :s,
    )

    @test df[begin, :period] ≈ 6000.000 (atol = 1e-3)

    df = design_sun_sync_ground_repeating_orbit(
        5,
        5;
        minimum_altitude = 750e3,
        maximum_altitude = 760e3,
        time_unit = :h,
    )

    @test df[begin, :period] ≈ 1.66667 (atol = 1e-5)

    # Errors
    # ==========================================================================

    @test_throws ArgumentError design_sun_sync_ground_repeating_orbit(0, 10)
    @test_throws ArgumentError design_sun_sync_ground_repeating_orbit(3, -3)
    @test_throws ArgumentError design_sun_sync_ground_repeating_orbit(3, 2)
    @test_throws ArgumentError design_sun_sync_ground_repeating_orbit(1, 5; e = 1)
end
