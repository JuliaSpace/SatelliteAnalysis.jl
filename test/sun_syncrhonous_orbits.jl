# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#   Tests related to the Sun-synchronous orbit functions.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# File: ./src/sun_synchronous_orbits.jl
# ==========================================================================================

# Function: design_sun_sync_ground_repeating_orbits
# ------------------------------------------------------------------------------------------

@testset "Designing Sun-syncrhonous, ground-repeating orbits" begin

    # Default
    # ======================================================================================

    df = design_sun_sync_ground_repeating_orbit(1, 1)

    @test size(df) == (5, 7)

    @test df[begin, :semi_major_axis]       ≈ 6382.409  (atol = 1e-3)
    @test df[begin, :altitude]              ≈ 4.272     (atol = 1e-3)
    @test df[begin, :inclination]           ≈ 95.6949   (atol = 1e-4)
    @test df[begin, :period]                ≈ 84.706    (atol = 1e-3)
    @test df[begin, :rev_per_days]         == "17"
    @test df[begin, :adjacent_gt_distance]  ≈ 2345.163  (atol = 1e-3)
    @test df[begin, :adjacent_gt_angle]     ≈ 169.0503  (atol = 1e-4)

    @test df[end,   :semi_major_axis]       ≈ 7635.252  (atol = 1e-3)
    @test df[end,   :altitude]              ≈ 1257.115  (atol = 1e-3)
    @test df[end,   :inclination]           ≈ 100.7057  (atol = 1e-4)
    @test df[end,   :period]                ≈ 110.769   (atol = 1e-3)
    @test df[end,   :rev_per_days]         == "13"
    @test df[end,   :adjacent_gt_distance]  ≈ 3024.567  (atol = 1e-3)
    @test df[end,   :adjacent_gt_angle]     ≈ 92.4444   (atol = 1e-4)

    # Altitude filter
    # ======================================================================================

    df = design_sun_sync_ground_repeating_orbit(
        5,
        5;
        minimum_altitude = 750e3,
        maximum_altitude = 760e3
    )

    @test size(df) == (1, 7)
    @test df[begin, :semi_major_axis]       ≈ 7130.984   (atol = 1e-3)
    @test df[begin, :altitude]              ≈ 752.847    (atol = 1e-3)
    @test df[begin, :inclination]           ≈ 98.4106    (atol = 1e-4)
    @test df[begin, :period]                ≈ 100.000    (atol = 1e-3)
    @test df[begin, :rev_per_days]         == "14 + ²/₅"
    @test df[begin, :adjacent_gt_distance]  ≈ 550.597    (atol = 1e-3)
    @test df[begin, :adjacent_gt_angle]     ≈ 39.872     (atol = 1e-3)

    # Revolution per days
    # ======================================================================================

    df = design_sun_sync_ground_repeating_orbit(
        5,
        5;
        minimum_altitude = 750e3,
        maximum_altitude = 760e3,
        pretify_rev_per_days = false
    )

    @test size(df) == (1, 7)
    @test df[begin, :semi_major_axis]       ≈ 7130.984     (atol = 1e-3)
    @test df[begin, :altitude]              ≈ 752.847      (atol = 1e-3)
    @test df[begin, :inclination]           ≈ 98.4106      (atol = 1e-4)
    @test df[begin, :period]                ≈ 100.000      (atol = 1e-3)
    @test df[begin, :rev_per_days]         == (14, 2 // 5)
    @test df[begin, :adjacent_gt_distance]  ≈ 550.597      (atol = 1e-3)
    @test df[begin, :adjacent_gt_angle]     ≈ 39.872       (atol = 1e-3)

    # Test units
    # ======================================================================================

    # Angle
    # --------------------------------------------------------------------------------------

    df = design_sun_sync_ground_repeating_orbit(
        5,
        5;
        minimum_altitude = 750e3,
        maximum_altitude = 760e3,
        angle_unit = :rad,
    )

    @test df[begin, :inclination]       ≈ 98.4106 |> deg2rad (atol = 2e-6)
    @test df[begin, :adjacent_gt_angle] ≈ 39.872  |> deg2rad (atol = 8e-5)

    # Distance
    # --------------------------------------------------------------------------------------

    df = design_sun_sync_ground_repeating_orbit(
        5,
        5;
        minimum_altitude = 750e3,
        maximum_altitude = 760e3,
        distance_unit = :m,
    )

    @test df[begin, :semi_major_axis]      ≈ 7130.984e3 (atol = 1)
    @test df[begin, :altitude]             ≈ 752.847e3  (atol = 1)
    @test df[begin, :adjacent_gt_distance] ≈ 550.597e3  (atol = 1)

    # Time
    # --------------------------------------------------------------------------------------

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
    # ======================================================================================

    @test_throws ArgumentError design_sun_sync_ground_repeating_orbit(0, 10)
    @test_throws ArgumentError design_sun_sync_ground_repeating_orbit(3, -3)
    @test_throws ArgumentError design_sun_sync_ground_repeating_orbit(3, 2)
    @test_throws ArgumentError design_sun_sync_ground_repeating_orbit(1, 5; e = 1)
end

# Function: sun_sync_orbit_from_angular_velocity
# ------------------------------------------------------------------------------------------

@testset "Function sun_sync_orbit_from_angular_velocity" begin
    # Float64
    # ======================================================================================

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

    # Float32
    # ======================================================================================

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

    # Test When Algorithm Did Not Converged
    # ======================================================================================

    a, i, c = (@test_logs(
        (:warn,),
        sun_sync_orbit_from_angular_velocity(0.06 |> deg2rad; max_iterations = 3)
    ))
    @test c == false

    # Test When the Orbit Is Not Valid
    # ======================================================================================

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
    @test_throws ArgumentError sun_sync_orbit_from_angular_velocity(+0.06  |> deg2rad, -0.01)
    @test_throws ArgumentError sun_sync_orbit_from_angular_velocity(+0.06  |> deg2rad, -1.1)
    @test_throws ArgumentError sun_sync_orbit_from_angular_velocity(+0.004 |> deg2rad, -1.1)
end

# Function: sun_sync_orbit_inclination
# ------------------------------------------------------------------------------------------

@testset "Function sun_sync_orbit_inclination" begin
    # Float64
    # ======================================================================================

    i, c = sun_sync_orbit_inclination(7130.982e3)
    @test i isa Float64
    @test i ≈ 1.7175896973066611
    @test c == true

    i, c = sun_sync_orbit_inclination(7130.982e3, 0.05)
    @test i isa Float64
    @test i ≈ 1.716851774960272
    @test c == true

    # Float32
    # ======================================================================================

    i, c = sun_sync_orbit_inclination(7130.982f3)
    @test i isa Float32
    @test i ≈ 1.7175897f0
    @test c == true

    i, c = sun_sync_orbit_inclination(7130.982f3, 0.05f0)
    @test i isa Float32
    @test i ≈ 1.7168517f0
    @test c == true

    # Test When Algorithm Did Not Converge
    # ======================================================================================

    a, c = (@test_logs(
        (:warn,),
        sun_sync_orbit_inclination(7130.982e3, 0.05; max_iterations = 2)
    ))
    @test c == false
end

@testset "Function sun_sync_orbit_inclination [ERRORS]" begin
    @test_throws ArgumentError sun_sync_orbit_inclination(7130.982e3,  1.1)
    @test_throws ArgumentError sun_sync_orbit_inclination(7130.982e3, -0.1)
    @test_throws ArgumentError sun_sync_orbit_inclination(7130.982e3,  0.4)

    # Test When the Orbit Is Not Valid
    # ======================================================================================

    @test_throws ArgumentError sun_sync_orbit_inclination(15_000e3)
end

# Function: sun_sync_orbit_semi_major_axis
# ------------------------------------------------------------------------------------------

@testset "Function sun_sync_orbit_semi_major_axis" begin
    # Float64
    # ======================================================================================

    a, c = sun_sync_orbit_semi_major_axis(98.410 |> deg2rad)
    @test a isa Float64
    @test a ≈ 7.130827866508739e6
    @test c == true

    a, c = sun_sync_orbit_semi_major_axis(98.410 |> deg2rad, 0.05)
    @test a isa Float64
    @test a ≈ 7.141033706031902e6
    @test c == true

    # Float32
    # ======================================================================================

    a, c = sun_sync_orbit_semi_major_axis(98.410f0 |> deg2rad)
    @test a isa Float32
    @test a ≈ 7.1308275f6
    @test c == true

    a, c = sun_sync_orbit_semi_major_axis(98.410f0 |> deg2rad, 0.05f0)
    @test a isa Float32
    @test a ≈ 7.1410335f6
    @test c == true

    # Test When Algorithm Did Not Converge
    # ======================================================================================

    a, c = (@test_logs(
        (:warn,),
        sun_sync_orbit_semi_major_axis(98.410 |> deg2rad, 0.05; max_iterations = 2)
    ))
    @test c == false

    # Test When the Orbit Is Not Valid
    # ======================================================================================

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
end
