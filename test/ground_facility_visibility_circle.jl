## Description #############################################################################
#
# Tests related to the ground facility visibility circle.
#
############################################################################################

# == File: ./src/ground_facilities/is_ground_facility_visible.jl ===========================

@testset "Function is_ground_facility_visible" verbose = true begin
    ret = is_ground_facility_visible([7000e3, 0, 0], [6378e3, 0, 0], 10 |> deg2rad)
    @test ret == true

    ret = is_ground_facility_visible([6500e3, 700e3, 0], [6378e3, 0, 0], 10 |> deg2rad)
    @test ret == false

    # == Method With the Ground Facility Position and Local Vertical =======================

    # This method must provide the same result as the one that receives the geodetic
    # coordinates for satellites around the entire Earth.
    for (gf_lat, gf_lon, gf_h) in ((0.0, 0.0, 0.0), (-0.4, 2.5, 800.0), (1.3, -1.0, 50.0))
        gf_r_e  = geodetic_to_ecef(gf_lat, gf_lon, gf_h)
        gf_up_e = [cos(gf_lat) * cos(gf_lon), cos(gf_lat) * sin(gf_lon), sin(gf_lat)]

        num_visible = 0

        for lat in -1.5:0.1:1.5, lon in -3.1:0.1:3.1, θ in (0.0, 0.2, 0.5)
            sat_r_e = geodetic_to_ecef(lat, lon, 700e3)

            expected = is_ground_facility_visible(sat_r_e, gf_lat, gf_lon, gf_h, θ)
            result   = is_ground_facility_visible(sat_r_e, gf_r_e, gf_up_e, θ)

            num_visible += expected

            @test result == expected
        end

        @test num_visible > 0
    end

    # == Method With the Satellite Position in the NED Reference Frame =====================

    @test is_ground_facility_visible([0.0, 0.0, -700e3], 10 |> deg2rad) == true
    @test is_ground_facility_visible([700e3, 0.0, -10e3], 10 |> deg2rad) == false
end

# == File: ./src/ground_facilities/ground_facility_visibility_circle.jl ====================

# -- Function: ground_facility_visibility_circle -------------------------------------------

@testset "Function ground_facility_visibility_circle" verbose = true begin
    gfv = ground_facility_visibility_circle((0, 0, 0), 7000e3; azimuth_step = π / 2)

    @test length(gfv) == 6

    @test gfv[1][1] ≈ -0.2842464765581056
    @test gfv[1][2] ≈ 0.0 atol = 1e-10
    @test gfv[2][1] ≈ 0.0 atol = 1e-10
    @test gfv[2][2] ≈ -0.28260412678209323
    @test gfv[3][1] ≈ +0.28424647655810553
    @test gfv[3][2] ≈ 0.0 atol = 1e-10
    @test gfv[4][1] ≈ 0.0 atol = 1e-10
    @test gfv[4][2] ≈ +0.28260412678209323

    # The last two are equal to the first two to complete the circle.
    @test gfv[5][1] ≈ gfv[1][1] atol = 1e-10
    @test gfv[5][2] ≈ gfv[1][2] atol = 1e-10
    @test gfv[6][1] ≈ gfv[2][1] atol = 1e-10
    @test gfv[6][2] ≈ gfv[2][2] atol = 1e-10

    # The satellite must be above the ground facility.
    @test_throws ArgumentError ground_facility_visibility_circle((0, 0, 0), 6000e3)
    @test_throws ArgumentError ground_facility_visibility_circle(
        (0, 0, 0), EARTH_EQUATORIAL_RADIUS
    )
end
