## Description #############################################################################
#
# Tests related to the ground facility visibility circle.
#
############################################################################################

# == File: ./src/ground_facilities/ground_facility_visibility_circle.jl ====================

# -- Function: ground_facility_visibility_circle -------------------------------------------

@testset "Function ground_facility_visibility_circle" verbose = true begin
    gfv = ground_facility_visibility_circle(
        (0, 0, 0),
        7000e3;
        azimuth_step = π / 2
    )

    @test length(gfv) == 6

    @test gfv[1][1] ≈ -0.2842464765581056
    @test gfv[1][2] ≈ -3.556079170773591e-17
    @test gfv[2][1] ≈  1.4150606804235496e-17
    @test gfv[2][2] ≈ -0.28260412678209323
    @test gfv[3][1] ≈  0.28424647655810553
    @test gfv[3][2] ≈  0.0
    @test gfv[4][1] ≈  1.4150606804235496e-17
    @test gfv[4][2] ≈  0.28260412678209323

    # The last two are equal to the first two to complete the circle.
    @test gfv[5][1] ≈ gfv[1][1] atol = 1e-10
    @test gfv[5][2] ≈ gfv[1][2] atol = 1e-10
    @test gfv[6][1] ≈ gfv[2][1] atol = 1e-10
    @test gfv[6][2] ≈ gfv[2][2] atol = 1e-10
end
