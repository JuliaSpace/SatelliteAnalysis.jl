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
 
    @test gfv[1][1] ≈ -0.7002878543070394
    @test gfv[1][2] ≈ -1.0276006774869804e-16
    @test gfv[2][1] ≈ 3.2561543045931737e-17
    @test gfv[2][2] ≈ -0.6981317007977318
    @test gfv[3][1] ≈ 0.7002878543070393
    @test gfv[3][2] ≈ 0.0
    @test gfv[4][1] ≈ 3.2561543045931737e-17
    @test gfv[4][2] ≈ 0.6981317007977318

    # The last two are equal to the first two to complete the circle.
    @test gfv[5][1] ≈ gfv[1][1] atol = 1e-10
    @test gfv[5][2] ≈ gfv[1][2] atol = 1e-10
    @test gfv[6][1] ≈ gfv[2][1] atol = 1e-10
    @test gfv[6][2] ≈ gfv[2][2] atol = 1e-10
end
