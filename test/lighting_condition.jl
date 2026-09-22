## Description #############################################################################
#
# Tests related to lighting condition analysis.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# == File: ./src/lighting_condition.jl =====================================================

# -- Function: lighting_condition ----------------------------------------------------------

@testset "Function lighting_condition" begin
    jd = date_to_jd(2021, 12, 8, 0, 0, 25)
    s_i = sun_position_mod(jd)

    r_i = SVector(1.06746e6, 3.22132e6, -6.2841e6)
    lc = lighting_condition(r_i, s_i)
    @test lc == :sunlight

    r_i = SVector(2.52449e6, 4.89401e6, -4.54648e6)
    lc = lighting_condition(r_i, s_i)
    @test lc == :penumbra

    r_i = SVector(3.9575e6, 5.91988e6, -489577.0)
    lc = lighting_condition(r_i, s_i)
    @test lc == :umbra

    lc = lighting_condition(s_i / norm(s_i) * 7000e3, s_i)
    @test lc == :sunlight

    # == Shadow Boundaries =================================================================

    # Regression test for the shadow boundaries at 7000 km behind the Earth as a function of
    # the distance `ρ` to the shadow axis. The reference values were obtained using the
    # algorithm with trigonometric functions: the umbra ends between 6345.2 km and
    # 6345.4 km, and the penumbra ends between 6411.6 km and 6411.8 km.
    jd  = date_to_jd(2021, 1, 1)
    s_i = sun_position_mod(jd)
    s̄_i = s_i / norm(s_i)
    p̄_i = normalize(SVector(0.0, 0.0, 1.0) - s̄_i[3] * s̄_i)

    lc(ρ) = lighting_condition(-7000e3 * s̄_i + ρ * p̄_i, s_i)

    @test lc(0.0) == :umbra
    @test lc(6345.2e3) == :umbra
    @test lc(6345.4e3) == :penumbra
    @test lc(6411.6e3) == :penumbra
    @test lc(6411.8e3) == :sunlight
    @test lc(8000e3) == :sunlight
end
