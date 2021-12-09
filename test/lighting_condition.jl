# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==============================================================================
#
#   Tests related to lighting condition analysis.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# File: ./src/lighting_condition.jl
# ==============================================================================

# Function: lighting_condition
# ------------------------------------------------------------------------------

@testset "Function lighting_condition" begin
    jd = date_to_jd(2021, 12, 8, 0, 0, 25)
    s_i = sun_position_i(jd)

    r_i = SVector(1.06746e6, 3.22132e6, -6.2841e6)
    lc = lighting_condition(r_i, s_i)
    @test lc === :sunlight

    r_i = SVector(2.52449e6, 4.89401e6, -4.54648e6)
    lc = lighting_condition(r_i, s_i)
    @test lc === :penumbra

    r_i = SVector(3.9575e6, 5.91988e6, -489577.0)
    lc = lighting_condition(r_i, s_i)
    @test lc === :umbra

    lc = lighting_condition(s_i / norm(s_i) * 7000e3, s_i)
    @test lc === :sunlight
end
