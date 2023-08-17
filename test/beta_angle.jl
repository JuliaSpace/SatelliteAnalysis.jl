# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#   Tests related to the beta angle analysis.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# File: ./src/beta_angle.jl
# ==========================================================================================

# Function: beta_angle
# ------------------------------------------------------------------------------------------

@testset "Function beta_angle" begin
    jd₀ = SatelliteAnalysis.date_to_jd(2021, 1, 1, 0, 0, 0)
    orb = KeplerianElements(
        jd₀,
        7130.982e3,
        0.001111,
        98.405 |> deg2rad,
        ltdn_to_raan(10.5, jd₀),
        90     |> deg2rad,
        0
    )

    β = beta_angle(orb, 5)
    @test β ≈ 0.44026044605171655

    β_max = beta_angle.(orb, 0:1:364) |> maximum
    @test β_max ≈ 0.4746822398784829
end
