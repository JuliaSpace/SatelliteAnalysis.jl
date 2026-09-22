## Description #############################################################################
#
# Tests related to the beta angle analysis.
#
############################################################################################

# == File: ./src/beta_angle.jl =============================================================

# -- Function: beta_angle ------------------------------------------------------------------

@testset "Function beta_angle" begin
    jd₀ = SatelliteAnalysis.date_to_jd(2021, 1, 1, 0, 0, 0)
    orb = KeplerianElements(
        jd₀,
        7130.982e3,
        0.001111,
        98.405 |> deg2rad,
        ltdn_to_raan(10.5, jd₀),
        90 |> deg2rad,
        0,
    )

    β = beta_angle(orb, 5)
    @test β ≈ 0.44026044605171655

    β_max = beta_angle.(orb, 0:1:364) |> maximum
    @test β_max ≈ 0.4746822398784829

    # == Regression Values =================================================================

    # Those values were obtained with the algorithm that computed the orbit normal using a
    # DCM. They cover all the perturbation models, polar orbits, and equatorial orbits.
    jd₁ = SatelliteAnalysis.date_to_jd(2021, 3, 7, 5, 30, 0)

    for (i, Ω, Δjd, perturbation, expected) in (
        (0.3, 1.1, 17.5, :J0, +0.2857419956622166),
        (1.7175, 3.0, 200, :J2, -0.1472939222910801),
        (2.9, 5.9, 364, :J4, -0.0658908106853518),
        (π / 2, 0.0, 0, :J2, +0.2100868619654301),
        (0.0, 1.1, 200, :J4, -0.0028774182089317),
    )
        orb = KeplerianElements(jd₁, 7130.982e3, 0.01, i, Ω, 1.0, 0.5)
        @test beta_angle(orb, Δjd; perturbation = perturbation) ≈ expected atol = 1e-14
    end
end
