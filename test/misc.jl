# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==============================================================================
#
#   Tests related to the miscellaneous functions.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# File: ./src/misc/miscellaneous.jl
# ==============================================================================

# Function: find_crossing
# ------------------------------------------------------------------------------

@testset "Function find_crossing" begin
    f(t, a) = sin(t - a) ≥ 0

    for a in (-0.5 .+ 0.5rand(3))
        t = SatelliteAnalysis.find_crossing(
            f,
            a - 0.1,
            a + 0.1,
            false,
            true,
            a;
            Δ = 1e-20
        )

        @test t ≈ a
    end
end

@testset "Function find_crossing (errors)" begin
    f(t, a) = :test

    @test_throws Exception SatelliteAnalysis.find_crossing(
        f,
        -0.1,
        +0.1,
        false,
        true,
        a;
        Δ = 1e-20
    )
end

