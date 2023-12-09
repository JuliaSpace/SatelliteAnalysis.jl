## Description #############################################################################
#
# Tests related to the ground tracks.
#
############################################################################################

# == File: ./src/ground_trakc.jl ===========================================================

# -- Function: ground_track ----------------------------------------------------------------

@testset "Function ground_track" begin
    jd₀ = SatelliteAnalysis.date_to_jd(2024, 1, 1)
    orb = KeplerianElements(
        jd₀,
        7130.982e3,
        0.001111,
        98.405 |> deg2rad,
        ltdn_to_raan(10.5, jd₀),
        90     |> deg2rad,
        0
    )

    # == Scenario 01 =======================================================================

    gt = ground_track(orbp; duration = 30)

    @test length(gt) == 2
    @test gt[1][1] ≈  1.4249660489382514
    @test gt[1][2] ≈ -1.9629472118375995
    @test gt[2][1] ≈  1.4239401788695614
    @test gt[2][2] ≈ -2.0832752921442603

    # == Scenario 02 =======================================================================

    gt = ground_track(orbp; step = 1000, duration = 1000)

    @test length(gt) == 3
    @test gt[1][1] ≈  1.4249660489382514
    @test gt[1][2] ≈ -1.9629472118375995
    @test isnan(gt[2][1]) == true
    @test isnan(gt[2][2]) == true
    @test gt[3][1] ≈ 0.5180892580826623
    @test gt[3][2] ≈ 2.76053859491629

    gt = ground_track(orbp; step = 1000, duration = 1000, add_nans = false)

    @test length(gt) == 2
    @test gt[1][1] ≈  1.4249660489382514
    @test gt[1][2] ≈ -1.9629472118375995
    @test gt[2][1] ≈ 0.5180892580826623
    @test gt[2][2] ≈ 2.76053859491629

    gt = ground_track(orbp; step = 1000, duration = 1000, track_types = :ascending)

    @test isempty(gt) == true
end

