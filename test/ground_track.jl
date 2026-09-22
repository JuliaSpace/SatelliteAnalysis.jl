## Description #############################################################################
#
# Tests related to the ground tracks.
#
############################################################################################

# == File: ./src/ground_track.jl ===========================================================

# -- Function: ground_track ----------------------------------------------------------------

@testset "Function ground_track" begin
    jd₀ = SatelliteAnalysis.date_to_jd(2024, 1, 1)
    orb = KeplerianElements(
        jd₀,
        7130.982e3,
        0.001111,
        98.405 |> deg2rad,
        ltdn_to_raan(10.5, jd₀),
        90 |> deg2rad,
        0,
    )
    orbp = Propagators.init(Val(:J2), orb)

    # == Scenario 01 =======================================================================

    gt = ground_track(orbp; duration = 30)

    @test length(gt) == 2
    @test gt[1][1] ≈ 1.4249660489382514
    @test gt[1][2] ≈ -1.9629472118375995
    @test gt[2][1] ≈ 1.4239401788695614
    @test gt[2][2] ≈ -2.0832752921442603

    # == Scenario 02 =======================================================================

    gt = ground_track(orbp; step = 1000, duration = 1000)

    @test length(gt) == 3
    @test gt[1][1] ≈ 1.4249660489382514
    @test gt[1][2] ≈ -1.9629472118375995
    @test isnan(gt[2][1]) == true
    @test isnan(gt[2][2]) == true
    @test gt[3][1] ≈ 0.5180892580826623
    @test gt[3][2] ≈ 2.76053859491629

    gt = ground_track(orbp; step = 1000, duration = 1000, add_nans = false)

    @test length(gt) == 2
    @test gt[1][1] ≈ 1.4249660489382514
    @test gt[1][2] ≈ -1.9629472118375995
    @test gt[2][1] ≈ 0.5180892580826623
    @test gt[2][2] ≈ 2.76053859491629

    gt = ground_track(orbp; step = 1000, duration = 1000, track_types = :ascending)

    @test isempty(gt) == true

    @test_throws ArgumentError ground_track(orbp; track_types = :unknown)
    @test_throws ArgumentError ground_track(orbp; step = 0)
    @test_throws ArgumentError ground_track(orbp; step = -10)

    # == Keywords initial_time and f_eci_to_ecef ===========================================

    # The ground track computed from `initial_time` must be equal to the corresponding part
    # of the ground track computed from the epoch.
    gt_ref = ground_track(orbp; step = 500, duration = 3000, add_nans = false)
    gt_ini = ground_track(
        orbp; step = 500, duration = 2000, initial_time = 1000, add_nans = false
    )

    @test length(gt_ini) == 5
    @test all(gt_ini[k] == gt_ref[k + 2] for k in 1:5)

    # A user-defined conversion from ECI to ECEF must be used. If it is the identity, the
    # longitude drifts with respect to the default conversion, but the latitude is almost
    # the same.
    gt_id = ground_track(
        orbp; step = 500, duration = 3000, add_nans = false, f_eci_to_ecef = (r, jd) -> r
    )

    @test length(gt_id) == length(gt_ref)
    @test all(abs(gt_id[k][1] - gt_ref[k][1]) < 0.01 for k in eachindex(gt_ref))
    @test !all(abs(gt_id[k][2] - gt_ref[k][2]) < 0.01 for k in eachindex(gt_ref))

    # == Passage Separation in Low-Inclination Orbits ======================================

    # If the inclination is lower than 45°, the latitude difference between the end of a
    # passage and the beginning of the next one is lower than 90°. Hence, the passages must
    # be separated by NaNs using the information that points were skipped. In this case,
    # the latitude must be monotonic between two consecutive valid points.
    orb_li = KeplerianElements(
        jd₀, 7130.982e3, 0.001111, 30 |> deg2rad, 0, 90 |> deg2rad, 0
    )
    orbp_li = Propagators.init(Val(:J2), orb_li)

    for (track_types, cmp) in ((:ascending, >), (:descending, <))
        gt = ground_track(orbp_li; track_types = track_types)

        @test count(p -> isnan(p[1]), gt) >= 13

        @test all(
            isnan(gt[k][1]) || isnan(gt[k - 1][1]) || cmp(gt[k][1], gt[k - 1][1]) for
            k in 2:length(gt)
        )

        # Without NaNs, the vector must contain exactly the same valid points.
        gt_no_nans = ground_track(orbp_li; track_types = track_types, add_nans = false)

        @test gt_no_nans == filter(p -> !isnan(p[1]), gt)
    end

    # The output element type follows promoted floating-point time inputs rather than being
    # fixed to Float64.
    gt = ground_track(orbp; step = big"1000", duration = big"1000")
    @test eltype(gt) == NTuple{2, BigFloat}
end

# -- Function: ground_track_inclination ----------------------------------------------------

@testset "Function ground_track_inclination" begin
    jd₀ = SatelliteAnalysis.date_to_jd(2021, 1, 1)
    orb = KeplerianElements(
        jd₀,
        7130.982e3,
        0.001111,
        98.410 |> deg2rad,
        ltdn_to_raan(10.5, jd₀),
        90 |> deg2rad,
        0,
    )

    i_gt = ground_track_inclination(orb) |> rad2deg

    @test i_gt ≈ 102.30052101658998

    # Test if the parameters selection is working properly by comparing the J₀ case with the
    # J₄ case but setting J₂ and J₄ to 0.
    i_gt_J₀ = ground_track_inclination(orb; perturbation = :J0)
    i_gt    = ground_track_inclination(orb; perturbation = :J4, J2 = 0, J4 = 0)
    @test i_gt == i_gt_J₀
end
