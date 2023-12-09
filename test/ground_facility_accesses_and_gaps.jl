## Description #############################################################################
#
# Tests related to the ground facility accesses and gaps.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# == File: ./src/ground_facilities/ground_facility_accesses_and_gaps.jl ====================

# -- Function: ground_facility_accesses ----------------------------------------------------

@testset "Function ground_facility_accesses" begin
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
    orbp = Propagators.init(Val(:J2), orb)

    # == Single Facility ===================================================================

    df = ground_facility_accesses(orbp, (0, 0, 0) , 1 * 86400, TOD(), PEF())

    @test size(df) == (2, 3)

    @test df.access_beginning[1] == DateTime("2021-01-01T10:20:03.163")
    @test df.access_beginning[2] == DateTime("2021-01-01T22:49:55.898")

    @test df.access_end[1] == DateTime("2021-01-01T10:30:02.985")
    @test df.access_end[2] == DateTime("2021-01-01T22:59:23.524")

    exp_duration = getfield.(df.access_end .- df.access_beginning, :value) ./ 1000

    @test df.duration[1] ≈ exp_duration[1]
    @test df.duration[2] ≈ exp_duration[2]

    @test metadata(df, "Description") == "Accesses to the ground facilities."
    @test colmetadata(df, :duration, "Unit") == :s

    # == Multiple Facilities ===============================================================

    df = ground_facility_accesses(
        orbp,
        [(0, 0, 0), (0, π / 4, 0)],
        1 * 86400,
        TOD(),
        PEF()
    )

    @test size(df) == (4, 3)

    @test df.access_beginning[1] == DateTime("2021-01-01T07:00:36.953")
    @test df.access_beginning[2] == DateTime("2021-01-01T10:20:03.163")
    @test df.access_beginning[3] == DateTime("2021-01-01T19:29:57.333")
    @test df.access_beginning[4] == DateTime("2021-01-01T22:49:55.898")

    @test df.access_end[1] == DateTime("2021-01-01T07:10:04.095")
    @test df.access_end[2] == DateTime("2021-01-01T10:30:02.985")
    @test df.access_end[3] == DateTime("2021-01-01T19:39:56.981")
    @test df.access_end[4] == DateTime("2021-01-01T22:59:23.524")

    exp_duration = getfield.(df.access_end .- df.access_beginning, :value) ./ 1000

    @test df.duration[1] ≈ exp_duration[1]
    @test df.duration[2] ≈ exp_duration[2]
    @test df.duration[3] ≈ exp_duration[3]
    @test df.duration[4] ≈ exp_duration[4]

    @test metadata(df, "Description") == "Accesses to the ground facilities."
    @test colmetadata(df, :duration, "Unit") == :s

    # If we change the reduction function to AND instead of OR, we should have not access.
    df = ground_facility_accesses(
        orbp,
        [(0, 0, 0), (0, π / 4, 0)],
        1 * 86400,
        TOD(),
        PEF();
        reduction = v -> (&)(v...)
    )

    @test size(df) == (0, 3)

    @test metadata(df, "Description") == "Accesses to the ground facilities."
    @test colmetadata(df, :duration, "Unit") == :s

    # == Duration Units ====================================================================

    # -- Minutes ---------------------------------------------------------------------------

    df = ground_facility_accesses(orbp, (0, 0, 0) , 1 * 86400, TOD(), PEF(); unit = :m)

    @test size(df) == (2, 3)

    exp_duration = getfield.(df.access_end .- df.access_beginning, :value) ./ 1000 / 60

    @test df.duration[1] ≈ exp_duration[1]
    @test df.duration[2] ≈ exp_duration[2]

    @test metadata(df, "Description") == "Accesses to the ground facilities."
    @test colmetadata(df, :duration, "Unit") == :m

    # -- Hours -----------------------------------------------------------------------------

    df = ground_facility_accesses(orbp, (0, 0, 0) , 1 * 86400, TOD(), PEF(); unit = :h)

    @test size(df) == (2, 3)

    exp_duration = getfield.(df.access_end .- df.access_beginning, :value) ./ 1000 / 3600

    @test df.duration[1] ≈ exp_duration[1]
    @test df.duration[2] ≈ exp_duration[2]

    @test metadata(df, "Description") == "Accesses to the ground facilities."
    @test colmetadata(df, :duration, "Unit") == :h

    # -- Unknown Symbol --------------------------------------------------------------------

    df = ground_facility_accesses(orbp, (0, 0, 0) , 1 * 86400, TOD(), PEF(); unit = :not_know)

    @test size(df) == (2, 3)

    exp_duration = getfield.(df.access_end .- df.access_beginning, :value) ./ 1000

    @test df.duration[1] ≈ exp_duration[1]
    @test df.duration[2] ≈ exp_duration[2]

    @test metadata(df, "Description") == "Accesses to the ground facilities."
    @test colmetadata(df, :duration, "Unit") == :s
end

@testset "Function ground_facility_gaps" begin
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
    orbp = Propagators.init(Val(:J2), orb)

    # == Single Facility ===================================================================

    df = ground_facility_gaps(orbp, (0, 0, 0) , 1 * 86400, TOD(), PEF())

    @test size(df) == (3, 3)

    @test df.gap_beginning[1] == DateTime("2021-01-01T00:00:00.000")
    @test df.gap_beginning[2] == DateTime("2021-01-01T10:30:02.985")
    @test df.gap_beginning[3] == DateTime("2021-01-01T22:59:23.524")

    @test df.gap_end[1] == DateTime("2021-01-01T10:20:03.163")
    @test df.gap_end[2] == DateTime("2021-01-01T22:49:55.898")
    @test df.gap_end[3] == DateTime("2021-01-02T00:00:00.000")

    exp_duration = getfield.(df.gap_end .- df.gap_beginning, :value) ./ 1000

    @test df.duration[1] ≈ exp_duration[1]
    @test df.duration[2] ≈ exp_duration[2]
    @test df.duration[3] ≈ exp_duration[3]

    @test metadata(df, "Description") == "Gaps to the ground facilities."
    @test colmetadata(df, :duration, "Unit") == :s

    # == Multiple Facilities ===============================================================

    df = ground_facility_gaps(
        orbp,
        [(0, 0, 0), (0, π / 4, 0)],
        1 * 86400,
        TOD(),
        PEF()
    )

    @test size(df) == (5, 3)

    @test df.gap_beginning[1] == DateTime("2021-01-01T00:00:00.000")
    @test df.gap_beginning[2] == DateTime("2021-01-01T07:10:04.095")
    @test df.gap_beginning[3] == DateTime("2021-01-01T10:30:02.985")
    @test df.gap_beginning[4] == DateTime("2021-01-01T19:39:56.981")
    @test df.gap_beginning[5] == DateTime("2021-01-01T22:59:23.524")

    @test df.gap_end[1] == DateTime("2021-01-01T07:00:36.953")
    @test df.gap_end[2] == DateTime("2021-01-01T10:20:03.163")
    @test df.gap_end[3] == DateTime("2021-01-01T19:29:57.333")
    @test df.gap_end[4] == DateTime("2021-01-01T22:49:55.898")
    @test df.gap_end[5] == DateTime("2021-01-02T00:00:00.000")

    exp_duration = getfield.(df.gap_end .- df.gap_beginning, :value) ./ 1000

    @test df.duration[1] ≈ exp_duration[1]
    @test df.duration[2] ≈ exp_duration[2]
    @test df.duration[3] ≈ exp_duration[3]
    @test df.duration[4] ≈ exp_duration[4]
    @test df.duration[5] ≈ exp_duration[5]

    @test metadata(df, "Description") == "Gaps to the ground facilities."
    @test colmetadata(df, :duration, "Unit") == :s

    # If we change the reduction function to AND instead of OR, we should have not access.
    # Hence, the gap must be equal to the interval of analysis.
    df = ground_facility_gaps(
        orbp,
        [(0, 0, 0), (0, π / 4, 0)],
        1 * 86400,
        TOD(),
        PEF();
        reduction = v -> (&)(v...)
    )

    @test size(df) == (1, 3)

    @test df.gap_beginning[1] == DateTime("2021-01-01T00:00:00.000")

    @test df.gap_end[1] == DateTime("2021-01-02T00:00:00.000")

    @test df.duration[1] ≈ 86400.0

    @test metadata(df, "Description") == "Gaps to the ground facilities."
    @test colmetadata(df, :duration, "Unit") == :s

    # == Duration Units ====================================================================

    # -- Minutes ---------------------------------------------------------------------------

    df = ground_facility_gaps(orbp, (0, 0, 0) , 1 * 86400, TOD(), PEF(); unit = :m)

    @test size(df) == (3, 3)

    @test df.gap_beginning[1] == DateTime("2021-01-01T00:00:00.000")
    @test df.gap_beginning[2] == DateTime("2021-01-01T10:30:02.985")
    @test df.gap_beginning[3] == DateTime("2021-01-01T22:59:23.524")

    @test df.gap_end[1] == DateTime("2021-01-01T10:20:03.163")
    @test df.gap_end[2] == DateTime("2021-01-01T22:49:55.898")
    @test df.gap_end[3] == DateTime("2021-01-02T00:00:00.000")

    exp_duration = getfield.(df.gap_end .- df.gap_beginning, :value) ./ 1000 / 60

    @test df.duration[1] ≈ exp_duration[1]
    @test df.duration[2] ≈ exp_duration[2]
    @test df.duration[3] ≈ exp_duration[3]

    @test metadata(df, "Description") == "Gaps to the ground facilities."
    @test colmetadata(df, :duration, "Unit") == :m

    # -- Hours -----------------------------------------------------------------------------

    df = ground_facility_gaps(orbp, (0, 0, 0) , 1 * 86400, TOD(), PEF(); unit = :h)

    @test size(df) == (3, 3)

    @test df.gap_beginning[1] == DateTime("2021-01-01T00:00:00.000")
    @test df.gap_beginning[2] == DateTime("2021-01-01T10:30:02.985")
    @test df.gap_beginning[3] == DateTime("2021-01-01T22:59:23.524")

    @test df.gap_end[1] == DateTime("2021-01-01T10:20:03.163")
    @test df.gap_end[2] == DateTime("2021-01-01T22:49:55.898")
    @test df.gap_end[3] == DateTime("2021-01-02T00:00:00.000")

    exp_duration = getfield.(df.gap_end .- df.gap_beginning, :value) ./ 1000 / 3600

    @test df.duration[1] ≈ exp_duration[1]
    @test df.duration[2] ≈ exp_duration[2]
    @test df.duration[3] ≈ exp_duration[3]

    @test metadata(df, "Description") == "Gaps to the ground facilities."
    @test colmetadata(df, :duration, "Unit") == :h

    # -- Unknown Symbol --------------------------------------------------------------------

    df = ground_facility_gaps(orbp, (0, 0, 0) , 1 * 86400, TOD(), PEF(); unit = :not_known)

    @test size(df) == (3, 3)

    @test df.gap_beginning[1] == DateTime("2021-01-01T00:00:00.000")
    @test df.gap_beginning[2] == DateTime("2021-01-01T10:30:02.985")
    @test df.gap_beginning[3] == DateTime("2021-01-01T22:59:23.524")

    @test df.gap_end[1] == DateTime("2021-01-01T10:20:03.163")
    @test df.gap_end[2] == DateTime("2021-01-01T22:49:55.898")
    @test df.gap_end[3] == DateTime("2021-01-02T00:00:00.000")

    exp_duration = getfield.(df.gap_end .- df.gap_beginning, :value) ./ 1000

    @test df.duration[1] ≈ exp_duration[1]
    @test df.duration[2] ≈ exp_duration[2]
    @test df.duration[3] ≈ exp_duration[3]

    @test metadata(df, "Description") == "Gaps to the ground facilities."
    @test colmetadata(df, :duration, "Unit") == :s
end
