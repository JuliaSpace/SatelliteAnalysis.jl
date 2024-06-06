## Description #############################################################################
#
# Tests related to issues.
#
############################################################################################

@testset "Issue #3 - Gap function is not using the correct step" begin
    epoch = DateTime("2024-01-01")
    gs = (0.0, 0.0, 0.0)
    orb = KeplerianElements(
        epoch |> datetime2julian,
        7878.14e3,
        0.0,
        deg2rad(60),
        0.0,
        0.0,
        0.0
    )

    orbp = Propagators.init(Val(:J2), orb)

    acc = ground_facility_accesses(
        orbp,
        gs;
        minimum_elevation = deg2rad(60),
        duration = 24 * 3600,
        step = 120
    )

    gap = ground_facility_gaps(
        orbp,
        gs;
        minimum_elevation = deg2rad(60),
        duration = 24 * 3600,
        step = 120
    )

    @test nrow(acc) == 1
    @test nrow(gap) == 2
    @test acc[1, 1] == gap[1, 2]
    @test acc[1, 2] == gap[2, 1]
end
