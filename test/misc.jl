## Description #############################################################################
#
# Tests related to the miscellaneous functions.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# == File: ./src/misc/miscellaneous.jl =====================================================

# -- Function: find_crossing ---------------------------------------------------------------

@testset "Function find_crossing" begin
    f(t, a) = sin(t - a) ≥ 0

    for a in (-0.5, 0, 0.5)
        t = SatelliteAnalysis.find_crossing(f, a - 0.1, a + 0.1, false, true, a; Δ = 1e-20)

        @test t ≈ a
    end
end

@testset "Function find_crossing [ERRORS]" begin
    f(t, a) = :test

    @test_throws Exception SatelliteAnalysis.find_crossing(
        f, -0.1, +0.1, false, true, 10; Δ = 1e-20
    )
end

# == File: ./src/misc/units.jl =============================================================

@testset "Unit Factors" begin
    @test SatelliteAnalysis._angle_unit_factor(:rad) == 1
    @test SatelliteAnalysis._angle_unit_factor(:deg) ≈ 180 / π
    @test_throws ArgumentError SatelliteAnalysis._angle_unit_factor(:unknown)

    @test SatelliteAnalysis._distance_unit_factor(:m) == 1
    @test SatelliteAnalysis._distance_unit_factor(:km) ≈ 1e-3
    @test_throws ArgumentError SatelliteAnalysis._distance_unit_factor(:unknown)

    @test SatelliteAnalysis._time_unit_factor(:s) == 1
    @test SatelliteAnalysis._time_unit_factor(:m) ≈ 1 / 60
    @test SatelliteAnalysis._time_unit_factor(:h) ≈ 1 / 3600
    @test SatelliteAnalysis._time_unit_factor(:d, (:s, :d, :y)) ≈ 1 / 86400
    @test SatelliteAnalysis._time_unit_factor(:y, (:s, :d, :y)) ≈ 1 / (365.25 * 86400)

    # The days and years are not valid by default.
    @test_throws ArgumentError SatelliteAnalysis._time_unit_factor(:d)
    @test_throws ArgumentError SatelliteAnalysis._time_unit_factor(:y)
    @test_throws ArgumentError SatelliteAnalysis._time_unit_factor(:unknown)
    @test_throws ArgumentError SatelliteAnalysis._time_unit_factor(:m, (:s, :d, :y))
end
