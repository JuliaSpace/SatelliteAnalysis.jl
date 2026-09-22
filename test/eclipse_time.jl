## Description #############################################################################
#
# Tests related to the eclipse time analysis.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# == File: ./src/eclipse_time.jl ===========================================================

# -- Functions: eclipse_time_summary -------------------------------------------------------

@testset "Function eclipse_time_summary" verbose = true begin
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

    orbp = Propagators.init(Val(:J2), orb)

    # == Seconds ===========================================================================

    df = eclipse_time_summary(orbp; num_days = 5)

    @test size(df) == (5, 4)

    @test sum(df.sunlight) / 5 ≈ 3974.7845507883844
    @test sum(df.penumbra) / 5 ≈ 20.45694465948482
    @test sum(df.umbra) / 5 ≈ 2004.756976375114

    @test metadata(df, "Description") == "Eclipse time PER ORBIT computed at each day."

    @test colmetadata(df, :sunlight, "Unit") == :s
    @test colmetadata(df, :penumbra, "Unit") == :s
    @test colmetadata(df, :umbra, "Unit") == :s

    # All the metadata must use the style `:note` to propagate through DataFrame
    # transformations.
    @test all(k -> metadata(df, k; style = true)[2] == :note, metadatakeys(df))
    @test all(
        colmetadata(df, col, k; style = true)[2] == :note for
        (col, keys) in colmetadatakeys(df) for k in keys
    )

    # == Minutes ===========================================================================

    df = eclipse_time_summary(orbp; num_days = 5, time_unit = :min)

    @test size(df) == (5, 4)

    @test sum(df.sunlight) / 5 ≈ 3974.7845507883844 / 60
    @test sum(df.penumbra) / 5 ≈ 20.45694465948482 / 60
    @test sum(df.umbra) / 5 ≈ 2004.756976375114 / 60

    @test metadata(df, "Description") == "Eclipse time PER ORBIT computed at each day."

    @test colmetadata(df, :sunlight, "Unit") == :min
    @test colmetadata(df, :penumbra, "Unit") == :min
    @test colmetadata(df, :umbra, "Unit") == :min

    # == Hours =============================================================================

    df = eclipse_time_summary(orbp; num_days = 5, time_unit = :h)

    @test size(df) == (5, 4)

    @test sum(df.sunlight) / 5 ≈ 3974.7845507883844 / 3600
    @test sum(df.penumbra) / 5 ≈ 20.45694465948482 / 3600
    @test sum(df.umbra) / 5 ≈ 2004.756976375114 / 3600

    @test metadata(df, "Description") == "Eclipse time PER ORBIT computed at each day."

    @test colmetadata(df, :sunlight, "Unit") == :h
    @test colmetadata(df, :penumbra, "Unit") == :h
    @test colmetadata(df, :umbra, "Unit") == :h

    # == Regions Inside One Step ===========================================================

    # If the step is larger than the time the satellite stays in the penumbra, this region
    # can be entirely inside one step. In this case, the algorithm must find both edges.
    df_coarse = eclipse_time_summary(orbp; num_days = 5, step = 60)
    df_fine   = eclipse_time_summary(orbp; num_days = 5, step = 1)

    @test all(df_coarse.penumbra .> 15)
    @test df_coarse.sunlight ≈ df_fine.sunlight atol = 1e-2
    @test df_coarse.penumbra ≈ df_fine.penumbra atol = 1e-2
    @test df_coarse.umbra ≈ df_fine.umbra atol = 1e-2

    # == Invalid Inputs ====================================================================

    @test_throws ArgumentError eclipse_time_summary(orbp; num_days = 0)
    @test_throws ArgumentError eclipse_time_summary(orbp; num_days = 5, step = 0)
    @test_throws ArgumentError eclipse_time_summary(orbp; num_days = 5, step = -1)
    @test_throws ArgumentError eclipse_time_summary(orbp; num_days = 5, step = 7000)

    # == Unknown Symbol ====================================================================

    @test_throws ArgumentError eclipse_time_summary(orbp; num_days = 5, time_unit = :bad)
end
