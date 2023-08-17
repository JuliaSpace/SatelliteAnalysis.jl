# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#   Tests related to the eclipse time analysis.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# File: ./src/eclipse_time.jl
# ==========================================================================================

# Functions: eclipse_time_summary
# ------------------------------------------------------------------------------------------

@testset "Function eclipse_time_summary" verbose = true begin
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

    # Seconds
    # ======================================================================================

    df = eclipse_time_summary(orbp; num_days = 5)

    @test size(df) == (5, 4)

    @test sum(df.sunlight) / 5 ≈ 3974.7845507883844
    @test sum(df.penumbra) / 5 ≈ 20.45694465948482
    @test sum(df.umbra)    / 5 ≈ 2004.756976375114

    @test metadata(df, "Description") == "Eclipse time PER ORBIT computed at each day."

    @test colmetadata(df, :sunlight, "Unit") == :s
    @test colmetadata(df, :penumbra, "Unit") == :s
    @test colmetadata(df, :umbra,    "Unit") == :s

    # Minutes
    # ======================================================================================

    df = eclipse_time_summary(orbp; num_days = 5, unit = :m)

    @test size(df) == (5, 4)

    @test sum(df.sunlight) / 5 ≈ 3974.7845507883844 / 60
    @test sum(df.penumbra) / 5 ≈ 20.45694465948482  / 60
    @test sum(df.umbra)    / 5 ≈ 2004.756976375114  / 60

    @test metadata(df, "Description") == "Eclipse time PER ORBIT computed at each day."

    @test colmetadata(df, :sunlight, "Unit") == :m
    @test colmetadata(df, :penumbra, "Unit") == :m
    @test colmetadata(df, :umbra,    "Unit") == :m

    # Hours
    # ======================================================================================

    df = eclipse_time_summary(orbp; num_days = 5, unit = :h)

    @test size(df) == (5, 4)

    @test sum(df.sunlight) / 5 ≈ 3974.7845507883844 / 3600
    @test sum(df.penumbra) / 5 ≈ 20.45694465948482  / 3600
    @test sum(df.umbra)    / 5 ≈ 2004.756976375114  / 3600

    @test metadata(df, "Description") == "Eclipse time PER ORBIT computed at each day."

    @test colmetadata(df, :sunlight, "Unit") == :h
    @test colmetadata(df, :penumbra, "Unit") == :h
    @test colmetadata(df, :umbra,    "Unit") == :h

    # Unknown Symbol
    # ======================================================================================

    df = eclipse_time_summary(orbp; num_days = 5, unit = :not_known)

    @test size(df) == (5, 4)

    @test sum(df.sunlight) / 5 ≈ 3974.7845507883844
    @test sum(df.penumbra) / 5 ≈ 20.45694465948482
    @test sum(df.umbra)    / 5 ≈ 2004.756976375114

    @test metadata(df, "Description") == "Eclipse time PER ORBIT computed at each day."

    @test colmetadata(df, :sunlight, "Unit") == :s
    @test colmetadata(df, :penumbra, "Unit") == :s
    @test colmetadata(df, :umbra,    "Unit") == :s
end
