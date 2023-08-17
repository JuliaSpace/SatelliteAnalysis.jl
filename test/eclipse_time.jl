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
    orbp = Propagators.init(Val(:J2), orb)

    @testset "Without IO" begin
        s, p, u = eclipse_time_summary(orbp; num_days = 5)

        @test length(s) == 5
        @test length(p) == 5
        @test length(u) == 5

        @test sum(s) / 5 ≈ 3974.7845507883844
        @test sum(p) / 5 ≈ 20.45694465948482
        @test sum(u) / 5 ≈ 2004.756976375114
    end

    @testset "With IO" begin
        expected = """
                             Day │ Sunlight time [s]  Penumbra time [s]  Umbra time [s]
            ─────────────────────┼──────────────────────────────────────────────────────
             2021-01-01T00:00:00 │           3972.63            20.4117         2006.96
             2021-01-02T00:00:00 │           3973.85            20.4376         2005.71
             2021-01-03T00:00:00 │           3974.77            20.4575         2004.77
             2021-01-04T00:00:00 │           3975.74            20.4758         2003.79
             2021-01-05T00:00:00 │           3976.94            20.5022         2002.55
            ─────────────────────┼──────────────────────────────────────────────────────
                         Maximum │           3976.94            20.5022         2006.96
                            Mean │           3974.78            20.4569         2004.76
                         Minimum │           3972.63            20.4117         2002.55
            """

        buf = IOBuffer()
        eclipse_time_summary(buf, orbp; num_days = 5)
        result = String(take!(buf))
        @test result == expected

        expected = """
                             Day │ Sunlight time [min]  Penumbra time [min]  Umbra time [min]
            ─────────────────────┼────────────────────────────────────────────────────────────
             2021-01-01T00:00:00 │             66.2105             0.340195           33.4493
             2021-01-02T00:00:00 │             66.2308             0.340627           33.4285
             2021-01-03T00:00:00 │             66.2461             0.340958           33.4129
             2021-01-04T00:00:00 │             66.2623             0.341263           33.3964
             2021-01-05T00:00:00 │             66.2824             0.341704           33.3759
            ─────────────────────┼────────────────────────────────────────────────────────────
                         Maximum │             66.2824             0.341704           33.4493
                            Mean │             66.2464             0.340949           33.4126
                         Minimum │             66.2105             0.340195           33.3759
            """

        buf = IOBuffer()
        eclipse_time_summary(buf, orbp; num_days = 5, unit = :m)
        result = String(take!(buf))
        @test result == expected

        expected = """
                             Day │ Sunlight time [h]  Penumbra time [h]  Umbra time [h]
            ─────────────────────┼──────────────────────────────────────────────────────
             2021-01-01T00:00:00 │           1.10351         0.00566991        0.557489
             2021-01-02T00:00:00 │           1.10385         0.00567712        0.557142
             2021-01-03T00:00:00 │           1.1041          0.00568263        0.556882
             2021-01-04T00:00:00 │           1.10437         0.00568771        0.556607
             2021-01-05T00:00:00 │           1.10471         0.00569506        0.556265
            ─────────────────────┼──────────────────────────────────────────────────────
                         Maximum │           1.10471         0.00569506        0.557489
                            Mean │           1.10411         0.00568248        0.556877
                         Minimum │           1.10351         0.00566991        0.556265
            """

        buf = IOBuffer()
        eclipse_time_summary(buf, orbp; num_days = 5, unit = :h)
        result = String(take!(buf))
        @test result == expected

        # Warning
        # ==================================================================================

        buf = IOBuffer()
        @test_logs(
            (:warn, "We can only use the pager if printing to `stdout`."),
            eclipse_time_summary(buf, orbp; use_pager = true)
        )
    end
end
