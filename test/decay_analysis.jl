## Description #############################################################################
#
# Tests related to the orbital decay analysis.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# == File: ./src/decay_analysis.jl =========================================================

# -- Function: decay_analysis --------------------------------------------------------------

@testset "Extension Loading" begin
    ext = Base.get_extension(SatelliteAnalysis, :SatelliteAnalysisDecayExt)
    @test !isnothing(ext)
end

@testset "Function decay_analysis" begin
    jd₀ = date_to_jd(2024, 1, 1)

    orb = KeplerianElements(
        jd₀,
        EARTH_EQUATORIAL_RADIUS + 300e3,
        0.001,
        98.0    |> deg2rad,
        ltdn_to_raan(10.5, jd₀),
        90.0    |> deg2rad,
        0.0
    )

    gm = GravityModels.load(IcgemFile, fetch_icgem_file(:EGM2008))

    df, sol = decay_analysis(
        orb;
        satellite_mass      = 100.0,
        satellite_mean_area = 1.0,
        gravity_model       = gm,
        F107                = 140.0,
        Ap                  = 15.0,
        return_solution     = true
    )

    # == Schema and Metadata ===============================================================

    @test df isa DataFrame
    @test names(df) == [
        "date",
        "semi_major_axis",
        "eccentricity",
        "inclination",
        "raan",
        "argument_of_perigee",
        "mean_anomaly",
        "perigee_altitude",
    ]

    @test eltype(df.date)                == DateTime
    @test eltype(df.semi_major_axis)     == Float64
    @test eltype(df.eccentricity)        == Float64
    @test eltype(df.inclination)         == Float64
    @test eltype(df.raan)                == Float64
    @test eltype(df.argument_of_perigee) == Float64
    @test eltype(df.mean_anomaly)        == Float64
    @test eltype(df.perigee_altitude)    == Float64

    @test metadata(df, "Description") ==
        "Mean orbital element evolution during the orbital decay."

    @test colmetadata(df, :semi_major_axis,     "Unit") == :m
    @test colmetadata(df, :eccentricity,        "Unit") == :dimensionless
    @test colmetadata(df, :inclination,         "Unit") == :rad
    @test colmetadata(df, :raan,                "Unit") == :rad
    @test colmetadata(df, :argument_of_perigee, "Unit") == :rad
    @test colmetadata(df, :mean_anomaly,        "Unit") == :rad
    @test colmetadata(df, :perigee_altitude,    "Unit") == :m

    # == Values ============================================================================

    @test sol.retcode == ReturnCode.Terminated

    # The analysis must start at the orbit epoch.
    @test df[begin, :date] == julian2datetime(jd₀)

    # The termination must happen when the mean perigee altitude reaches 120 km.
    @test df[end, :perigee_altitude] ≈ 120e3 atol = 1e-3

    # Regression test for the estimated lifetime.
    lifetime = datetime2julian(df[end, :date]) - jd₀
    @test lifetime ≈ 31.328780 rtol = 1e-4

    # == Keyword terminate_altitude ========================================================

    df_200 = decay_analysis(
        orb;
        satellite_mass      = 100.0,
        satellite_mean_area = 1.0,
        gravity_model       = gm,
        F107                = 140.0,
        Ap                  = 15.0,
        terminate_altitude  = 200e3
    )

    @test df_200 isa DataFrame
    @test df_200[end, :perigee_altitude] ≈ 200e3 atol = 1e-3
    @test df_200[end, :date] < df[end, :date]
end
