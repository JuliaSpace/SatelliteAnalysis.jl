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
        "time",
        "f107",
        "ap",
        "mean_elements",
        "apogee_altitude",
        "perigee_altitude",
    ]

    @test eltype(df.date)             == DateTime
    @test eltype(df.time)             == Float64
    @test eltype(df.f107)             == Float64
    @test eltype(df.ap)               == Float64
    @test eltype(df.mean_elements)    == KeplerianElements{Float64, Float64}
    @test eltype(df.apogee_altitude)  == Float64
    @test eltype(df.perigee_altitude) == Float64

    @test metadata(df, "Description") ==
        "Mean orbital element evolution during the orbital decay."

    @test colmetadata(df, :date,             "Unit") == :UTC
    @test colmetadata(df, :time,             "Unit") == :y
    @test colmetadata(df, :f107,             "Unit") == :sfu
    @test colmetadata(df, :ap,               "Unit") == :dimensionless
    @test colmetadata(df, :mean_elements,    "Unit") == :SI
    @test colmetadata(df, :apogee_altitude,  "Unit") == :km
    @test colmetadata(df, :perigee_altitude, "Unit") == :km

    # == Values ============================================================================

    @test sol.retcode == ReturnCode.Terminated

    # The analysis must start at the orbit epoch.
    @test df[begin, :date] == julian2datetime(jd₀)
    @test df[begin, :time] == 0.0
    @test df[end,   :time] ≈ (datetime2julian(df[end, :date]) - jd₀) / 365.25 atol = 1e-8

    # The space indices provided by the user must be recorded in the output.
    @test all(df.f107 .== 140.0)
    @test all(df.ap   .== 15.0)

    # The mean element epochs must match the point dates, and the apogee must be above the
    # perigee.
    @test julian2datetime(df[end, :mean_elements].t) == df[end, :date]
    @test all(df.apogee_altitude .>= df.perigee_altitude)

    # The termination must happen when the mean perigee altitude reaches 120 km.
    @test df[end, :perigee_altitude] ≈ 120.0 atol = 1e-6

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
    @test df_200[end, :perigee_altitude] ≈ 200.0 atol = 1e-6
    @test df_200[end, :date] < df[end, :date]

    # == Keywords time_unit and distance_unit =============================================

    df_u = decay_analysis(
        orb;
        satellite_mass      = 100.0,
        satellite_mean_area = 1.0,
        gravity_model       = gm,
        F107                = 140.0,
        Ap                  = 15.0,
        distance_unit       = :m,
        time_unit           = :s
    )

    @test df_u.time             ≈ df.time .* (365.25 * 86400)
    @test df_u.apogee_altitude  ≈ df.apogee_altitude .* 1000
    @test df_u.perigee_altitude ≈ df.perigee_altitude .* 1000

    @test colmetadata(df_u, :time,             "Unit") == :s
    @test colmetadata(df_u, :apogee_altitude,  "Unit") == :m
    @test colmetadata(df_u, :perigee_altitude, "Unit") == :m

    df_d = decay_analysis(
        orb;
        satellite_mass      = 100.0,
        satellite_mean_area = 1.0,
        gravity_model       = gm,
        F107                = 140.0,
        Ap                  = 15.0,
        time_unit           = :d
    )

    @test df_d.time ≈ df.time .* 365.25
    @test colmetadata(df_d, :time, "Unit") == :d

    # Unknown unit symbols must fall back to the defaults (years and kilometers).
    df_f = decay_analysis(
        orb;
        satellite_mass      = 100.0,
        satellite_mean_area = 1.0,
        gravity_model       = gm,
        F107                = 140.0,
        Ap                  = 15.0,
        distance_unit       = :unknown,
        time_unit           = :unknown
    )

    @test df_f.time             ≈ df.time
    @test df_f.perigee_altitude ≈ df.perigee_altitude
    @test colmetadata(df_f, :time,             "Unit") == :y
    @test colmetadata(df_f, :apogee_altitude,  "Unit") == :km
    @test colmetadata(df_f, :perigee_altitude, "Unit") == :km
end

@testset "Validation Against a Cowell Reference" begin
    # Propagate the full osculating dynamics (Cowell formulation) using the same force
    # model implementations of the extension, and compare the mean elements after seven
    # days with the averaged model. The tolerances below quantify the couplings neglected
    # by the averaged formulation (e.g. short-period drag-density coupling and J₂ × J₃
    # terms). They were obtained from a reference run with a 2x margin.
    ext = Base.get_extension(SatelliteAnalysis, :SatelliteAnalysisDecayExt)

    duration = 7 * 86400.0
    F107     = 140.0
    Ap       = 15.0
    mass     = 100.0
    area     = 1.0
    C_d      = 2.2
    C_r      = 1.25

    jd₀ = date_to_jd(2024, 1, 1)

    orb = KeplerianElements(
        jd₀,
        EARTH_EQUATORIAL_RADIUS + 500e3,
        0.001,
        98.0 |> deg2rad,
        ltdn_to_raan(10.5, jd₀),
        90.0 |> deg2rad,
        0.0
    )

    gm = GravityModels.load(IcgemFile, fetch_icgem_file(:EGM2008))

    # == Cowell Reference ==================================================================

    function cowell_dynamics!(du, u, params, t)
        jd_utc = params.jd₀_utc + t / 86400

        r_tod = u[1:3]
        v_tod = u[4:6]

        # Frames and third-body positions, following the same procedure of the averaged
        # dynamics.
        rsun_mod  = sun_position_mod(jd_utc)
        rmoon_mod = moon_position_mod(jd_utc, Val(:Vallado))
        D_tod_mod = r_eci_to_eci(MOD(), jd_utc, TOD(), jd_utc)
        rsun_tod  = D_tod_mod * rsun_mod
        rmoon_tod = D_tod_mod * rmoon_mod

        D_pef_tod = r_eci_to_ecef(TOD(), PEF(), jd_utc)
        D_tod_pef = D_pef_tod'
        ω_pef     = @SVector [0.0, 0.0, EARTH_ANGULAR_SPEED]

        r_pef = D_pef_tod * r_tod
        v_pef = D_pef_tod * v_tod - ω_pef × r_pef

        # Total gravity: central term plus zonal harmonics up to degree 7, as in the
        # averaged model.
        a_grav_pef = GravityModels.gravitational_acceleration(
            params.gm, r_pef, jd_utc; max_degree = 7, max_order = 0
        )

        # Atmospheric drag using the same routine and space indices.
        a_drag_pef = ext._atmospheric_drag_acceleration(
            jd_utc, r_pef, v_pef, area, mass, C_d; F107 = F107, Ap = Ap
        )

        # Third-body point masses using the same routine.
        a_3b_tod = ext._point_mass_acceleration(r_tod, rsun_tod, ext._μ_SUN) +
            ext._point_mass_acceleration(r_tod, rmoon_tod, ext._μ_MOON)

        # Solar radiation pressure with the same shadow rule.
        lc = lighting_condition(r_tod, rsun_tod)
        ν  = lc == :sunlight ? 1.0 : (lc == :penumbra ? 0.5 : 0.0)
        a_srp_tod = ν * ext._solar_radiation_acceleration(r_tod, rsun_tod, area, mass, C_r)

        a_tod = D_tod_pef * (a_grav_pef + a_drag_pef) + a_3b_tod + a_srp_tod

        du[1:3] .= v_tod
        du[4:6] .= a_tod

        return nothing
    end

    r₀, v₀ = kepler_to_rv(orb)
    u₀     = vcat(Vector(r₀), Vector(v₀))
    prob   = ODEProblem(cowell_dynamics!, u₀, (0.0, duration), (gm = gm, jd₀_utc = jd₀))
    sol    = solve(prob, Vern7(); reltol = 1e-11, abstol = 1e-6)

    @test sol.retcode == ReturnCode.Success

    # Convert the final osculating state to mean elements using the same convention of the
    # extension.
    rf = sol.u[end][1:3]
    vf = sol.u[end][4:6]
    ke_mean, ~ = fit_j2osc_mean_elements(
        [0.0], [rf], [vf]; max_iterations = 50, verbose = false
    )
    M_mean = true_to_mean_anomaly(ke_mean.e, ke_mean.f)

    # == Averaged Model ====================================================================

    df = decay_analysis(
        orb;
        satellite_mass      = mass,
        satellite_mean_area = area,
        gravity_model       = gm,
        C_d                 = C_d,
        C_r                 = C_r,
        F107                = F107,
        Ap                  = Ap,
        tf                  = duration
    )

    # == Comparison ========================================================================

    wrap(x) = mod(x + π, 2π) - π

    ke_beg = df[begin, :mean_elements]
    ke_end = df[end,   :mean_elements]
    M_end  = true_to_mean_anomaly(ke_end.e, ke_end.f)

    # The averaged model must show significant motion, otherwise the comparison is
    # meaningless.
    @test ke_end.a - ke_beg.a < -200
    @test abs(wrap(ke_end.Ω - ke_beg.Ω)) > 0.1

    @test abs(ke_end.a - ke_mean.a)          < 150
    @test abs(ke_end.e - ke_mean.e)          < 2e-5
    @test abs(wrap(ke_end.i - ke_mean.i))    < 2e-6
    @test abs(wrap(ke_end.Ω - ke_mean.Ω))    < 4e-4
    @test abs(wrap(ke_end.ω - ke_mean.ω))    < 2e-4
    @test abs(wrap(M_end - M_mean))          < 2e-2
end
