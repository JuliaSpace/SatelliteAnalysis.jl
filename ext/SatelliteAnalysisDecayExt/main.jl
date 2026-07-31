## Description #############################################################################
#
# Main function to compute the decay analysis of a satellite in orbit.
#
############################################################################################

function SatelliteAnalysis.decay_analysis(
    orb::KeplerianElements;
    # Required keywords.
    satellite_mass::Number,
    satellite_mean_area::Number,
    # Optional keywords.
    gravity_model::Union{AbstractGravityModel, Nothing} = nothing,
    num_sampling_points_per_orbit::Int = 17,
    abstol::Number = 1e-6,
    Ap::Union{Nothing, Number} = nothing,
    C_d::Number = 2.2,
    C_r::Number = 1.25,
    distance_unit::Symbol = :km,
    F107::Union{Nothing, Number} = nothing,
    reltol::Number = 1e-6,
    return_solution::Bool = false,
    solver = VCABM(),
    terminate_altitude::Number = 120e3,
    tf::Number = 30 * 365.25 * 86400.0,
    time_unit::Symbol = :y,
)
    gm = isnothing(gravity_model) ?
        GravityModels.load(IcgemFile, fetch_icgem_file(:EGM2008)) :
        gravity_model

    # The keyword `gravity_model` is abstractly typed. This function barrier ensures the
    # gravity model is concretely typed inside the numerical integration, avoiding dynamic
    # dispatch at every right-hand-side evaluation.
    return _decay_analysis(
        orb,
        gm;
        satellite_mass                = satellite_mass,
        satellite_mean_area           = satellite_mean_area,
        num_sampling_points_per_orbit = num_sampling_points_per_orbit,
        abstol                        = abstol,
        Ap                            = Ap,
        C_d                           = C_d,
        C_r                           = C_r,
        distance_unit                 = distance_unit,
        F107                          = F107,
        reltol                        = reltol,
        return_solution               = return_solution,
        solver                        = solver,
        terminate_altitude            = terminate_altitude,
        tf                            = tf,
        time_unit                     = time_unit
    )
end

############################################################################################
#                                    Private Functions                                     #
############################################################################################

function _decay_analysis(
    orb::KeplerianElements,
    gm::AbstractGravityModel;
    satellite_mass::Number,
    satellite_mean_area::Number,
    num_sampling_points_per_orbit::Int,
    abstol::Number,
    Ap::Union{Nothing, Number},
    C_d::Number,
    C_r::Number,
    distance_unit::Symbol,
    F107::Union{Nothing, Number},
    reltol::Number,
    return_solution::Bool,
    solver,
    terminate_altitude::Number,
    tf::Number,
    time_unit::Symbol,
)
    M = true_to_mean_anomaly(orb.e, orb.f)
    ā, ē, ī, Ω̄, ω̄, M̄ = _osculating_to_mean_elements(orb.a, orb.e, orb.i, orb.Ω, orb.ω, M)
    u = Vector(_classical_to_equinoctial(ā, ē, ī, Ω̄, ω̄, M̄))

    tspan = (0.0, tf)

    # Pre-allocate the J2 osculating propagator used for the mean-to-osculating
    # conversion inside the drag quadrature, avoiding one propagator allocation per
    # sampling point.
    j2osc_prop = Propagators.init(
        Val(:J2osc),
        KeplerianElements(0.0, ā, ē, ī, Ω̄, ω̄, mean_to_true_anomaly(ē, M̄))
    )

    params = (
        satellite_mean_area           = satellite_mean_area,
        satellite_mass                = satellite_mass,
        num_sampling_points_per_orbit = num_sampling_points_per_orbit,
        Ap                            = Ap,
        C_d                           = C_d,
        C_r                           = C_r,
        F107                          = F107,
        gm                            = gm,
        jd₀_utc                       = orb.t,
        j2osc_prop                    = j2osc_prop,
        terminate_altitude            = terminate_altitude,
    )

    cbset = CallbackSet(ContinuousCallback(_cb_altitude_condition, _cb_altitude_affect!))

    prob = ODEProblem(_dynamics!, u, tspan, params)
    sol = solve(
        prob,
        solver;
        dt       = 24.0 * 60 * 60,
        reltol   = reltol,
        abstol   = abstol,
        callback = cbset,
        maxiters = 1e8
    )

    # == Assemble the Output ===============================================================

    num_points = length(sol.t)

    date             = Vector{DateTime}(undef, num_points)
    time             = Vector{Float64}(undef, num_points)
    f107             = Vector{Float64}(undef, num_points)
    ap               = Vector{Float64}(undef, num_points)
    mean_elements    = Vector{KeplerianElements{Float64, Float64}}(undef, num_points)
    apogee_altitude  = Vector{Float64}(undef, num_points)
    perigee_altitude = Vector{Float64}(undef, num_points)

    @inbounds for k in 1:num_points
        aₖ, eₖ, iₖ, Ωₖ, ωₖ, Mₖ = _equinoctial_to_classical(sol.u[k])

        jdₖ = orb.t + sol.t[k] / 86400

        date[k] = julian2datetime(jdₖ)
        time[k] = sol.t[k]

        # Record the space indices used by the dynamics at each instant.
        f107[k] = isnothing(F107) ?
            Float64(space_index(Val(:F10obs_avg_center81), jdₖ)) :
            Float64(F107)
        ap[k] = isnothing(Ap) ? Float64(space_index(Val(:Ap_daily), jdₖ)) : Float64(Ap)

        mean_elements[k] = KeplerianElements(
            jdₖ, aₖ, eₖ, iₖ, Ωₖ, ωₖ, mean_to_true_anomaly(eₖ, Mₖ)
        )

        apogee_altitude[k]  = aₖ * (1 + eₖ) - EARTH_EQUATORIAL_RADIUS
        perigee_altitude[k] = aₖ * (1 - eₖ) - EARTH_EQUATORIAL_RADIUS
    end

    # Convert the time and altitude columns to the selected units.
    if time_unit == :m
        time ./= 60
    elseif time_unit == :h
        time ./= 3600
    elseif time_unit == :d
        time ./= 86400
    elseif time_unit != :s
        # Julian year, consistent with the default `tf` of 30 years. If the symbol is not
        # known, we must use the default unit (years).
        time ./= 365.25 * 86400
        time_unit = :y
    end

    if distance_unit != :m
        # If the symbol is not known, we must use the default unit (kilometers).
        apogee_altitude  ./= 1000
        perigee_altitude ./= 1000
        distance_unit = :km
    end

    df = DataFrame(;
        date             = date,
        time             = time,
        f107             = f107,
        ap               = ap,
        mean_elements    = mean_elements,
        apogee_altitude  = apogee_altitude,
        perigee_altitude = perigee_altitude,
    )

    metadata!(df, "Description", "Mean orbital element evolution during the orbital decay.")
    colmetadata!(df, :date,             "Unit", :UTC)
    colmetadata!(df, :time,             "Unit", time_unit)
    colmetadata!(df, :f107,             "Unit", :sfu)
    colmetadata!(df, :ap,               "Unit", :dimensionless)
    colmetadata!(df, :mean_elements,    "Unit", :SI)
    colmetadata!(df, :apogee_altitude,  "Unit", distance_unit)
    colmetadata!(df, :perigee_altitude, "Unit", distance_unit)

    return_solution && return df, sol

    return df
end

function _cb_altitude_condition(u, t, integrator)
    a_mean, e_mean, i_mean, Ω_mean, ω_mean, M_mean = _equinoctial_to_classical(u)

    # Terminate on the perigee altitude of the mean orbit. This condition is a smooth,
    # monotonically decreasing function of the mean elements, whereas the instantaneous
    # altitude oscillates between the perigee and apogee altitudes within a single
    # integrator step, which can make the callback root finder miss the crossing.
    perigee_altitude = a_mean * (1 - e_mean) - EARTH_EQUATORIAL_RADIUS

    return perigee_altitude - integrator.p.terminate_altitude
end

function _cb_altitude_affect!(integrator)
    terminate!(integrator)
    return nothing
end
