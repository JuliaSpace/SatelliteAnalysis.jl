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
    num_sampling_points_per_orbit::Int = 33,
    Ap::Union{Nothing, Number} = nothing,
    C_d::Number = 2.2,
    C_r::Number = 1.25,
    F107::Union{Nothing, Number} = nothing,
    terminate_altitude::Number = 120e3,
    tf::Number = 30 * 365.25 * 86400.0,
)
    M = true_to_mean_anomaly(orb.e, orb.f)
    ā, ē, ī, Ω̄, ω̄, M̄ = _osculating_to_mean_elements(orb.a, orb.e, orb.i, orb.Ω, orb.ω, M)
    u = Vector(_classical_to_equinoctial(ā, ē, ī, Ω̄, ω̄, M̄))
    # u = Vector(_classical_to_equinoctial(orb.a, orb.e, orb.i, orb.Ω, orb.ω, M))

    tspan = (0.0, tf)

    gm = isnothing(gravity_model) ?
        GravityModels.load(IcgemFile, fetch_icgem_file(:EGM2008)) :
        gravity_model

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
        terminate_altitude            = terminate_altitude,
    )

    cbset = CallbackSet(ContinuousCallback(_cb_altitude_condition, _cb_altitude_affect!))

    prob = ODEProblem(_dynamics!, u, tspan, params)
    sol = solve(
        prob,
        Tsit5();
        dt       = 24.0 * 60 * 60,
        reltol   = 1e-8,
        abstol   = 1e-8,
        callback = cbset,
        maxiters = 1e8
    )

    return sol
end

############################################################################################
#                                    Private Functions                                     #
############################################################################################

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
