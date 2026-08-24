## Description #############################################################################
#
# Dynamics function to compute the satellite decay.
#
############################################################################################

"""
    _dynamics(u::AbstractVector{T}, params::NamedTuple, t::Real) where T <: Number -> SVector{6, T}

Compute the time derivatives of the mean equinoctial orbital elements `u` at the time `t`
[s] after the epoch using Gauss variational equations with averaged perturbations.

This routine evaluates instantaneous forces on the mean-reference orbit for conservative
perturbations (Earth gravity and third bodies) and adds temporally averaged rates for
non-conservative perturbations (atmospheric drag and solar radiation pressure), together
with closed-form J₂² rates.

!!! note

    - Conservative perturbations are evaluated on the mean orbit at evenly spaced mean
      anomalies.
    - Non-conservative perturbations are added via averaged routines (drag and SRP).
    - Units must be consistent: meters, seconds, kilograms.

# Arguments

- `u::AbstractVector{T}`: Mean equinoctial orbital elements `[a, ψ, e_x, e_y, i_x, i_y]`.
- `params::NamedTuple`: Named tuple with the integration parameters. See
    `_decay_analysis`.
- `t::Real`: Time since the orbit epoch [s].

# Returns

- `SVector{6, T}`: Time derivatives of the mean equinoctial orbital elements
    `[da, dψ, de_x, de_y, di_x, di_y]`.

# References

- **[1]** Battin, R. H. (1999). An Introduction to the Mathematics and Methods of
    Astrodynamics. Revised ed. AIAA Education Series, Reston, VA.
"""
function _dynamics(u::AbstractVector{T}, params, t::Real) where T <: Number
    # Gravity model.
    gm = params.gm
    μ  = GravityModels.gravity_constant(gm)

    # Time and ephemerides.
    jd₀_utc = params.jd₀_utc
    jd_utc  = jd₀_utc + t / 86400.0

    # Unpack mean orbital elements.
    ā, ē, ī, Ω̄, ω̄, M̄ = _equinoctial_to_classical(u)

    # Clamp the mean elements to a physically meaningful region. The adaptive integrator
    # can evaluate trial stages with unphysical states (e.g. e > 1) near the decay end.
    # The right-hand side must remain finite for those states so that the error control
    # can reject the step instead of throwing a domain error.
    ā = max(ā, T(0.9) * T(EARTH_EQUATORIAL_RADIUS))
    ē = clamp(ē, T(1e-6), max(T(1e-6), 1 - T(EARTH_EQUATORIAL_RADIUS) / ā))
    ī = max(ī, T(1e-6))

    # Auxiliaries with mean elements.
    ē² = ē^2
    η² = 1 - ē²
    p̄  = ā * η²
    h̄  = √(μ * p̄)
    n̄  = √(μ / ā^3)
    η̄  = √η²

    # Third-body constants and positions.
    rsun_mod  = sun_position_mod(jd_utc)
    rmoon_mod = moon_position_mod(jd_utc, Val(:Vallado))

    # Convert the vectors to TOD.
    D_tod_mod = r_eci_to_eci(MOD(), jd_utc, TOD(), jd_utc)
    rsun_tod  = D_tod_mod * rsun_mod
    rmoon_tod = D_tod_mod * rmoon_mod

    # Convert the vectors to PEF. Both rotations are orthonormal DCMs, so the inverse
    # conversion is obtained by transposition instead of a second frame reduction.
    D_pef_tod = r_eci_to_ecef(TOD(), PEF(), jd_utc)
    D_tod_pef = D_pef_tod'

    # Initialization of averaged variation rates.
    ∂C_avg = @SVector zeros(T, 6)

    N = params.num_sampling_points_per_orbit

    # NOTE: Theoretically, we need to update the Moon and Sun positions at each sampling
    # point, but for efficiency, we assume they are constant over one orbit.
    for k in 0:(N - 1)
        # Sampling point along the mean orbit.
        M̄k = mod(M̄ + (2π * k / N), 2π)
        f̄k = mean_to_true_anomaly(ē, M̄k)

        # Position and velocity from mean orbital elements.
        rk_tod, vk_tod = _coe_to_rv(ā, ē, ī, Ω̄, ω̄, f̄k)
        rk = norm(rk_tod)

        rk_pef = D_pef_tod * rk_tod

        # Conservative perturbations in TOD.
        δak_tod =
            (D_tod_pef * _perturbational_gravity_acceleration(
                gm, jd_utc, rk_pef; P = params.gravity_P, dP = params.gravity_dP
            )) +
            _point_mass_acceleration(rk_tod, rsun_tod, _μ_SUN) +
            _point_mass_acceleration(rk_tod, rmoon_tod, _μ_MOON)

        # Perturbations acceleration in Hill frame
        Dk_hill_tod = _r_eci_to_hill(rk_tod, vk_tod)
        δak_hill    = Dk_hill_tod * δak_tod

        # Instantaneous Gauss equations for the equinoctial elements (with mean
        # parameters).
        Ak, Bk = _equinoctial_gauss_variational_matrices(
            ā, ē, ī, Ω̄, ω̄, f̄k, rk, p̄, h̄, η̄, n̄
        )

        # Accumulation
        ∂C_avg = ∂C_avg + (Ak * δak_hill + Bk)
    end

    # Average of conservative perturbations.
    ∂C_total = ∂C_avg / N

    # The J₂² rates are expressed in classical elements. Convert them to equinoctial
    # rates using the Jacobian of the transformation.
    Jec    = _classical_to_equinoctial_jacobian(ē, ī, Ω̄, ω̄)
    ∂u_J₂² = Jec * J₂²_variational_rates(ā, ē, ī, ω̄, gm)

    # Non-conservative perturbations (averaged)
    ∂u_drag, ∂u_srp = _atmospheric_drag_and_solar_radiation_pressure_variational_rates(
        jd_utc,
        ā,
        ē,
        ī,
        Ω̄,
        ω̄,
        rsun_tod,
        params
    )

    ∂C_total = ∂C_total + ∂u_drag + ∂u_srp + ∂u_J₂²

    return ∂C_total
end
