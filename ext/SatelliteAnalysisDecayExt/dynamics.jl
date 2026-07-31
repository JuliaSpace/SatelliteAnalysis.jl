## Description #############################################################################
#
# Dynamics function to compute the satellite decay.
#
############################################################################################

"""
    _dynamics!(
        du::AbstractVector{T},
        state_mean::AbstractVector{T},
        params::Dict,
        t::Real
    ) -> Nothing

Compute Gauss mean-element derivatives with averaged perturbations.

This routine computes the time derivatives of the mean orbital elements using Gauss'
variational equations. It evaluates instantaneous forces on the mean-reference orbit for
conservative perturbations (Earth gravity, third bodies) and adds temporally averaged rates
for non-conservative perturbations (Atmospheric Drag and Solar Radiation Pressure).

!!! note

    - Conservative perturbations are evaluated on the mean orbit at evenly spaced mean anomalies.
    - Non-conservative perturbations are added via averaged routines (drag and SRP).
    - Units must be consistent: meters, seconds, kilograms.

# Keywords

- `du::AbstractVector{T}`: Output derivative vector `[da, de, di, dΩ, dω, dM]`.
- `state_mean::AbstractVector{T}`: Mean orbital elements `[a, e, i, Ω, ω, M]`.
- `params::Dict`: Dictionary with environment/model parameters (e.g., `jd0_utc`, `F107`, `Ap`).
- `t::Real`: Time since epoch [s].

# References

- **[1]** Battin, R. H. (1999). An Introduction to the Mathematics and Methods of
    Astrodynamics. Revised ed. AIAA Education Series, Reston, VA.
"""
function _dynamics!(
    du::AbstractVector{T},
    u::AbstractVector{T},
    params,
    t::Real
) where T <: Number
    # Gravity model.
    gm = params.gm
    μ  = GravityModels.gravity_constant(gm)

    # Time and ephemerides.
    jd₀_utc = params.jd₀_utc
    jd_utc  = jd₀_utc + t / 86400.0

    # Unpack mean orbital elements.
    ā, ē, ī, Ω̄, ω̄, M̄ = _equinoctial_to_classical(u)
    ē = max(ē, 1e-6)
    ī = max(ī, 1e-6)

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

    # Convert the vectors to PEF.
    D_pef_tod = r_eci_to_ecef(TOD(), PEF(), jd_utc)
    D_tod_pef = r_ecef_to_eci(PEF(), TOD(), jd_utc)

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
        rk_tod, vk_tod = _coe_to_rv(ā, ē, ī, Ω̄, ω̄, M̄k)
        rk = norm(rk_tod)

        rk_pef = D_pef_tod * rk_tod

        # Conservative perturbations in TOD.
        δak_tod =
            (D_tod_pef * _perturbational_gravity_acceleration(gm, jd_utc, rk_pef)) +
            _point_mass_acceleration(rk_tod, rsun_tod, _μ_SUN) +
            _point_mass_acceleration(rk_tod, rmoon_tod, _μ_MOON)

        # Perturbations acceleration in Hill frame
        Dk_hill_tod = _r_eci_to_hill(rk_tod, vk_tod)
        δak_hill    = Dk_hill_tod * δak_tod

        # Instantaneous Gauss equations (with mean parameters)
        Ak, Bk = _gauss_variational_matrices(ā, ē, ī, ω̄, f̄k, rk, p̄, h̄, η̄, n̄)

        # Accumulation
        ∂C_avg = ∂C_avg + (Ak * δak_hill + Bk)
    end

    # Average of conservative perturbations.
    ∂C_total = ∂C_avg / N

    ∂u_J₂² = J₂²_variational_rates(ā, ē, ī, ω̄, gm)

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

    Jec = _coe_to_rv_jacobian(ē, ī, Ω̄, ω̄)

    du[1:6] .= Jec * ∂C_total

    return nothing
end
