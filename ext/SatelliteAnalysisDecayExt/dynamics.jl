## Description #############################################################################
#
# Dynamics function to compute the satellite decay.
#
############################################################################################

"""
    _dynamics(u::AbstractVector{T}, params::NamedTuple, t::Real) where T <: Number -> SVector{6, T}

Compute the time derivatives of the mean equinoctial orbital elements `u` at the time `t`
[s] after the epoch using Gauss variational equations with averaged perturbations.

This routine averages the rates caused by the conservative perturbations (Earth gravity and
third bodies) and by the non-conservative perturbations (atmospheric drag and solar
radiation pressure) over one orbit, and adds the closed-form J₂² rates.

!!! note

    - All the perturbations are averaged using a single quadrature with points equally
      spaced in the true anomaly and the temporal weighting `r² / h`. Hence, the sampling
      is denser near the perigee, where the perturbations are stronger.
    - Conservative perturbations are evaluated on the mean orbit.
    - The accelerations of the non-conservative perturbations are evaluated on the
      osculating orbit obtained from the mean elements, but they are mapped to the rates
      using the Gauss matrix, the Hill frame, and the weighting of the mean orbit (see
      `_atmospheric_drag_and_solar_radiation_pressure_rates`).
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
function _dynamics(u::AbstractVector{T}, params, t::Real) where {T <: Number}
    # Gravity model and its constants, hoisted into `params` by `_decay_analysis`.
    gm = params.gm
    μ  = params.μ

    # Time and ephemerides.
    jd₀_utc = params.jd₀_utc
    jd_utc  = jd₀_utc + t / 86400.0

    # Unpack mean orbital elements.
    ā, ē, ī, Ω̄, ω̄, _ = _equinoctial_to_classical(u)

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

    # Resolve the space indices once per evaluation since `jd_utc` is constant here,
    # avoiding one interpolation per sampling point.
    space_indices = params.space_indices(jd_utc)

    # Initialization of the weighted sums of the variation rates.
    ∂C_sum  = @SVector zeros(T, 6)
    ∂u_drag = @SVector zeros(T, 6)
    ∂u_srp  = @SVector zeros(T, 6)
    Wsum    = zero(T)

    N = params.num_sampling_points_per_orbit

    # Quadrature in f ∈ [0, 2π) [rad]. The time average is obtained using the temporal
    # weighting `dt = (r² / h) df`.
    #
    # NOTE: Theoretically, we need to update the Moon and Sun positions at each sampling
    # point, but for efficiency, we assume they are constant over one orbit.
    for k in 0:(N - 1)
        # Sampling point along the mean orbit.
        f̄k = T(2π) * k / N

        # Position and velocity from mean orbital elements.
        rk_tod, vk_tod = _coe_to_rv(ā, ē, ī, Ω̄, ω̄, f̄k)
        rk² = dot(rk_tod, rk_tod)
        rk = √rk²

        rk_pef = D_pef_tod * rk_tod

        # Conservative perturbations in TOD.
        δak_tod =
            (
                D_tod_pef * _perturbational_gravity_acceleration(
                    gm, jd_utc, rk_pef; μ = μ, workspace = params.gravity_workspace
                )
            ) +
            _point_mass_acceleration(rk_tod, rsun_tod, _μ_SUN) +
            _point_mass_acceleration(rk_tod, rmoon_tod, _μ_MOON)

        # Perturbations acceleration in Hill frame
        Dk_hill_tod = _r_eci_to_hill(rk_tod, vk_tod)
        δak_hill    = Dk_hill_tod * δak_tod

        # Instantaneous Gauss equations for the equinoctial elements (with mean
        # parameters).
        Ak = _equinoctial_gauss_variational_matrices(ā, ē, ī, Ω̄, ω̄, f̄k, rk, p̄, h̄, η̄)

        # Temporal weighting.
        w_t = rk² / h̄

        # Accumulation.
        ∂C_sum += w_t * (Ak * δak_hill)
        Wsum   += w_t

        # Non-conservative perturbations at the same sampling point, using the same Gauss
        # matrix, Hill frame, and temporal weighting.
        ∂u_drag_k, ∂u_srp_k = _atmospheric_drag_and_solar_radiation_pressure_rates(
            jd_utc,
            ā,
            ē,
            ī,
            Ω̄,
            ω̄,
            f̄k,
            Ak,
            Dk_hill_tod,
            rsun_tod,
            D_pef_tod,
            D_tod_pef,
            space_indices,
            params,
        )

        ∂u_drag += w_t * ∂u_drag_k
        ∂u_srp  += w_t * ∂u_srp_k
    end

    # Average of the conservative perturbations plus the constant Kepler term of the
    # mean longitude rate. The sum of weights is strictly positive since
    # `w_t = rk² / h̄ > 0`.
    ∂C_total = ∂C_sum / Wsum + @SVector T[0, n̄, 0, 0, 0, 0]

    # Average of the non-conservative perturbations.
    ∂u_drag = ∂u_drag / Wsum
    ∂u_srp  = ∂u_srp / Wsum

    # The J₂² rates are expressed in classical elements. Convert them to equinoctial
    # rates using the Jacobian of the transformation.
    Jec    = _classical_to_equinoctial_jacobian(ē, ī, Ω̄, ω̄)
    ∂u_J₂² = Jec * J₂²_variational_rates(ā, ē, ī, ω̄, μ, params.Re, params.J₂)

    ∂C_total = ∂C_total + ∂u_drag + ∂u_srp + ∂u_J₂²

    return ∂C_total
end
