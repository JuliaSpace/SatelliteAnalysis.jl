## Description #############################################################################
#
# Compute variational rates due to non-conservative perturbations.
#
############################################################################################

"""
    _atmospheric_drag_and_solar_radiation_pressure_variational_rates(
        jd_utc::T,
        ā::T,
        ē::T,
        ī::T,
        Ω̄::T,
        ω̄::T,
        rsun_tod::SVector{3, T},
        params::NamedTuple
    ) where T<:Number -> SVector{6, T}, SVector{6, T}

Compute averaged Gauss variation rates due to atmospheric drag and solar radiation pressure
**[1]**. This routine calculates the averaged rates of change of the orbital elements under
atmospheric drag using quadrature equally spaced in true anomaly and temporal weighting
(`r²/h`).

# Arguments

- `jd_utc::T`: Julian date [UTC] at which the rates are computed.
- `ā::T`: Mean semi-major axis [m].
- `ē::T`: Mean eccentricity [-].
- `ī::T`: Mean inclination [rad].
- `Ω̄::T`: Mean right ascension of ascending node [rad].
- `ω̄::T`: Mean argument of perigee [rad].
- `rsun_tod::SVector{3, T}`: Sun position vector [m] in TOD frame.
- `params::NamedTuple`: Named tuple containing environment parameters:
    - `gm::AbstractGravityModel`: Gravity model for Earth.
    - `num_sampling_points_per_orbit::Int`: Number of quadrature points for averaging.
    - `satellite_mass::Number`: Spacecraft mass [kg].
    - `satellite_mean_area::Number`: Effective cross-sectional area [m²].
    - `j2osc_prop::OrbitPropagatorJ2Osculating`: Pre-allocated J2 osculating propagator
        used for the mean-to-osculating conversion.
    - `Ap::Function`: A function to retrieve the geomagnetic index [-]. It must have the
        signature `Ap(jd_utc::Number) -> Number`.
    - `C_d::Number`: Drag coefficient [-].
    - `C_r::Number`: Reflectivity coefficient [-].
    - `F107::Function`: A function to retrieve the solar flux index [sfu]. It must have the
        signature `F107(jd_utc::Number) -> Number`.

# Returns

- `SVector{6, T}`: Averaged Gauss rates due to atmospheric drag.
- `SVector{6, T}`: Averaged Gauss rates due to solar radiation pressure.

# References

- **[1]** Battin, R. H. (1999). An Introduction to the Mathematics and Methods of
    Astrodynamics. Revised ed. AIAA Education Series, Reston, VA.
"""
function _atmospheric_drag_and_solar_radiation_pressure_variational_rates(
    jd_utc::T,
    ā::T,
    ē::T,
    ī::T,
    Ω̄::T,
    ω̄::T,
    rsun_tod::SVector{3, T},
    params::NamedTuple
) where T<:Number
    C_d       = params.C_d
    C_r       = params.C_r
    N         = params.num_sampling_points_per_orbit
    gm        = params.gm
    mass      = params.satellite_mass
    mean_area = params.satellite_mean_area
    orbp      = params.j2osc_prop

    # Resolve the space indices once per evaluation since `jd_utc` is constant here,
    # avoiding one interpolation per sampling point.
    F107 = Float64(params.F107(jd_utc))
    Ap   = Float64(params.Ap(jd_utc))

    μ = GravityModels.gravity_constant(gm)

    # Auxiliaries with mean elements. See `_dynamics!` for the rationale of the clamps.
    ā  = max(ā, 0.9 * EARTH_EQUATORIAL_RADIUS)
    ē  = clamp(ē, 1e-6, max(1e-6, 1 - EARTH_EQUATORIAL_RADIUS / ā))
    ī  = max(ī, 1e-6)
    ē² = ē^2
    η̄² = 1 - ē²
    p̄  = ā * η̄²
    h̄  = √(μ * p̄)
    n̄  = √(μ / ā^3)
    η̄  = √η̄²

    # The rotation between TOD and PEF is constant within one evaluation. Hoisting it out
    # of the quadrature loop avoids one frame reduction per sampling point. The PEF frame
    # rotates with the Earth, so the frame angular velocity must be accounted for when
    # converting velocity vectors.
    D_tod_pef = r_ecef_to_eci(PEF(), TOD(), jd_utc)
    D_pef_tod = D_tod_pef'
    ω_pef     = @SVector T[0, 0, EARTH_ANGULAR_SPEED]

    # Initialization.
    ∂u_drag = @SVector zeros(T, 6)
    ∂u_srp  = @SVector zeros(T, 6)
    Wsum    = zero(T)

    # Quadrature in f ∈ [0, 2π) [rad].
    for k in 0:(N - 1)
        f̄k = 2π * k / N
        M̄k = mod(true_to_mean_anomaly(ē, f̄k), 2π)

        # Obtain osculating elements.
        a, e, i, Ω, ω, M = _mean_to_osculating_elements(ā, ē, ī, Ω̄, ω̄, M̄k, orbp)

        # Position and velocity in TOD.
        rk_tod, vk_tod = _coe_to_rv(a, e, i, Ω, ω, M)
        rk² = dot(rk_tod, rk_tod)
        rk  = √rk²

        # Matrix to convert TOD to Hill frame.
        D_hill_tod = _r_eci_to_hill(rk_tod, vk_tod)

        # Position and velocity in PEF to compute the atmospheric drag acceleration.
        rk_pef = D_pef_tod * rk_tod
        vk_pef = D_pef_tod * vk_tod - ω_pef × rk_pef

        # Drag acceleration in PEF.
        adrag_pef = _atmospheric_drag_acceleration(
            jd_utc,
            rk_pef,
            vk_pef,
            mean_area,
            mass,
            Ap,
            C_d,
            F107
        )

        # Drag acceleration in TOD.
        adrag_tod = D_tod_pef * adrag_pef

        # Solar radiation pressure acceleration in TOD, gated by the Earth shadow at
        # the sampling point: full acceleration under direct sunlight, half in the
        # penumbra, and none in the umbra.
        lc = lighting_condition(rk_tod, rsun_tod)
        ν  = lc == :sunlight ? T(1) : (lc == :penumbra ? T(1 // 2) : T(0))
        asrp_tod =
            ν * _solar_radiation_acceleration(rk_tod, rsun_tod, mean_area, mass, C_r)

        # Compute accelerations in Hill frame.
        adrag_hill = D_hill_tod * adrag_tod
        asrp_hill  = D_hill_tod * asrp_tod

        # Gauss equations with mean parameters. The Kepler term `B` is not added here
        # because it is already accounted for in the conservative average, avoiding
        # counting the mean motion multiple times in the mean anomaly rate.
        Ak, _ = _equinoctial_gauss_variational_matrices(
            ā, ē, ī, Ω̄, ω̄, f̄k, rk, p̄, h̄, η̄, n̄
        )

        # Temporal weighting.
        w_t = rk² / h̄

        ∂u_drag = ∂u_drag + w_t * (Ak * adrag_hill)
        ∂u_srp  = ∂u_srp  + w_t * (Ak * asrp_hill)

        Wsum += w_t
    end

    # Normalize.
    if Wsum > 0
        ∂u_drag = ∂u_drag / Wsum
        ∂u_srp  = ∂u_srp  / Wsum
    else
        ∂u_drag = @SVector zeros(T, 6)
        ∂u_srp  = @SVector zeros(T, 6)
    end

    return ∂u_drag, ∂u_srp
end

"""
    J₂²_variational_rates(
        ā::T,
        ē::T,
        ī::T,
        ω̄::T,
        gm::AbstractGravityModel
    ) where T <: Number -> SVector{6, T}

Compute the averaged variation rates of the mean classical orbital elements due to the
second-order zonal harmonic effects (J₂²) for an orbit with mean semi-major axis `ā` [m],
mean eccentricity `ē` [-], mean inclination `ī` [rad], and mean argument of perigee `ω̄`
[rad]. The J₂ coefficient is obtained from the gravity model `gm`.

The returned rates contain **only** the terms proportional to J₂². The first-order secular
rates are already captured by numerically averaging the zonal gravitational acceleration
over the mean orbit, so they must not be included here. The secular J₂² rates of `Ω`, `ω`,
and `M` follow the same analytical theory used by the J4 orbit propagator of
**SatelliteToolbox.jl** **[1]**, and the long-period rates of `e` and `i` follow **[2]**.
The long-period J₂² rates of the angular elements are neglected since they do not affect
the decay evolution.

# Returns

- `SVector{6, T}`: Averaged rates `[∂a, ∂e, ∂i, ∂Ω, ∂ω, ∂M]` due to J₂² effects, in
    [m/s; 1/s; rad/s; rad/s; rad/s; rad/s].

# References

- **[1]** Kozai, Y (1959). The Motion of a Close Earth Satellite. The Astronomical
    Journal, v. 64, no. 1274, pp. 367 -- 377.
- **[2]** Vallado, D. A (2013). *Fundamentals of Astrodynamics and Applications*. 4th ed.
    Microcosm Press, Hawthorne, CA, sec. 9.6.
"""
function J₂²_variational_rates(
    ā::T,
    ē::T,
    ī::T,
    ω̄::T,
    gm::AbstractGravityModel
) where T <: Number
    μ    = GravityModels.gravity_constant(gm)
    Re   = GravityModels.radius(gm)
    C₂_₀ = GravityModels.coefficients(gm, 2, 0) |> first
    J₂   = -C₂_₀ * √5

    # Regularization.
    ē = max(ē, 1e-6)
    ī = max(ī, 1e-6)

    # Auxiliary variables.
    ē²  = ē^2
    η̄²  = 1 - ē²
    η̄   = √η̄²
    p̄   = ā * η̄²
    n₀  = √(μ / ā^3)
    R̄e² = (Re / p̄)^2
    R̄e⁴ = R̄e²^2
    J₂² = J₂^2

    sin_ī, cos_ī = sincos(ī)
    sin_ī²  = sin_ī^2
    sin_ī⁴  = sin_ī²^2
    cos_ī⁴  = cos_ī^4
    sin_2ω̄  = sin(2ω̄)

    kn₂  = J₂  * R̄e²
    kn₂₂ = J₂² * R̄e⁴

    # Perturbed mean motion considering the J₂ and J₂² secular terms [1].
    n̄ = n₀ * (
        1 +
        (3 // 4) * kn₂ * η̄ * (2 - 3sin_ī²) +
        (3 // 128) * kn₂₂ * η̄ * (
            120 + 64η̄ - 40η̄² +
            (-240 - 192η̄ + 40η̄²) * sin_ī² +
            (105 + 144η̄ + 25η̄²) * sin_ī⁴
        )
    )

    # == Secular J₂² Rates =================================================================
    #
    # The expressions below are the full (J₂ + J₂²) secular rates [1] minus the first-order
    # J₂ secular rates evaluated with the Kepler mean motion of the mean elements, which is
    # exactly the part captured by the numerical averaging of the zonal acceleration.

    ∂M = (n̄ - n₀) - (3 // 4) * n₀ * kn₂ * η̄ * (2 - 3sin_ī²)

    ∂ω = (3 // 4) * (n̄ - n₀) * kn₂ * (4 - 5sin_ī²) +
        (3 // 128) * n̄ * kn₂₂ * (
            384 + 96ē² - 384η̄ +
            (-824 - 116ē² + 1056η̄) * sin_ī² +
            (430 - 5ē² - 720η̄) * sin_ī⁴
        ) -
        (15 // 16) * n₀ * kn₂₂ * ē² * cos_ī⁴

    ∂Ω = -(3 // 2) * (n̄ - n₀) * kn₂ * cos_ī +
        (3 // 32) * n̄ * kn₂₂ * cos_ī * (
            -36 - 4ē² + 48η̄ + (40 - 5ē² - 72η̄) * sin_ī²
        )

    # == Long-Period J₂² Rates =============================================================
    #
    # The pair (∂e, ∂i) satisfies the exact invariant of zonal fields H = √(μ p̄) cos ī,
    # i.e., ∂i = -(ē ∂e cos ī) / (η̄² sin ī).

    ∂e = -(3 // 32) * n₀ * kn₂₂ * sin_ī² * (14 - 15sin_ī²) * ē * η̄² * sin_2ω̄

    ∂i = +(3 // 64) * n₀ * kn₂₂ * sin(2ī) * (14 - 15sin_ī²) * ē² * sin_2ω̄

    ∂a = zero(T)

    return @SVector T[∂a, ∂e, ∂i, ∂Ω, ∂ω, ∂M]
end
