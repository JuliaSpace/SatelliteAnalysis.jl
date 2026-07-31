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
    - `Ap::Number`: Geomagnetic index.
    - `C_d::Number`: Drag coefficient [-].
    - `C_r::Number`: Reflectivity coefficient [-].
    - `F107::Number`: Solar flux index.

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
    Ap        = params.Ap
    C_d       = params.C_d
    C_r       = params.C_r
    F107      = params.F107
    N         = params.num_sampling_points_per_orbit
    gm        = params.gm
    mass      = params.satellite_mass
    mean_area = params.satellite_mean_area

    μ = GravityModels.gravity_constant(gm)

    # Auxiliaries with mean elements.
    ē  = max(ē, 1e-6)
    ī  = max(ī, 1e-6)
    ē² = ē^2
    η̄² = 1 - ē²
    p̄  = ā * η̄²
    h̄  = √(μ * p̄)
    n̄  = √(μ / ā^3)
    η̄  = √η̄²

    D_tod_pef = r_ecef_to_eci(PEF(), TOD(), jd_utc)

    # Initialization.
    ∂u_drag = @SVector zeros(T, 6)
    ∂u_srp  = @SVector zeros(T, 6)
    Wsum    = zero(T)

    # Quadrature in f ∈ [0, 2π) [rad].
    for k in 0:(N - 1)
        f̄k = 2π * k / N
        M̄k = mod(true_to_mean_anomaly(ē, f̄k), 2π)

        # Obtain osculating elements.
        a, e, i, Ω, ω, M = _mean_to_osculating_elements(ā, ē, ī, Ω̄, ω̄, M̄k)

        # Position and velocity in TOD.
        rk_tod, vk_tod = _coe_to_rv(a, e, i, Ω, ω, M)
        rk² = dot(rk_tod, rk_tod)
        rk  = √rk²

        # Matrix to convert TOD to Hill frame.
        D_hill_tod = _r_eci_to_hill(rk_tod, vk_tod)

        # Position and velocity in PEF to compute the atmospheric drag acceleration.
        svk_tod = OrbitStateVector(jd_utc, rk_tod, vk_tod)
        svk_pef = sv_eci_to_ecef(svk_tod, TOD(), PEF(), jd_utc)
        rk_pef  = svk_pef.r
        vk_pef  = svk_pef.v

        # Drag acceleration in PEF.
        adrag_pef = _atmospheric_drag_acceleration(
            jd_utc,
            rk_pef,
            vk_pef,
            mean_area,
            mass,
            C_d;
            F107 = F107,
            Ap = Ap
        )

        # Drag acceleration in TOD.
        adrag_tod = D_tod_pef * adrag_pef

        # Solar radiation pressure acceleration in TOD.
        asrp_tod = _solar_radiation_acceleration(rk_tod, rsun_tod, mean_area, mass, C_r)

        # Compute accelerations in Hill frame.
        adrag_hill = D_hill_tod * adrag_tod
        asrp_hill  = D_hill_tod * asrp_tod

        # Gauss equations with mean parameters
        Ak, Bk = _gauss_variational_matrices(ā, ē, ī, ω̄, f̄k, rk, p̄, h̄, η̄, n̄)

        # Temporal weighting.
        w_t = rk² / h̄

        ∂u_drag = ∂u_drag + w_t * (Ak * adrag_hill + Bk)
        ∂u_srp  = ∂u_srp  + w_t * (Ak * asrp_hill  + Bk)

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
        ā::T,
        ē::T,
        ī::T,
        ω̄::T,
        gm::AbstractGravityModel
    ) where T<:Number -> SVector{6, T}

Compute averaged Gauss variation rates due to Earth's J₂² mean effects.

This routine calculates the mean (averaged) rates of change of the orbital elements under
second-order zonal harmonic effects (J2 squared), using the closed-form expressions for the
long-term (secular and long-period) variations.

# Arguments

- `ā::T`: Mean semi-major axis [m].
- `ē::T`: Mean eccentricity [-].
- `ī::T`: Mean inclination [rad].
- `ω̄::T`: Mean argument of perigee [rad].
- `gm::AbstractGravityModel`: Gravity model for Earth.

# Returns

- `SVector{6, T}`: Averaged Gauss rates due to J₂² effects.
"""
function J₂²_variational_rates(
    ā::T,
    ē::T,
    ī::T,
    ω̄::T,
    gm::AbstractGravityModel
) where T<:Number
    μ    = GravityModels.gravity_constant(gm)
    Re   = GravityModels.radius(gm)
    C₂_₀ = GravityModels.coefficients(gm, 2, 0) |> first
    J₂   = -C₂_₀ * √5

    # Regularization.
    ē = max(ē, 1e-6)
    ī = max(ī, 1e-6)

    # Auxiliary variables.
    ē²    = ē^2
    ē⁴    = ē²^2
    η̄²    = 1 - ē²
    η̄     = √η̄²
    p̄     = ā * η̄²
    n̄     = √(μ / ā^3)
    R̄e    = Re / p̄
    R̄e²   = R̄e^2
    R̄e⁴   = R̄e²^2
    J₂²   = J₂^2
    k₁     = n̄ * J₂² * R̄e⁴

    sin_ī, cos_ī   = sincos(ī)
    sin_ī²         = sin_ī^2
    sin_ī⁴         = sin_ī²^2
    sin_2ω̄, cos_2ω̄ = sincos(2ω̄)

    # Mean J₂² rates.
    ∂a = T(0)

    ∂e = (-3 // 2) * k₁ * sin_ī² * (14 - 15 * sin_ī^2) * ē * η̄² * sin_2ω̄

    ∂i = (+3 // 64) * k₁ * sin(2ī) * (14 - 15 * sin_ī²) * ē² * sin_2ω̄

    ∂Ω = (-3 //  2) * k₁ * cos_ī * (
        ((9 // 4) + (3 // 2) * η̄) -
        ((5 // 2) + (9 // 4) * η̄) * sin_ī² +
        (1 // 4) * (1 + (5 // 4) * sin_ī²) * ē²
    ) + (-3 // 16) * n̄ * J₂² * R̄e⁴ * cos_ī * (7 - 15sin_ī²) * ē² * cos_2ω̄

    ∂ω = (+3 // 4) * k₁ * (
        12 - (103 // 4) * sin_ī² + (215 // 16) * sin_ī⁴ + (
            (7 // 4) - (9 // 8) * sin_ī² - (45 // 32) * sin_ī⁴
        ) * ē^2 + (3 // 2) * (1 - (3 // 2) * sin_ī²) * (4 - 5 * sin_ī²) * η̄
    ) + (3 // 64) * k₁ * (
            -2 * (14 - 15 * sin_ī²) * sin_ī² + (28 - 158 * sin_ī² + 135 * sin_ī⁴) * ē^2
        ) * cos_2ω̄

    ∂M = (+3 // 2) * n̄ * R̄e² * (1 - (3 // 2) * sin_ī²) * η̄² * (-15 // 8) * k₁ * (
        (-1 + (5 // 2) * sin_ī² - (13 // 8) * sin_ī⁴) +
        (1 // 2) * ē² * (-1 + sin_ī² + (5 // 8) * sin_ī⁴)
    ) + (+3 //  2) * k₁ * (1 - (3 // 2) * sin_ī²)^2 * η̄² +
        (-9 // 64) * k₁ * sin_ī² * (14 - 15 * sin_ī²) * ē² * η̄ * cos_2ω̄ +
        (+3 // 32) * k₁ * sin_ī² * (14 - 15 * sin_ī²) * η̄^3 * cos_2ω̄ +
        (+9 //  8) * k₁ / η̄ * (
            (3 - (15 // 2) * sin_ī² + (47 // 8) * sin_ī⁴) +
            ((3 // 2) - 5 * sin_ī² + (117 // 16) * sin_ī⁴) * ē² +
            (1 // 8) * (1 + 5 * sin_ī² - (101 // 8) * sin_ī⁴) * ē² +
            (1 // 24) * sin_ī² * (
                (70 - 123 * sin_ī²) * ē² +
                2 * (28 - 33 * sin_ī²) * ē⁴
            ) * cos_2ω̄ +
            (9 // 128) * ē⁴ * sin_ī⁴ * cos(4ω̄)
        )

    return @SVector T[∂a, ∂e, ∂i, ∂Ω, ∂ω, ∂M]
end
