## Description #############################################################################
#
# Compute variational rates due to non-conservative perturbations.
#
############################################################################################

"""
    _atmospheric_drag_and_solar_radiation_pressure_rates(
        jd_utc::T,
        ā::T,
        ē::T,
        ī::T,
        Ω̄::T,
        ω̄::T,
        f̄::T,
        A::StaticMatrix{6, 3, T},
        D_hill_tod::StaticMatrix{3, 3, T},
        rsun_tod::SVector{3, T},
        D_pef_tod::StaticMatrix{3, 3, T},
        D_tod_pef::StaticMatrix{3, 3, T},
        space_indices::NamedTuple,
        params::NamedTuple
    ) where T <: Number -> SVector{6, T}, SVector{6, T}

Compute the Gauss variation rates due to the atmospheric drag and the solar radiation
pressure **[1]** at the sampling point of the orbit-averaging quadrature with mean true
anomaly `f̄` [rad].

The accelerations are computed using the osculating position and velocity obtained from the
mean elements because the atmospheric density, the velocity relative to the atmosphere, and
the Earth shadow depend on the actual satellite position. However, the variational
equations are written for the mean elements. Hence, the accelerations are projected onto
the Hill frame of the mean orbit and mapped to the rates using the Gauss matrix of the mean
orbit, which are computed by the caller, consistently with the conservative perturbations.

The mean elements and the auxiliaries passed to this function must already be clamped to
the physically meaningful region and derived consistently, as done by `_dynamics`.

# Arguments

- `jd_utc::T`: Julian date [UTC] at which the rates are computed.
- `ā::T`: Mean semi-major axis [m].
- `ē::T`: Mean eccentricity [-].
- `ī::T`: Mean inclination [rad].
- `Ω̄::T`: Mean right ascension of ascending node [rad].
- `ω̄::T`: Mean argument of perigee [rad].
- `f̄::T`: Mean true anomaly of the sampling point [rad], which must be in `[0, 2π)`.
- `A::StaticMatrix{6, 3, T}`: Gauss matrix of the mean orbit at the sampling point,
    computed by `_equinoctial_gauss_variational_matrices`.
- `D_hill_tod::StaticMatrix{3, 3, T}`: DCM that rotates vectors from the TOD frame to the
    Hill frame of the mean orbit at the sampling point.
- `rsun_tod::SVector{3, T}`: Sun position vector [m] in TOD frame.
- `D_pef_tod::StaticMatrix{3, 3, T}`: DCM that rotates vectors from the TOD frame to the
    PEF frame at `jd_utc`.
- `D_tod_pef::StaticMatrix{3, 3, T}`: DCM that rotates vectors from the PEF frame to the
    TOD frame at `jd_utc`.
- `space_indices::NamedTuple`: Space indices at `jd_utc` required by the atmospheric model.
- `params::NamedTuple`: Named tuple containing environment parameters:
    - `atmospheric_model::Any`: Callable object that computes the atmospheric density
        [kg/m³] at a given location and time considering a set of space indices. It must
        have the signature
        `(jd_utc::Number, lat::Number, lon::Number, alt::Number, space_indices::NamedTuple) -> Number`.
    - `satellite_mass::Number`: Spacecraft mass [kg].
    - `satellite_mean_area::Number`: Effective cross-sectional area [m²].
    - `j2osc_prop::OrbitPropagatorJ2Osculating`: Pre-allocated J2 osculating propagator
        used for the mean-to-osculating conversion.
    - `C_d::Number`: Drag coefficient [-].
    - `C_r::Number`: Solar radiation pressure coefficient [-].

# Returns

- `SVector{6, T}`: Gauss rates due to atmospheric drag.
- `SVector{6, T}`: Gauss rates due to solar radiation pressure.

# References

- **[1]** Battin, R. H. (1999). An Introduction to the Mathematics and Methods of
    Astrodynamics. Revised ed. AIAA Education Series, Reston, VA.
"""
function _atmospheric_drag_and_solar_radiation_pressure_rates(
    jd_utc::T,
    ā::T,
    ē::T,
    ī::T,
    Ω̄::T,
    ω̄::T,
    f̄::T,
    A::StaticMatrix{6, 3, T},
    D_hill_tod::StaticMatrix{3, 3, T},
    rsun_tod::SVector{3, T},
    D_pef_tod::StaticMatrix{3, 3, T},
    D_tod_pef::StaticMatrix{3, 3, T},
    space_indices::NamedTuple,
    params::NamedTuple,
) where {T <: Number}
    atmospheric_model = params.atmospheric_model
    C_d               = params.C_d
    C_r               = params.C_r
    mass              = params.satellite_mass
    mean_area         = params.satellite_mean_area
    orbp              = params.j2osc_prop

    # The PEF frame rotates with the Earth, so the frame angular velocity must be accounted
    # for when converting velocity vectors.
    ω_pef = @SVector T[0, 0, EARTH_ANGULAR_SPEED]

    # Position and velocity of the osculating orbit in TOD.
    r_tod, v_tod = _mean_to_osculating_rv(ā, ē, ī, Ω̄, ω̄, f̄, orbp)

    # Position and velocity in PEF to compute the atmospheric drag acceleration.
    r_pef = D_pef_tod * r_tod
    v_pef = D_pef_tod * v_tod - ω_pef × r_pef

    # Drag acceleration in PEF.
    adrag_pef = _atmospheric_drag_acceleration(
        atmospheric_model, jd_utc, r_pef, v_pef, mean_area, mass, C_d, space_indices
    )

    # Drag acceleration in TOD.
    adrag_tod = D_tod_pef * adrag_pef

    # Solar radiation pressure acceleration in TOD, gated by the Earth shadow at the
    # sampling point: full acceleration under direct sunlight, half in the penumbra, and
    # none in the umbra.
    lc = lighting_condition(r_tod, rsun_tod)
    ν  = lc == :sunlight ? T(1) : (lc == :penumbra ? T(1 // 2) : T(0))

    asrp_tod = ν * _solar_radiation_acceleration(r_tod, rsun_tod, mean_area, mass, C_r)

    # Compute accelerations in the Hill frame of the mean orbit.
    adrag_hill = D_hill_tod * adrag_tod
    asrp_hill  = D_hill_tod * asrp_tod

    # Gauss equations with mean parameters. The Kepler term is not added here because it is
    # already accounted for in `_dynamics`, avoiding counting the mean motion multiple times
    # in the mean anomaly rate.
    return A * adrag_hill, A * asrp_hill
end

"""
    J₂²_variational_rates(
        ā::T,
        ē::T,
        ī::T,
        ω̄::T,
        μ::Number,
        Re::Number,
        J₂::Number
    ) where T <: Number -> SVector{6, T}

Compute the averaged variation rates of the mean classical orbital elements due to the
second-order zonal harmonic effects (J₂²) for an orbit with mean semi-major axis `ā` [m],
mean eccentricity `ē` [-], mean inclination `ī` [rad], and mean argument of perigee `ω̄`
[rad]. The gravitational constant `μ` [m³/s²], the reference radius `Re` [m], and the J₂
coefficient `J₂` [-] must be consistent with the gravity model used by the caller.

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
    ā::T, ē::T, ī::T, ω̄::T, μ::Number, Re::Number, J₂::Number
) where {T <: Number}
    # Regularization.
    ē = max(ē, T(1e-6))
    ī = max(ī, T(1e-6))

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
    sin_ī² = sin_ī^2
    sin_ī⁴ = sin_ī²^2
    cos_ī⁴ = cos_ī^4
    sin_2ω̄ = sin(2ω̄)

    kn₂  = J₂ * R̄e²
    kn₂₂ = J₂² * R̄e⁴

    # Perturbed mean motion considering the J₂ and J₂² secular terms [1].
    n̄ =
        n₀ * (
            1 +
            (3 // 4) * kn₂ * η̄ * (2 - 3sin_ī²) +
            (3 // 128) *
            kn₂₂ *
            η̄ *
            (
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

    ∂ω =
        (3 // 4) * (n̄ - n₀) * kn₂ * (4 - 5sin_ī²) +
        (3 // 128) *
        n̄ *
        kn₂₂ *
        (384 + 96ē² - 384η̄ + (-824 - 116ē² + 1056η̄) * sin_ī² + (430 - 5ē² - 720η̄) * sin_ī⁴) -
        (15 // 16) * n₀ * kn₂₂ * ē² * cos_ī⁴

    ∂Ω =
        -(3 // 2) * (n̄ - n₀) * kn₂ * cos_ī +
        (3 // 32) * n̄ * kn₂₂ * cos_ī * (-36 - 4ē² + 48η̄ + (40 - 5ē² - 72η̄) * sin_ī²)

    # == Long-Period J₂² Rates =============================================================
    #
    # The pair (∂e, ∂i) satisfies the exact invariant of zonal fields H = √(μ p̄) cos ī,
    # i.e., ∂i = -(ē ∂e cos ī) / (η̄² sin ī).

    ∂e = -(3 // 32) * n₀ * kn₂₂ * sin_ī² * (14 - 15sin_ī²) * ē * η̄² * sin_2ω̄

    ∂i = +(3 // 64) * n₀ * kn₂₂ * sin(2ī) * (14 - 15sin_ī²) * ē² * sin_2ω̄

    ∂a = zero(T)

    return @SVector T[∂a, ∂e, ∂i, ∂Ω, ∂ω, ∂M]
end
