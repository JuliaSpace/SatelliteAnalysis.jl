## Description #############################################################################
#
# Compute accelerations due to various perturbative forces acting on a satellite in orbit.
#
############################################################################################

"""
    _atmospheric_drag_acceleration(
        jd_utc::Number,
        r_ecef::AbstractVector{T},
        v_ecef::AbstractVector{T},
        area::Number,
        mass::Number,
        Ap::Number,
        Cd::Number,
        F107::Number
    ) where T <: Number -> SVector{3, T}

Compute the acceleration [m/s²] due to atmospheric drag in the ECEF frame using the
NRLMSISE-00 model.

# Arguments

- `jd_utc::Number`: Julian date [UTC] in which the atmospheric drag will be computed.
- `r_ecef::AbstractVector{T}`: Satellite position vector [m] in ECEF frame.
- `v_ecef::AbstractVector{T}`: Satellite velocity vector [m/s] in ECEF frame.
- `area::Number`: Effective cross-sectional area [m²] exposed to atmosphere.
- `mass::Number`: Spacecraft mass [kg].
- `Ap::Number`: Geomagnetic index [-].
- `Cd::Number`: Drag coefficient [-].
- `F107::Number`: Solar flux index [sfu].

# Returns

- `SVector{3, T}`: Drag acceleration vector [m/s²] represented in the ECEF frame.
"""
function _atmospheric_drag_acceleration(
    jd_utc::Number,
    r_ecef::AbstractVector{T},
    v_ecef::AbstractVector{T},
    area::Number,
    mass::Number,
    Ap::Number,
    Cd::Number,
    F107::Number
) where T <: Number

    lat, lon, h = ecef_to_geodetic(r_ecef)

    # The altitude can be negative when the integrator evaluates trial stages with
    # unphysical states near the decay end. Clamp it to keep the atmospheric model valid so
    # that the error control can reject the step.
    h = max(h, 0.0)

    atmos = AtmosphericModels.nrlmsise00(jd_utc, h, lat, lon, F107, F107, Ap)
    ρ     = atmos.total_density

    a_drag_ecef = -(1 // 2) * T(Cd) * (T(area) / T(mass)) * T(ρ) * norm(v_ecef) .* v_ecef

    return a_drag_ecef
end

"""
    _perturbational_gravity_acceleration(
        gm::AbstractGravityModel,
        jd_utc::Number,
        r_ecef::AbstractVector;
        kwargs...
    ) -> SVector{3, Number}

Compute the perturbational acceleration [m/s²] due to Earth's gravity field in the ECI
frame. This algorithm obtains the acceleration from the gravity model `gm` at the ECEF
position `r_ecef` [m] and the Julian date `jd_utc` [UTC]. This algorithm obtains the
perturbed acceleration by computing the gravitational acceleration using `max_degree` and
`max_order` and subtracting the central term (degree 0, order 0) to isolate the
perturbational part.

# Keywords

- `max_degree::Int`: Maximum degree of the gravity model to consider.
    (**Default**: 7).
- `max_order::Int`: Maximum order of the gravity model to consider.
    (**Default**: 0).

# Returns

- `SVector{3, Number}`: Perturbed acceleration vector [m/s²] represented in the ECEF frame.
    The vector type is the same as the input gravity model.
"""
function _perturbational_gravity_acceleration(
    gm::AbstractGravityModel,
    jd_utc::Number,
    r_ecef::AbstractVector;
    max_degree::Int = 7,
    max_order::Int = 0,
)
    # Total acceleration.
    a_ecef_total = GravityModels.gravitational_acceleration(
        gm,
        r_ecef,
        jd_utc;
        max_degree = max_degree,
        max_order = max_order
    )

    # Central term (degree 0, order 0), computed analytically as -μ r / r³ to avoid a
    # second gravity model evaluation.
    μ  = GravityModels.gravity_constant(gm)
    r² = dot(r_ecef, r_ecef)
    a_ecef_central = -μ / (r² * √r²) * r_ecef

    # Perturbational acceleration.
    a_ecef_pert = a_ecef_total - a_ecef_central

    return a_ecef_pert
end

"""
    _point_mass_acceleration(
        rsat_eci::AbstractVector{T},
        rbody_eci::AbstractVector{T},
        body_μ::Number
    ) where T <: Number -> SVector{3, T}

Compute the acceleration [m/s²] due to a third-body point mass (e.g. Sun, Moon, or another
planet) on a satellite in the Earth-centered inertial (ECI) frame. `rsat_eci` is the
satellite position vector [m] in ECI, `rbody_eci` is the point mass position vector [m] in
ECI, and `body_μ` is the gravitational parameter of the point mass [m³/s²].

# Returns

- `SVector{3, T}`: Perturbational acceleration vector [m/s²].
"""
function _point_mass_acceleration(
    rsat_eci::AbstractVector{T},
    rbody_eci::AbstractVector{T},
    body_μ::Number
) where T <: Number

    # Relative position vector of satellite w.r.t. point mass.
    Δr_eci = rsat_eci .- rbody_eci

    a_eci = -T(body_μ) * (Δr_eci / norm(Δr_eci)^3 + rbody_eci / norm(rbody_eci)^3)

    return a_eci
end

"""
    _solar_radiation_acceleration(
        rsat_eci::AbstractVector{T},
        rsun_eci::AbstractVector{T},
        area::Number,
        mass::Number,
        C_r::Number
    ) where T <: Number -> SVector{3, T}

Compute the acceleration due to Solar Radiation Pressure **[1]**.

The acceleration is computed as the effect of solar radiation pressure on the spacecraft,
scaled by the inverse-square law with respect to the Sun–spacecraft distance. The worst-case
assumption is used, where the effective area is always normal to the Sun direction.

# Arguments

- `rsat_eci::AbstractVector{T}`: Spacecraft position vector [m] in an ECI frame.
- `rsun_eci::AbstractVector{T}`: Sun position vector [m] in an ECI frame.
- `area::Number`: Effective cross-sectional area [m²] exposed to Sun.
- `mass::Number`: Spacecraft mass [kg].
- `C_r::Number`: Solar radiation pressure coefficient [-].

# Returns

- `SVector{3, T}`: Perturbational acceleration vector [m/s²].

# Extended help

## References

- **[1]** Vallado, D. A. (2013). *Fundamentals of Astrodynamics and Applications*. 4th ed.
    Microcosm Press, Hawthorne, CA.
"""
function _solar_radiation_acceleration(
    rsat_eci::AbstractVector{T},
    rsun_eci::AbstractVector{T},
    area::Number,
    mass::Number,
    C_r::Number
) where T <: Number

    # Relative position vector of spacecraft w.r.t. Sun
    Δr_eci = rsat_eci .- rsun_eci

    # Acceleration due to solar radiation pressure
    a_eci = C_r * T(area) / T(mass) *
        _SOLAR_PRESSURE_1AU *
        (ASTRONOMICAL_UNIT^2) *
        (Δr_eci / norm(Δr_eci)^3)

    return a_eci
end
