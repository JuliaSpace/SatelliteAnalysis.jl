## Description #############################################################################
#
# Functions related to ground repeating orbits.
#
## Remarks #################################################################################
#
# A ground repeating orbit is any orbit that the number of revolutions per day is a rational
# number. Hence, this type of orbit repeats its ground trace after a finite number of days.
#
############################################################################################

export ground_repeating_orbit_adjacent_track_angle
export ground_repeating_orbit_adjacent_track_distance

"""
    ground_repeating_orbit_adjacent_track_angle(a::T1, e::T2, i::T3, orbit_cycle::Integer; kwargs...) where {T1 <: Number, T2 <: Number, T3 <: Number}

Compute the adjacent track angle [rad] at Equator in a ground repeating orbit measured from
the satellite position. The orbit is described by its semi-major axis `a` [m], eccentricity
[ ], inclination `i` [rad], and orbit cycle `orbit_cyle` [day].

!!! warning

    The code does not check if the orbit is ground-repeating with `orbit_cycle` [day].

!!! note

    Internally, this function uses the precision obtained by promoting `T1`, `T2`, and `T3`
    to a float-pointing number `T`.

# Keywords

- `perturbation::Symbol`: Symbol to select the perturbation terms that will be used. It can
    be `:J0`, `:J2`, or `:J4`.
    (**Default**: `:J2`)
- `m0::Number`: Standard gravitational parameter for Earth [m³ / s²].
    (**Default**: `GM_EARTH`)
- `J2::Number`: J₂ perturbation term.
    (**Default**: `EGM_2008_J2`)
- `J4::Number`: J₄ perturbation term.
    (**Default**: `EGM_2008_J4`)
- `R0::Number`: Earth's equatorial radius [m].
    (**Default**: `EARTH_EQUATORIAL_RADIUS`)
- `we::Number`: Earth's angular speed [rad / s].
    (**Default**: `EARTH_ANGULAR_SPEED`)

# Extended help

A ground repeating orbit is any orbit that the number of revolutions per day is a rational
number. Hence, this type of orbit repeats its ground trace after a finite number of days.

The information `orbit_cyle` is redundant given that we have `a`, `e`, and `i`. However,
it is necessary to improve the algorithm precision. Otherwise, the `orbit_cycle` must be
obtained by computing the orbit period using `a`, `e`, and `i` and then converting it to a
rational number, leading to numerical problems.
"""
function ground_repeating_orbit_adjacent_track_angle(
    a::T1,
    e::T2,
    i::T3,
    orbit_cycle::Integer;
    perturbation::Symbol = :J2,
    J2::Number = EGM_2008_J2,
    J4::Number = EGM_2008_J4,
    R0::Number = EARTH_EQUATORIAL_RADIUS,
    m0::Number = GM_EARTH,
    we::Number = EARTH_ANGULAR_SPEED
) where {T1 <: Number, T2 <: Number, T3 <: Number}
    T   = float(promote_type(T1, T2, T3))
    R₀  = T(R0)
    ω_e = T(EARTH_ANGULAR_SPEED)

    # Compute the orbital period [s].
    ΔT = orbital_period(
        a,
        e,
        i;
        perturbation = perturbation,
        J2 = J2,
        J4 = J4,
        R0 = R0,
        m0 = m0
    )

    # Compute the RAAN time derivative [rad / s].
    ∂Ω_∂t = raan_time_derivative(
        a,
        e,
        i;
        perturbation = perturbation,
        J2 = J2,
        J4 = J4,
        m0 = m0,
        R0 = R0
    )

    # Angle between one ground track and the middle of the region between the two adjacent
    # tracks in the Equator [rad]. This angle is measured from the Earth's center.
    θ = T(ΔT) * (ω_e - ∂Ω_∂t) / T(orbit_cycle) / 2
    sin_θ, cos_θ = sincos(θ)
    cot_θ = cos_θ / sin_θ

    # We need to compute the angle between one ground track and the middle of the region
    # between the two adjacent tracks measured from the Earth's center (β). Thus, first we
    # need to find the ground trace inclination [rad], which is a composition between the
    # Earth's rotation rate and the satellite speed.
    i_gt = ground_track_inclination(
        a,
        e,
        i;
        J2 = J2,
        J4 = J4,
        R0 = R0,
        m0 = m0,
        perturbation = perturbation,
        we = we
    )

    sin_i_gt, cos_i_gt = sincos(i_gt)
    cot_i_gt = cos_i_gt / sin_i_gt

    β = acot(cot_θ * sin_i_gt + cot_i_gt * cos_i_gt / sin_θ)

    # Compute the angle between the two ground tracks measured from the satellite. `a` is an
    # auxiliary distance and `γ` is the angle we are looking for.
    sin_β, cos_β = sincos(β)
    α = √(R₀^2 + a^2 - 2R₀ * a * cos_β)
    γ = asin(R₀ / α * sin_β)

    # Finally, the adjacent track distance is two times `γ`.
    return 2γ
end

"""
    ground_repeating_orbit_adjacent_track_distance(orbit_period::T1, i::T2, orbit_cycle::Integer; kwargs...) where {T1 <: Number, T2 <: Number} -> T

Compute the adjacent track distance [m] at Equator in a ground repeating orbit.  The orbit
is described by its orbital period `orbit_period` [s], inclination `i` [rad], and orbit
cycle `orbit_cycle` [day].

!!! note

    Internally, this function uses the precision obtained by promoting `T1` and `T2` to a
    float-pointing number `T`.

# Keywords

- `perturbation::Symbol`: Symbol to select the perturbation terms that will be used. It can
    be `:J0`, `:J2`, or `:J4`.
    (**Default**: `:J2`)
- `m0::Number`: Standard gravitational parameter for Earth [m³ / s²].
    (**Default**: `GM_EARTH`)
- `J2::Number`: J₂ perturbation term.
    (**Default**: `EGM_2008_J2`)
- `J4::Number`: J₄ perturbation term.
    (**Default**: `EGM_2008_J4`)
- `R0::Number`: Earth's equatorial radius [m].
    (**Default**: `EARTH_EQUATORIAL_RADIUS`)
- `we::Number`: Earth's angular speed [rad / s].
    (**Default**: `EARTH_ANGULAR_SPEED`)

# Extended help

A ground repeating orbit is any orbit that the number of revolutions per day is a rational
number. Hence, this type of orbit repeats its ground trace after a finite number of days.

The information `orbit_period` and `orbit_cyle` is redundant. However, they are necessary to
improve the algorithm precision. Otherwise, the `orbit_cycle` must be obtained by converting
the floating-point number `orbit_period` to a rational number, leading to numerical
problems.
"""
function ground_repeating_orbit_adjacent_track_distance(
    a::T1,
    e::T2,
    i::T3,
    orbit_cycle::Integer;
    perturbation::Symbol = :J2,
    J2::Number = EGM_2008_J2,
    J4::Number = EGM_2008_J4,
    R0::Number = EARTH_EQUATORIAL_RADIUS,
    m0::Number = GM_EARTH,
    we::Number = EARTH_ANGULAR_SPEED
) where {T1 <: Number, T2 <: Number, T3 <: Number}
    T   = float(promote_type(T1, T2, T3))
    R₀  = T(R0)
    ω_e = T(EARTH_ANGULAR_SPEED)

    # Compute the orbital period [s].
    ΔT = orbital_period(
        a,
        e,
        i;
        perturbation = perturbation,
        J2 = J2,
        J4 = J4,
        R0 = R0,
        m0 = m0
    )

    # Compute the RAAN time derivative [rad / s].
    ∂Ω_∂t = raan_time_derivative(
        a,
        e,
        i;
        perturbation = perturbation,
        J2 = J2,
        J4 = J4,
        m0 = m0,
        R0 = R0
    )

    # Angle between one ground track and the middle of the region between the two adjacent
    # tracks in the Equator [rad]. This angle is measured from the Earth's center.
    θ = T(ΔT) * (ω_e - ∂Ω_∂t) / T(orbit_cycle) / 2
    sin_θ, cos_θ = sincos(θ)
    cot_θ = cos_θ / sin_θ

    # We need to compute the angle between one ground track and the middle of the region
    # between the two adjacent tracks measured from the Earth's center (β). Thus, first we
    # need to find the ground trace inclination [rad], which is a composition between the
    # Earth's rotation rate and the satellite speed.
    i_gt = ground_track_inclination(
        a,
        e,
        i;
        J2 = J2,
        J4 = J4,
        R0 = R0,
        m0 = m0,
        perturbation = perturbation,
        we = we
    )

    sin_i_gt, cos_i_gt = sincos(i_gt)
    cot_i_gt = cos_i_gt / sin_i_gt

    β = acot(cot_θ * sin_i_gt + cot_i_gt * cos_i_gt / sin_θ)

    # Distance between two adjacent tracks on the Earth's surface.
    d = 2β * R₀

    return d
end
