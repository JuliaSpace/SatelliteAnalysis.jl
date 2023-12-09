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
    ground_repeating_orbit_adjacent_track_angle(h::T1, orbit_period::T2, i::T3, orbit_cycle::Integer) where {T1 <: Number, T2 <: Number, T3 <: Number} -> T

Compute the adjacent track angle [rad] at Equator in a ground repeating orbit measured from
the satellite position. The orbit is described by its altitude at Equator `h` [m], orbital
period `orbit_period` [s], inclination `i` [rad], and orbit cycle `orbit_cyle` [day].

!!! warning
    The code does not check if the altitude `h` is correct for the orbit with the period
    `orbit_period`.

!!! note
    Internally, this function uses the precision obtained by promoting `T1` and `T2` to a
    float-pointing number `T`.

# Extended help

A ground repeating orbit is any orbit that the number of revolutions per day is a rational
number. Hence, this type of orbit repeats its ground trace after a finite number of days.

The information `orbit_period` and `orbit_cyle` is redundant. However, they are necessary to
improve the algorithm precision. Otherwise, the `orbit_cycle` must be obtained by converting
the floating-point number `orbit_period` to a rational number, leading to numerical
problems.
"""
function ground_repeating_orbit_adjacent_track_angle(
    h::T1,
    orbit_period::T2,
    i::T3,
    orbit_cycle::Integer
) where {T1 <: Number, T2 <: Number, T3 <: Number}
    T   = float(promote_type(T1, T2, T3))
    R₀  = T(EARTH_EQUATORIAL_RADIUS)
    ω_e = T(EARTH_ANGULAR_SPEED)

    # Satellite mean angular velocity [rad / s].
    ω_s = T(2π) / T(orbit_period)

    # Angle between one ground track and the middle of the region between the two adjacent
    # tracks in the Equator [rad]. This angle is measured from the Earth's center.
    θ = T(orbit_period) / T(orbit_cycle) * T(π) / 43200 / 2
    sin_θ, cos_θ = sincos(θ)
    cot_θ = cos_θ / sin_θ

    # We need to compute the angle between one ground track and the middle of the region
    # between the two adjacent tracks measured from the Earth's center (β). Thus, first we
    # need to find the ground trace inclination, which is a composition between the Earth's
    # rotation rate and the satellite speed.
    sin_i, cos_i = sincos(i)
    cot_i = cos_i / sin_i

    i_gt = atan(ω_s * sin_i, ω_s * cos_i - ω_e)
    β    = acot(cot_θ * sin_i + cot_i * cos_i / sin_θ)

    # Compute the angle between the two ground tracks measured from the satellite. `a` is an
    # auxiliary distance and `γ` is the angle we are looking for.
    sin_β, cos_β = sincos(β)
    α = √(R₀^2 + (R₀ + T(h))^2 - 2R₀ * (R₀ + T(h)) * cos_β)
    γ = asin(R₀ / α * sin_β)

    # Finally, the adjacent track distance is two times `γ`.
    return 2γ
end

"""
    ground_repeating_orbit_adjacent_track_distance(orbit_period::T1, i::T2, orbit_cycle::Integer) where {T1 <: Number, T2 <: Number} -> T

Compute the adjacent track distance [m] at Equator in a ground repeating orbit.  The orbit
is described by its orbital period `orbit_period` [s], inclination `i` [rad], and orbit
cycle `orbit_cycle` [day].

!!! note
    Internally, this function uses the precision obtained by promoting `T1` and `T2` to a
    float-pointing number `T`.

# Extended help

A ground repeating orbit is any orbit that the number of revolutions per day is a rational
number. Hence, this type of orbit repeats its ground trace after a finite number of days.

The information `orbit_period` and `orbit_cyle` is redundant. However, they are necessary to
improve the algorithm precision. Otherwise, the `orbit_cycle` must be obtained by converting
the floating-point number `orbit_period` to a rational number, leading to numerical
problems.
"""
function ground_repeating_orbit_adjacent_track_distance(
    orbit_period::T1,
    i::T2,
    orbit_cycle::Integer
) where {T1 <: Number, T2 <: Number}
    T   = float(promote_type(T1, T2))
    R₀  = T(EARTH_EQUATORIAL_RADIUS)
    ω_e = T(EARTH_ANGULAR_SPEED)

    # Satellite mean angular velocity [rad / s].
    ω_s = T(2π) / T(orbit_period)

    # Angle between one ground track and the middle of the region between the two adjacent
    # tracks in the Equator [rad]. This angle is measured from the Earth's center.
    θ = T(orbit_period) / T(orbit_cycle) * T(π) / 43200 / 2
    sin_θ, cos_θ = sincos(θ)
    cot_θ = cos_θ / sin_θ

    # We need to compute the angle between one ground track and the middle of the region
    # between the two adjacent tracks measured from the Earth's center (β). Thus, first we
    # need to find the ground trace inclination, which is a composition between the Earth's
    # rotation rate and the satellite speed.
    sin_i, cos_i = sincos(i)
    cot_i = cos_i / sin_i

    i_gt = atan(ω_s * sin_i, ω_s * cos_i - ω_e)
    β    = acot(cot_θ * sin_i + cot_i * cos_i / sin_θ)

    # Distance between two adjacent tracks on the Earth's surface.
    d = 2β * R₀

    return d
end
