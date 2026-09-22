## Description #############################################################################
#
# Miscellaneous functions for the SatelliteAnalysisDecayExt package.
#
############################################################################################

"""
    _state_to_alternate_equinoctial(u::AbstractVector{T}, epoch::Number) where {T <: Number} -> AlternateEquinoctialElements

Convert the state vector `u = [a, ψ, e_x, e_y, i_x, i_y]` of the numerical integration to the
alternate equinoctial elements of **SatelliteToolboxBase.jl** with the Julian Day `epoch`
[UTC]. The elements are the same, but stored in a different order: `h = e_y`, `k = e_x`,
`p = i_y`, `q = i_x`, and `mean_longitude = ψ`.

The adaptive integrator can evaluate trial stages with unphysical states. If the
inclination elements satisfy `i_x² + i_y² > 1`, which does not represent an orbit, they are
scaled to the unit circle so that the conversion to Keplerian elements remains finite,
allowing the error control to reject the step.
"""
function _state_to_alternate_equinoctial(
    u::AbstractVector{T}, epoch::Number
) where {T <: Number}
    a, ψ, e_x, e_y, i_x, i_y = u

    sin_io2 = hypot(i_x, i_y)

    if sin_io2 > 1
        i_x /= sin_io2
        i_y /= sin_io2
    end

    return AlternateEquinoctialElements(epoch, a, e_y, e_x, i_y, i_x, ψ)
end

"""
    _state_to_keplerian(u::AbstractVector{T}, epoch::Number) where {T <: Number} -> KeplerianElements{MeanAnomaly}

Convert the state vector `u = [a, ψ, e_x, e_y, i_x, i_y]` of the numerical integration to
Keplerian elements storing the mean anomaly with the Julian Day `epoch` [UTC]. The angles
are returned in the interval `[0, 2π)` [rad], and the eccentricity is not clamped, so the
returned value is faithful to the input. Consumers that require `e > 0` must clamp it
themselves.
"""
function _state_to_keplerian(u::AbstractVector{T}, epoch::Number) where {T <: Number}
    aee = _state_to_alternate_equinoctial(u, epoch)
    return convert(KeplerianElements{MeanAnomaly}, aee)
end

"""
    _keplerian_to_state(ke::KeplerianElements{Tanomaly, Tepoch, T}) where {Tanomaly, Tepoch, T} -> SVector{6, T}

Convert the Keplerian elements `ke`, with any anomaly type, to the state vector
`u = [a, ψ, e_x, e_y, i_x, i_y]` of the numerical integration (see
[`_state_to_alternate_equinoctial`](@ref)).
"""
function _keplerian_to_state(
    ke::KeplerianElements{Tanomaly, Tepoch, T}
) where {Tanomaly <: AbstractAnomaly, Tepoch <: Number, T <: Number}
    aee = convert(AlternateEquinoctialElements, ke)

    return SVector{6, T}(
        aee.semi_major_axis, aee.mean_longitude, aee.k, aee.h, aee.q, aee.p
    )
end

"""
    _state_apsis_altitudes(u::AbstractVector{T}) where {T <: Number} -> T, T

Compute the mean perigee and apogee altitudes [m] above the Earth's equatorial radius from
the state vector `u = [a, ψ, e_x, e_y, i_x, i_y]` of the numerical integration. Only the
semi-major axis and the eccentricity are required. Hence, this function avoids the full
conversion to Keplerian elements in the callbacks of the integration.
"""
function _state_apsis_altitudes(u::AbstractVector{T}) where {T <: Number}
    a = u[1]
    e = hypot(u[3], u[4])

    perigee_altitude = a * (1 - e) - T(EARTH_EQUATORIAL_RADIUS)
    apogee_altitude  = a * (1 + e) - T(EARTH_EQUATORIAL_RADIUS)

    return perigee_altitude, apogee_altitude
end

"""
    _classical_to_equinoctial_jacobian(e::T, i::T, Ω::T, ω::T) where T <: Number -> SMatrix{6, 6, T}

Compute the Jacobian matrix of the transformation from Classical Orbital Elements (COE) to
Equinoctial Elements.

This routine computes the Jacobian matrix associated with the mapping between classical
orbital elements and equinoctial elements. The Jacobian provides the partial derivatives of
equinoctial elements with respect to classical orbital elements, which is useful for
sensitivity analysis, covariance propagation, and orbit determination.

# Arguments

- `e::Number`: Eccentricity (0 ≤ e < 1) [-].
- `i::Number`: Inclination [rad].
- `Ω::Number`: Right ascension of ascending node [rad].
- `ω::Number`: Argument of perigee [rad].

# Returns

- `SMatrix{6, 6, T}`: Jacobian matrix `J` of size `6x6`, where each entry corresponds to the
    partial derivative of an equinoctial element with respect to a classical orbital
    element.
"""
function _classical_to_equinoctial_jacobian(e::T, i::T, Ω::T, ω::T) where {T <: Number}
    ξ = Ω + ω

    sin_ξ, cos_ξ = sincos(ξ)
    sin_Ω, cos_Ω = sincos(Ω)
    cos_io2      = cos(i / 2)
    sin_io2      = sin(i / 2)

    # E1 (a) -> [1, 1]
    J₁₁ = T(1)

    # E2 (Ψ = M + Ω + ω) -> [2, 4], [2, 5], [2, 6]
    J₂₄ = T(1)                          # dΨ/dΩ
    J₂₅ = T(1)                          # dΨ/dω
    J₂₆ = T(1)                          # dΨ/dM

    # E3 (e_x = e cos(ξ)) -> [3, 2], [3, 4], [3, 5]
    J₃₂ = cos_ξ                         # de_x/de
    J₃₄ = -e * sin_ξ                    # de_x/dΩ
    J₃₅ = -e * sin_ξ                    # de_x/dω

    # E4 (e_y = e sin(ξ)) -> [4, 2], [4, 4], [4, 5]
    J₄₂ = sin_ξ                         # de_y/de
    J₄₄ = e * cos_ξ                     # de_y/dΩ
    J₄₅ = e * cos_ξ                     # de_y/dω

    # E5 (i_x = sin(i/2) cos(Ω)) -> [5, 3], [5, 4]
    J₅₃ = (1 // 2) * cos_io2 * cos_Ω    # di_x/di
    J₅₄ = -sin_io2 * sin_Ω              # di_x/dΩ

    # E6 (i_y = sin(i/2) sin(Ω)) -> [6, 3], [6, 4]
    J₆₃ = (1 // 2) * cos_io2 * sin_Ω    # di_y/di
    J₆₄ = sin_io2 * cos_Ω               # di_y/dΩ

    J = @SMatrix T[
        J₁₁   0    0    0   0   0
         0    0    0   J₂₄ J₂₅ J₂₆
         0   J₃₂   0   J₃₄ J₃₅  0
         0   J₄₂   0   J₄₄ J₄₅  0
         0    0   J₅₃  J₅₄  0   0
         0    0   J₆₃  J₆₄  0   0
    ]

    return J
end

"""
    _mean_to_osculating_rv(
        a::T,
        e::T,
        i::T,
        Ω::T,
        ω::T,
        f::T,
        orbp::OrbitPropagatorJ2Osculating
    ) where T <: Number -> SVector{3, T}, SVector{3, T}

Compute the osculating position and velocity vectors from the mean classical elements
`[a, e, i, Ω, ω, f]` under the J2 perturbation model **[1]** using the pre-allocated J2
osculating propagator `orbp`, which is re-initialized in place.

!!! note

    This function is called at every sampling point of the right-hand side. Hence, it does
    not normalize the inputs, which must be performed by the caller: the eccentricity and
    the inclination must not be lower than `1e-6`, and the angles must be in the interval
    `[0, 2π)` [rad].

# Arguments

- `a::T`: Mean semi-major axis [m].
- `e::T`: Mean eccentricity (0 ≤ e < 1) [-].
- `i::T`: Mean inclination [rad].
- `Ω::T`: Mean right ascension of ascending node [rad].
- `ω::T`: Mean argument of perigee [rad].
- `f::T`: Mean true anomaly [rad].
- `orbp::OrbitPropagatorJ2Osculating`: Pre-allocated J2 osculating propagator. It is
    re-initialized in place, avoiding one propagator allocation per call.

# Returns

- `SVector{3, T}`: Osculating position vector [m] in the frame of the input elements.
- `SVector{3, T}`: Osculating velocity vector [m/s] in the frame of the input elements.

# Extended help

# References

- **[1]** Vallado, D. A. (2013). *Fundamentals of Astrodynamics and Applications*. 4th ed.
    Microcosm Press, Hawthorne, CA.
"""
function _mean_to_osculating_rv(
    a::T, e::T, i::T, Ω::T, ω::T, f::T, orbp::OrbitPropagatorJ2Osculating
) where {T <: Number}
    orb_tod = KeplerianElements(0.0, a, e, i, Ω, ω, f)

    # Convert to osculating. The J2 osculating conversion is not total: for unphysical
    # states evaluated by the integrator in trial stages (mainly with loose tolerances),
    # the short-period corrections can push the osculating eccentricity outside [0, 1),
    # which throws inside the propagator. In this case, fall back to the mean elements so
    # the right-hand side remains finite and the step error control can act.
    r_tod, v_tod = try
        Propagators.init!(orbp, orb_tod)
        Propagators.propagate!(orbp, 0.0)
    catch err
        # The propagator and the element conversions throw `ArgumentError` or
        # `DomainError` for those unphysical states. Any other exception is a genuine
        # error and must propagate.
        err isa Union{ArgumentError, DomainError} || rethrow()
        kepler_to_rv(orb_tod)
    end

    return r_tod, v_tod
end

"""
    _osculating_to_mean_elements(
        a::T,
        e::T,
        i::T,
        Ω::T,
        ω::T,
        f::T
    ) where T <: Number -> NTuple{6, T}

Convert the osculating classical elements `[a, e, i, Ω, ω, f]` to the mean elements
`[a, e, i, Ω, ω, M]` under the J2 perturbation model **[1]**. Notice that the input contains
the true anomaly, as in the elements provided by the user, whereas the output contains the
mean anomaly, as in the state of the numerical integration, avoiding unnecessary
conversions between the anomalies.

# Arguments

- `a::T`: Semi-major axis [m].
- `e::T`: Eccentricity (0 ≤ e < 1) [-].
- `i::T`: Inclination [rad].
- `Ω::T`: Right ascension of ascending node [rad].
- `ω::T`: Argument of perigee [rad].
- `f::T`: True anomaly [rad].

# Returns

- `T`: Mean semi-major axis [m].
- `T`: Mean eccentricity [-].
- `T`: Mean inclination [rad].
- `T`: Mean right ascension of ascending node [rad].
- `T`: Mean argument of perigee [rad].
- `T`: Mean mean anomaly [rad].

# Extended help

# References

- **[1]** Vallado, D. A. (2013). *Fundamentals of Astrodynamics and Applications*. 4th ed.
    Microcosm Press, Hawthorne, CA.
"""
function _osculating_to_mean_elements(
    a::T, e::T, i::T, Ω::T, ω::T, f::T
) where {T <: Number}
    # Normalize.
    a, e, i, Ω, ω, f = _normalize_classical_elements(a, e, i, Ω, ω, f)
    e = max(e, T(1e-6))
    i = max(i, T(1e-6))

    # Convert to mean.
    orb_osc_tod  = KeplerianElements(0.0, a, e, i, Ω, ω, f)
    r_tod, v_tod = kepler_to_rv(orb_osc_tod)

    ke, _ = fit_j2osc_mean_elements(
        [0.0], [r_tod], [v_tod]; max_iterations = 50, verbose = false
    )

    ap = ke.a
    ep = ke.e
    ip = ke.i
    Ωp = ke.Ω
    ωp = ke.ω

    # The fitted elements already store the mean anomaly. Hence, this function does not
    # perform any conversion in this case.
    Mp = mean_anomaly(ke)

    # Return new state.
    return ap, ep, ip, Ωp, ωp, Mp
end

"""
    _normalize_classical_elements(
        a::T,
        e::T,
        i::T,
        Ω::T,
        ω::T,
        M::T
    ) where T <: Number -> NTuple{6, T}

Normalize classical orbital elements to valid ranges:

- Semi-major axis `a` is forced positive.
- Eccentricity `e` is clamped to non-negative values.
- Inclination `i` is clamped to [0, π] [rad].
- Angular elements (Ω, ω, M) are normalized to [0, 2π) [rad].

# Arguments

- `a::T`: Semi-major axis [m].
- `e::T`: Eccentricity [-].
- `i::T`: Inclination [rad].
- `Ω::T`: Right Ascension of Ascending Node [rad].
- `ω::T`: Argument of Perigee [rad].
- `M::T`: Mean Anomaly [rad].

# Returns

- `T`: Normalized semi-major axis [m].
- `T`: Normalized eccentricity [-].
- `T`: Normalized inclination [rad].
- `T`: Normalized right ascension of ascending node [rad].
- `T`: Normalized argument of perigee [rad].
- `T`: Normalized mean anomaly [rad].
"""
function _normalize_classical_elements(
    a::T, e::T, i::T, Ω::T, ω::T, M::T
) where {T <: Number}
    a_norm = abs(a)
    e_norm = max(e, zero(T))
    i_norm = clamp(i, zero(T), T(π))
    Ω_norm = mod(Ω, 2π)
    ω_norm = mod(ω, 2π)
    M_norm = mod(M, 2π)

    return a_norm, e_norm, i_norm, Ω_norm, ω_norm, M_norm
end
