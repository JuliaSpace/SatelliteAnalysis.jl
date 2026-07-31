## Description #############################################################################
#
# Miscellaneous functions for the SatelliteAnalysisDecayExt package.
#
############################################################################################

"""
    _classical_to_equinoctial(
        a::T,
        e::T,
        i::T,
        Ω::T,
        ω::T,
        M::T
    ) where T<:Number -> SVector{6, T}

Compute the equinoctial orbital elements from the classical orbital elements **[1]**:

- `a::T`: Semi-major axis [m].
- `e::T`: Eccentricity (0 ≤ e < 1) [-].
- `i::T`: Inclination [rad].
- `Ω::T`: Right Ascension of Ascending Node [rad].
- `ω::T`: Argument of Perigee [rad].
- `M::T`: Mean Anomaly [rad].

# Returns

- `SVector{6, T}`: Equinoctial orbital elements `[a, ψ, e_x, e_y, i_x, i_y]`.

# Extended help

## References

- **[1]** Battin, R. H. (1999). An Introduction to the Mathematics and Methods of
    Astrodynamics. Revised ed. AIAA Education Series, Reston, VA.
"""
function _classical_to_equinoctial(a::T, e::T, i::T, Ω::T, ω::T, M::T) where T<:Number
    ψ = M + Ω + ω
    ξ = Ω + ω

    sin_ξ, cos_ξ = sincos(ξ)
    sin_Ω, cos_Ω = sincos(Ω)
    sin_io2      = sin(i / 2)

    e_x = e * cos_ξ
    e_y = e * sin_ξ
    i_x = sin_io2 * cos_Ω
    i_y = sin_io2 * sin_Ω

    return @SVector [a, ψ, e_x, e_y, i_x, i_y]
end

"""
    _coe_to_rv(
        a::T,
        e::T,
        i::T,
        Ω::T,
        ω::T,
        M::T
    ) where T <: Number -> SVector{3, T}, SVector{3, T}

Compute position and velocity vectors from Classical Orbital Elements:

- `a::T`: Semi-major axis [m].
- `e::T`: Eccentricity (0 ≤ e < 1) [-].
- `i::T`: Inclination [rad].
- `Ω::T`: Right Ascension of Ascending Node [rad].
- `ω::T`: Argument of Perigee [rad].
- `M::T`: Mean Anomaly [rad].

# Returns

- `SVector{3, T}`: Position vector in ECI frame [m].
- `SVector{3, T}`: Velocity vector in ECI frame [m/s].
"""
function _coe_to_rv(a::T, e::T, i::T, Ω::T, ω::T, M::T) where T <: Number
    f  = mean_to_true_anomaly(e, M)
    ke = KeplerianElements(0.0, a, e, i, Ω, ω, f)

    return kepler_to_rv(ke)
end

"""
    _coe_to_rv_jacobian(e::T, i::T, Ω::T, ω::T) where T <: Number -> SMatrix{6, 6, T}

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

- `SMatrix{T, 6, 6}`: Jacobian matrix `J` of size `6x6`, where each entry corresponds to the
    partial derivative of an equinoctial element with respect to a classical orbital
    element.
"""
function _coe_to_rv_jacobian(e::T, i::T, Ω::T, ω::T) where T <: Number
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

    # E4 (e_y = e sin(ξ)) -> ₄₂, [4, 4], [4, 5]
    J₄₂ = sin_ξ                         # de_y/de
    J₄₄ = e * cos_ξ                     # de_y/dΩ
    J₄₅ = e * cos_ξ                     # de_y/dω

    # E5 (i_x = sin(i/2) cos(Ω)) -> [5, 3], [5, 4]
    J₅₃ = (1 // 2) * cos_io2 * cos_Ω    # dq/di
    J₅₄ = -sin_io2 * sin_Ω              # dq/dΩ

    # E6 (i_y = sin(i/2) sin(Ω)) -> [6, 3], [6, 4]
    J₆₃ = (1 // 2) * cos_io2 * sin_Ω    # dp/di
    J₆₄ = sin_io2 * cos_Ω               # dp/dΩ

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
    _equinoctial_to_classical(
        eq::AbstractVector{T}
    ) where T <: Number -> NTuple{6, T}

Conversion from Equinoctial Elements `equino` to Classical Orbital Elements (COE) **[1]**.
The input `equino` is a vector of equinoctial orbital elements `[a, ψ, e_x, e_y, i_x, i_y]`,
where:

- `a`    [m]   Semi-major axis
- `ψ`    [rad] Mean longitude
- `e_x`  [-]   Eccentricity vector component along cos(Ω+ω)
- `e_y`  [-]   Eccentricity vector component along sin(Ω+ω)
- `i_x`  [-]   Inclination half-angle component along cos(Ω)
- `i_y`  [-]   Inclination half-angle component along sin(Ω)

This routine ensures proper normalization of angles and handles potential singularities in
eccentricity and inclination.

!!! note

    - Eccentricity `e` is constrained to `e ≥ 1e-6` to avoid division by zero.
    - Inclination is reconstructed from half-angle parameters `(i_x, i_y)`.
    - All angular quantities are normalized to the range `[0, 2π)` [rad].

# Returns

- `T`: Semi-major axis [m].
- `T`: Eccentricity [-].
- `T`: Inclination [rad].
- `T`: Right ascension of ascending node [rad].
- `T`: Argument of perigee [rad].
- `T`: Mean anomaly [rad].

# Extended help

## References

- **[1]** Battin, R. H. (1999). An Introduction to the Mathematics and Methods of
    Astrodynamics. Revised ed. AIAA Education Series, Reston, VA.
"""
function _equinoctial_to_classical(
    eq::AbstractVector{T}
) where T <: Number
    # Unpack.
    a, ψ, e_x, e_y, i_x, i_y = eq

    e = max(√(e_x^2 + e_y^2), T(1e-6))

    i_half_sq = clamp(i_x^2 + i_y^2, 0, 1)
    i = 2asin(√i_half_sq)

    Ω = atan(i_y, i_x)
    ξ = atan(e_y, e_x)

    Ω = mod(Ω, 2π)
    ξ = mod(ξ, 2π)
    ω = mod(ξ - Ω, 2π)
    M = mod(ψ - ξ, 2π)

    a, e, i, Ω, ω, M = _normalize_classical_elements(a, e, i, Ω, ω, M)

    return a, e, i, Ω, ω, M
end

"""
    _mean_to_osculating_elements(
        a::T,
        e::T,
        i::T,
        Ω::T,
        ω::T,
        M::T
    ) where T <: Number -> NTuple{6, T}

Convert the mean classical elements `[a, e, i, Ω, ω, M]` to osculating elements under the J2
perturbation model **[1]**.

# Arguments

- `a::T`: Semi-major axis [m].
- `e::T`: Eccentricity (0 ≤ e < 1) [-].
- `i::T`: Inclination [rad].
- `Ω::T`: Right ascension of ascending node [rad].
- `ω::T`: Argument of perigee [rad].
- `M::T`: Mean anomaly [rad].

# Returns

- `T`: Osculating semi-major axis [m].
- `T`: Osculating eccentricity [-].
- `T`: Osculating inclination [rad].
- `T`: Osculating right ascension of ascending node [rad].
- `T`: Osculating argument of perigee [rad].
- `T`: Osculating mean anomaly [rad].

# Extended help

# References

- **[1]** Vallado, D. A. (2013). *Fundamentals of Astrodynamics and Applications*. 4th ed.
    Microcosm Press, Hawthorne, CA.
"""
function _mean_to_osculating_elements(a::T, e::T, i::T, Ω::T, ω::T, M::T) where T <: Number
    # Normalize.
    a, e, i, Ω, ω, M = _normalize_classical_elements(a, e, i, Ω, ω, M)
    e = max(e, 1e-6)
    i = max(i, 1e-6)

    # Convert to osculating.
    orb_tod      = KeplerianElements(0.0, a, e, i, Ω, ω, mean_to_true_anomaly(e, M))
    orbp         = Propagators.init(Val(:J2osc), orb_tod)
    r_tod, v_tod = Propagators.propagate!(orbp, 0.0)
    ke           = rv_to_kepler(r_tod, v_tod)

    ap = ke.a
    ep = ke.e
    ip = ke.i
    Ωp = ke.Ω
    ωp = ke.ω
    Mp = true_to_mean_anomaly(ke.e, ke.f)

    # Return new state.
    return ap, ep, ip, Ωp, ωp, Mp
end

"""
    _osculating_to_mean_elements(
        a::T,
        e::T,
        i::T,
        Ω::T,
        ω::T,
        M::T
    ) where T <: Number -> NTuple{6, T}

Convert the osculating classical elements `[a, e, i, Ω, ω, M]` to mean elements under the J2
perturbation model **[1]**.

# Arguments

- `a::T`: Semi-major axis [m].
- `e::T`: Eccentricity (0 ≤ e < 1) [-].
- `i::T`: Inclination [rad].
- `Ω::T`: Right ascension of ascending node [rad].
- `ω::T`: Argument of perigee [rad].
- `M::T`: Mean anomaly [rad].

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
function _osculating_to_mean_elements(a::T, e::T, i::T, Ω::T, ω::T, M::T) where T <: Number
    # Normalize.
    a, e, i, Ω, ω, M = _normalize_classical_elements(a, e, i, Ω, ω, M)
    e = max(e, 1e-6)
    i = max(i, 1e-6)

    # Convert to mean.
    orb_osc_tod  = KeplerianElements(0.0, a, e, i, Ω, ω, mean_to_true_anomaly(e, M))
    r_tod, v_tod = kepler_to_rv(orb_osc_tod)
    ke, _        = fit_j2osc_mean_elements([0], [r_tod], [v_tod]; max_iterations = 50, verbose = false)

    ap = ke.a
    ep = ke.e
    ip = ke.i
    Ωp = ke.Ω
    ωp = ke.ω
    Mp = true_to_mean_anomaly(ke.e, ke.f)

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
    a::T,
    e::T,
    i::T,
    Ω::T,
    ω::T,
    M::T
) where T <: Number
    a_norm = a > 0 ? a : abs(a)
    e_norm = e < 0 ? 0.0 : e
    i_norm = clamp(i, 0.0, π)
    Ω_norm = mod(Ω, 2π)
    ω_norm = mod(ω, 2π)
    M_norm = mod(M, 2π)

    return a_norm, e_norm, i_norm, Ω_norm, ω_norm, M_norm
end

"""
    _r_eci_to_hill(
        r_eci::AbstractVector{T},
        v_eci::AbstractVector{T}
    ) where T <: Number -> SMatrix{3, 3, T}

Compute the Direction Cosine Matrix (DCM) from an Earth Centered Inertial (ECI) frame to the
Hill frame (also known as then RTN or RSW frame), given the satellite position [m] and
satellite velocity [m/s] in the eci reference frame.

The Hill frame is defined as follows:

- X-axis: along the radial direction (from Earth to satellite).
- Y-axis: along the along-track direction (in the direction of motion).
- Z-axis: along the cross-track direction (perpendicular to the orbital plane).
"""
function _r_eci_to_hill(
    r_eci::AbstractVector{T},
    v_eci::AbstractVector{T}
) where T <: Number
    r̄_eci = normalize(r_eci)
    h_eci = cross(r_eci, v_eci)
    h̄_eci = normalize(h_eci)
    θ̄_eci = cross(h̄_eci, r̄_eci)

    return @SMatrix [
        r̄_eci[1] r̄_eci[2] r̄_eci[3]
        θ̄_eci[1] θ̄_eci[2] θ̄_eci[3]
        h̄_eci[1] h̄_eci[2] h̄_eci[3]
    ]
end

