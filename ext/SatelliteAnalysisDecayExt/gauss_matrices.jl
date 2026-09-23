## Description #############################################################################
#
# Functions to compute the Gauss variational matrix A(y) for orbital element propagation
# under perturbations.
#
############################################################################################

"""
    _equinoctial_gauss_variational_matrices(
        a::T,
        e::T,
        i::T,
        Ω::T,
        ω::T,
        f::T,
        r::T,
        p::T,
        h::T,
        η::T
    ) where T <: Number -> SMatrix{6, 3, T}

Compute the Gauss variational equations in matrix form for the alternate equinoctial
orbital elements **[1]**:

    u̇ = A(u) * ap + B

where `u = [a, h, k, p, q, λ]` are the alternate equinoctial elements, stored in the same
order as the fields of `AlternateEquinoctialElements`, and `ap = [u_r, u_θ, u_h]` is the
perturbation acceleration represented in the Hill frame (radial, along-track, cross-track).
The constant Kepler term `B = [0, 0, 0, 0, 0, n]`, where `n` is the mean motion, is not
returned and must be added once by the caller.

!!! note

    Inside this function, the eccentricity elements are named `e_x = k = e cos(Ω + ω)` and
    `e_y = h = e sin(Ω + ω)` because `h` denotes the specific angular momentum.

This formulation is obtained by analytically combining the classical Gauss variational
equations with the Jacobian of the classical-to-equinoctial transformation. All `1 / e`
terms cancel symbolically, so the matrix is well-conditioned for circular orbits. The
remaining singularity at `i = π` (retrograde equatorial orbits) is inherent to this
equinoctial element set.

# Arguments

- `a::T`: Semi-major axis [m].
- `e::T`: Eccentricity [-].
- `i::T`: Inclination [rad].
- `Ω::T`: Right ascension of the ascending node [rad].
- `ω::T`: Argument of perigee [rad].
- `f::T`: True anomaly [rad].
- `r::T`: Orbit radius at the true anomaly `f` [m].
- `p::T`: Semi-latus rectum [m].
- `h::T`: Specific angular momentum [m²/s].
- `η::T`: `√(1 - e²)` [-].

# Returns

- `SMatrix{6, 3, T}`: Matrix `A` multiplying the perturbation vector `[u_r, u_θ, u_h]`.

# References

- **[1]** Battin, R. H. (1999). An Introduction to the Mathematics and Methods of
    Astrodynamics. Revised ed. AIAA Education Series, Reston, VA.
"""
function _equinoctial_gauss_variational_matrices(
    a::T, e::T, i::T, Ω::T, ω::T, f::T, r::T, p::T, h::T, η::T
) where {T <: Number}
    ξ = Ω + ω
    L = f + ξ
    θ = f + ω

    sin_f, cos_f     = sincos(f)
    sin_L, cos_L     = sincos(L)
    sin_θ, cos_θ     = sincos(θ)
    sin_Ω, cos_Ω     = sincos(Ω)
    sin_ξ, cos_ξ     = sincos(ξ)
    sin_io2, cos_io2 = sincos(i / 2)
    tan_io2          = sin_io2 / cos_io2

    e_x = e * cos_ξ
    e_y = e * sin_ξ

    a²  = a^2
    psr = p + r
    k₁  = 2a² / h
    k₂  = r * sin_θ / h
    k₃  = 1 / (1 + η)

    A = @SMatrix T[
        (k₁ * e * sin_f)                      (k₁ * p / r)                   0
        (-p * cos_L / h)                      ((psr * sin_L + r * e_y) / h)  (+k₂ * tan_io2 * e_x)
        (p * sin_L / h)                       ((psr * cos_L + r * e_x) / h)  (-k₂ * tan_io2 * e_y)
        0  0  (r / h * (cos_io2 * sin_Ω * cos_θ + cos_Ω * sin_θ / cos_io2) / 2)
        0  0  (r / h * (cos_io2 * cos_Ω * cos_θ - sin_Ω * sin_θ / cos_io2) / 2)
        (-(p * e * cos_f * k₃ + 2r * η) / h)  (psr * e * sin_f * k₃ / h)     (k₂ * tan_io2)
    ]

    return A
end
