## Description #############################################################################
#
# Functions to compute the Gauss variational matrices A(y) and B for orbital element
# propagation under perturbations.
#
############################################################################################

"""
    _gauss_variational_matrices(
        a::T,
        e::T,
        i::T,
        ω::T,
        f::T,
        r::T,
        p::T,
        h::T,
        η::T,
        n::T
    ) where T <: Number -> Tuple{SMatrix{6, 3, T}, SVector{6, T}}

Compute Gauss variational matrices A(y) and B **[1]**.

This routine builds the matrix form of Gauss' variational equations:

    ẏ = A(y) * p + B

where `y = [a, e, i, Ω, ω, M]`, `p = [u_r, u_θ, u_h]`.

# Keywords

- `a::Number`: Semi-major axis [m].
- `e::Number`: Stabilized eccentricity (≥ 1e-6) [-].
- `i::Number`: Inclination [rad].
- `ω::Number`: Argument of perigee [rad].
- `f::Number`: True anomaly [rad].
- `r::Number`: Radius [m].
- `p::Number`: Semi-latus rectum [m].
- `h::Number`: Specific angular momentum [m²/s].
- `η::Number`: √(1 - e²) [-].
- `n::Number`: Mean motion [rad/s].

# Returns

- `A::Matrix{Float64}`: 6x3 matrix multiplying perturbation vector `[u_r, u_θ, u_h]`.
- `B::Vector{Float64}`: Constant vector `[0, 0, 0, 0, 0, n]`.

# Extended help

# References

- **[1]** Battin, R. H. (1999). An Introduction to the Mathematics and Methods of
    Astrodynamics. Revised ed. AIAA Education Series, Reston, VA.
"""
function _gauss_variational_matrices(
    a::T,
    e::T,
    i::T,
    ω::T,
    f::T,
    r::T,
    p::T,
    h::T,
    η::T,
    n::T
) where T <: Number
    sin_f, cos_f = sincos(f)
    sin_i, cos_i = sincos(i)

    θ  = f + ω
    sin_θ, cos_θ = sincos(θ)

    a²      = a^2
    h_sin_i = h * sin_i
    r_sin_θ = r * sin_θ
    k₁      = 2a² / h
    he      = h * e
    re      = r * e
    η_o_he  = η / he
    psr     = p + r

    A = @SMatrix T[
        (k₁ * e * sin_f)              (k₁ * p / r)              0
        (p / h * sin_f)               ((psr * cos_f + re) / h)  0
        0                             0                         (r * cos_θ / h)
        0                             0                         (r_sin_θ / h_sin_i)
        (-p / he * cos_f)             (psr / he * sin_f)      (-r_sin_θ * cos_i / h_sin_i)
        (η_o_he * (p * cos_f - 2re))  (-η_o_he * psr * sin_f)   0
    ]

    B = T[0, 0, 0, 0, 0, n]

    return A, B
end
