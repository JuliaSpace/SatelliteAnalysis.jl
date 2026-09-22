## Description #############################################################################
#
# Functions to design frozen orbits.
#
## References ##############################################################################
#
# [1] Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications. 4th ed.
#     Microcosm Press, Hawthorne, CA.
#
# [2] Kaula, W. M (1966). Theory of Satellite Geodesy. Blaisdell Publishing Company,
#     Waltham, MA.
#
# [3] Rosborough, G. W.; Ocampo, C. A (1991). Influence of higher degree zonals on the
#     frozen orbit geometry. Proceedings of the AAS/AIAA Astrodynamics Conference, Durango,
#     CO.
#
############################################################################################

export frozen_orbit

# Cache with the default gravity model (EGM96), avoiding fetching and parsing the ICGEM file
# at every call. The lock makes the lazy initialization thread-safe.
const _FROZEN_ORBIT_DEFAULT_GRAVITY_MODEL = Ref{Union{Nothing, AbstractGravityModel}}(
    nothing
)

const _FROZEN_ORBIT_DEFAULT_GRAVITY_MODEL_LOCK = ReentrantLock()

"""
    frozen_orbit(a::Number, i::Number; kwargs...) -> Float64, Float64

Compute the eccentricity [ ] and argument of perigee [rad] to obtain a frozen orbit when the
orbit has semi-major axis `a` [m] and inclination `i` [rad]. This function uses the theory
in **[1]**.

!!! note

    This function uses `BigFloat` internally to perform all computations, allowing very high
    degrees. However, the user must ensure that the default precision is enough for the
    required degree. Refer to the function `setprecision` for more information.

# Keywords

- `gravity_model::Union{Nothing, AbstractGravityModel}`: Gravity model used to compute the
    frozen eccentricity. Refer to the object `AbstractGravityModel` of the package
    `SatelliteToolboxGravityModels.jl` for more information. If it is `nothing`, the system
    will automatically fetch and load the EGM96 gravity model at the first call, keeping it
    in memory for the next ones.
    (**Default** = `nothing`)
- `max_degree::Int`: Maximum gravity model degree used to compute the frozen eccentricity.
    If it is equal to or lower than 0, the maximum degree in `gravity_model` will be used.
    Otherwise, if it is lower than 3 or higher than the `gravity_model` maximum degree, it
    will be clamped accordingly.
    (**Default** = 53)

# References

- **[1]** Rosborough, G. W.; Ocampo, C. A (1991). Influence of higher degree zonals on the
    frozen orbit geometry. Proceedings of the AAS/AIAA Astrodynamics Conference, Durango,
    CO.

# Extended Help

Due to the Earth's gravitational perturbation, the orbit of a satellite will experience
secular changes in the argument of perigee. Hence, the satellite mean altitude per latitude
will differ during the mission. This effect can be problematic, especially if we must
compare images by a camera onboard the satellite in different periods. The altitude
variation will change the resolution, leading to some problems when comparing the data.

We can avoid this problem if we compute an eccentricity and the argument of perigee that
yields theoretically:

    de      dω
    ── = 0, ── = 0
    dt      dt

This orbit is called **frozen**. Refer to **[1]** for more information.

## Throws

- `ArgumentError`: If the inclination `i` is not within the interval `(0, π)` [rad] because
    the frozen orbit is not defined for equatorial orbits.

## Examples

```julia-repl
julia> using SatelliteAnalysis

julia> frozen_orbit(7130.982e3, 98.410 |> deg2rad)
(0.0011641853028456078, 1.5707963267948966)

julia> jgm3 = GravityModels.load(IcgemFile, fetch_icgem_file(:JGM3))
[ Info: Downloading the ICGEM file 'JGM3.gfc' from 'http://icgem.gfz-potsdam.de/getmodel/gfc/a3375e01a717ac162962138a5e94f10
466b71aa4a130d7f7d5b18ab3d5f90c3d/JGM3.gfc'...
IcgemFile{Float64}:
      Product type : gravity_field
       Model name  : JGM3
  Gravity constant : 3.986004415e14
            Radius : 6.3781363e6
    Maximum degree : 70
            Errors : formal
       Tide system : unknown
              Norm : fully_normalized
         Data type : Float64

julia> frozen_orbit(7130.982e3, 98.410 |> deg2rad; gravity_model = jgm3)
(0.001163484769069545, 1.5707963267948966)
```
"""
function frozen_orbit(
    a::Number,
    i::Number;
    gravity_model::Union{Nothing, AbstractGravityModel} = nothing,
    max_degree::Int = 53,
)
    # The inclination functions are singular in equatorial orbits, in which the argument of
    # perigee is not defined.
    ((0 < i < π) && !iszero(sin(i))) ||
        throw(ArgumentError("The inclination must be within the interval (0, π) [rad]."))

    # Fetch EGM-96 gravity model if the user does not specify one.
    gm = isnothing(gravity_model) ? _frozen_orbit__default_gravity_model() : gravity_model

    gm_max_degree = GravityModels.maximum_degree(gm)

    # If `max_degree` is negative, use maximum degree available.
    if max_degree <= 0
        max_degree = gm_max_degree
    end

    # The maximum degree must be at least 3.
    max_degree = min(gm_max_degree, max(3, max_degree))

    # Obtain the maximum `p` given the maximum degree.
    p_max = floor(Int, (max_degree - 1) / 2)

    # The theory in [2] requires unnormalized coefficients. Hence, we need the normalization
    # of the coefficients in the gravity model to convert them.
    coefficient_norm = GravityModels.coefficient_norm(gm)

    # We must using `BigInt` and `BigFloat` to compute for high degree.
    ab  = big(a)
    ib  = big(i)
    num = zero(BigFloat)
    den = zero(BigFloat)
    R_e = big(GravityModels.radius(gm))

    sin_i, cos_i = sincos(ib)

    # Auxiliary variables to compute (R_e / a)^2p.
    fact = R_e / ab
    prod = fact * fact # ............................................... (R_e / a)^2 (p = 0)

    # Compute the numerator and denominator of the frozen eccentricity as indicated in Eq.
    # 38 [3, p. 1297].
    for p in 1:p_max
        C_2p_0  = -(GravityModels.coefficients(gm, 2p, 0) |> first) |> big
        C_2p1_0 = -(GravityModels.coefficients(gm, 2p + 1, 0) |> first) |> big

        C_2p_0  *= _frozen_orbit__zonal_unnormalization_factor(coefficient_norm, 2p)
        C_2p1_0 *= _frozen_orbit__zonal_unnormalization_factor(coefficient_norm, 2p + 1)

        F_2p_0_p, ∂F_2p_0_p = _F_and_∂F_l0p(2p, p, sin_i, cos_i)
        F_2p1_0_p, ~ = _F_and_∂F_l0p(2p + 1, p, sin_i, cos_i)

        den  += prod * (∂F_2p_0_p * cos_i - p * (2p + 1) * F_2p_0_p * sin_i) * C_2p_0
        prod *= fact # .................................................. (R_e / a)^(2p + 1)
        num  += prod * p * F_2p1_0_p * C_2p1_0 * sin_i
        prod *= fact # ............................................. (R_e / a)^(2 * (p + 1))
    end

    # Compute the eccentricity and convert it back to `Float64`.
    e = Float64(2 * num / den)

    # If `e` is positive, the argument of perigee must be at 90°. Otherwise, it must be at
    # 270°. Refer to [3, p. 1297] for more information.
    if e >= 0
        ω = π / 2
    else
        e = -e
        ω = 3π / 2
    end

    return e, ω
end

############################################################################################
#                                    Private Functions                                     #
############################################################################################

"""
    _frozen_orbit__default_gravity_model() -> AbstractGravityModel

Return the default gravity model used by [`frozen_orbit`](@ref) (EGM96). The model is loaded
at the first call and kept in memory. This function is thread-safe.
"""
function _frozen_orbit__default_gravity_model()
    return lock(_FROZEN_ORBIT_DEFAULT_GRAVITY_MODEL_LOCK) do
        if isnothing(_FROZEN_ORBIT_DEFAULT_GRAVITY_MODEL[])
            _FROZEN_ORBIT_DEFAULT_GRAVITY_MODEL[] = GravityModels.load(
                IcgemFile, fetch_icgem_file(:EGM96)
            )
        end

        return _FROZEN_ORBIT_DEFAULT_GRAVITY_MODEL[]
    end
end

"""
    _frozen_orbit__zonal_unnormalization_factor(::Val{:full}, l::Integer) -> Float64
    _frozen_orbit__zonal_unnormalization_factor(::Val{:schmidt}, l::Integer) -> Float64
    _frozen_orbit__zonal_unnormalization_factor(::Val{:unnormalized}, l::Integer) -> Float64

Return the factor that converts the zonal coefficient (order 0) of degree `l` to the
unnormalized coefficient given the normalization used in the gravity model, which is
returned by the function `GravityModels.coefficient_norm`.

The factor of the full normalization is `√(2l + 1)` for the zonal terms, whereas the Schmidt
quasi-normalization does not modify them.

# Extended help

## Throws

- `ArgumentError`: If the normalization is not supported.
"""
_frozen_orbit__zonal_unnormalization_factor(::Val{:full}, l::Integer) = √(2l + 1)
_frozen_orbit__zonal_unnormalization_factor(::Val{:schmidt}, l::Integer) = 1.0
_frozen_orbit__zonal_unnormalization_factor(::Val{:unnormalized}, l::Integer) = 1.0

function _frozen_orbit__zonal_unnormalization_factor(coefficient_norm, l::Integer)
    return throw(
        ArgumentError(
            "The coefficient normalization `$coefficient_norm` of the gravity model is not supported.",
        ),
    )
end

"""
    _F_and_∂F_l0p(l::Integer, p::Integer, sin_i::BigFloat, cos_i::BigFloat) -> BigFloat, BigFloat

Compute the inclination function `F_{l,0,p}(i)` and its derivative `∂F_{l,0,p} / ∂i` as
defined in **[1, p. 642]**. The inclination is specified by its sine `sin_i` [-] and cosine
`cos_i` [-], allowing the caller to compute them only once for all the degrees.

!!! note

    Internally, we must use `BigFloat` to allow calculations on high degrees.

# References

- **[1]** Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications. 4th ed.
    Microcosm Press, Hawthorne, CA.
"""
function _F_and_∂F_l0p(l::Integer, p::Integer, sin_i::BigFloat, cos_i::BigFloat)
    # We must using `BigInt` and `BigFloat` to compute for high degree.
    lb = big(l)
    pb = big(p)

    F  = zero(BigFloat)
    ∂F = zero(BigFloat)
    kb = div(lb, 2)

    sin²_i = sin_i * sin_i

    # `k_t` is the summation term of `F_l0p(i)` for a specific `t`. We will compute those
    # terms iteratively. The equation for `k_t` is:
    #
    #                              (2l - 2t)!
    #  k_t =  ───────────────────────────────────────────────────── ⋅ sin(i)^(l - 2t) ⋅ (-1)^{p - t - k} .
    #         t! ⋅ (l - t)! ⋅ (p - t)! ⋅ (l - p - t)! ⋅ 2^(2l - 2t)
    #
    # Thus, if we divide `k_t` / `k_{t - 1}` we get:
    #
    #    k_t         2 ⋅ (p - t + 1) ⋅ (l - p - t + 1)
    #  ───────── = - ───────────────────────────────── .
    #  k_{t - 1}       t ⋅ (2l - 2t + 1) ⋅ sin(i)^2
    #
    #  Finally, we can compute `k_t` iteratively using:
    #
    #          2 ⋅ (p - t + 1) ⋅ (l - p - t + 1)
    #  k_t = - ───────────────────────────────── ⋅ k_{t - 1} ,
    #             t ⋅ (2l - 2t + 1) ⋅ sin(i)^2
    #
    # where the initialization is:
    #
    #                    2l!
    #  k_0 =  ─────────────────────────── ⋅ sin(i)^(l) ⋅ (-1)^{p - k} .
    #         l! ⋅ p! ⋅ (l - p)! ⋅ 2^(2l)
    #
    #        ┌       ┐   ┌       ┐
    #        │  2l   │   │   l   │   sin(i)^(l)
    #  k_0 = │       │ ⋅ │       │ ⋅ ────────── ⋅ (-1)^{p - k} .
    #        │   l   │   │   p   │     2^(2l)
    #        └       ┘   └       ┘
    #
    # !!! note
    #   Thanks to @stevengj and @danielwe for the suggestions:
    #
    #       https://discourse.julialang.org/t/and-julia-keeps-amazing-me-high-precision-computation/107740/6
    #       https://discourse.julialang.org/t/and-julia-keeps-amazing-me-high-precision-computation/107740/10

    fact = abs(pb - kb) % 2 == 0 ? big(1) : big(-1)
    k_t  = binomial(2lb, lb) * binomial(lb, pb) * sin_i^l * fact / 2^(2lb)

    F  = k_t
    ∂F = k_t * lb / sin_i

    for tb in big(1):min(kb, pb)
        k_t *= -2 * (pb - tb + 1) * (lb - pb - tb + 1) / (tb * (2lb - 2tb + 1) * sin²_i)
        F   += k_t
        ∂F  += k_t * (lb - 2tb) / sin_i
    end

    return F, ∂F * cos_i
end
