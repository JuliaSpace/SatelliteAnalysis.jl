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
    will automatically fetch and load the EGM96 gravity model. However, loading a gravity
    model can significantly decrease the performance. Thus, it is advisable to pass a
    gravity model here.
    (**Default** = `nothing`)
- `max_degree`: Maximum gravity model degree used to compute the frozen eccentricity. If it
    is equal to or lower than 0, the maximum degree in `grav_model` will be used. Otherwise,
    if it is lower than 3 or higher than the `grav_model` maximum degree, it will be clamped
    accordingly.
    (**Default** = 53)

# References

- **[1]** Rosborough, G. W.; Ocampo, C. A (1991). Influence of higher degree zonals on the
    frozen orbit geometry. Proceedings of the AAS/AIAA Astrodynamics Conference, Durango,
    CO.

# Extended Help

Due to the Earth's gravitational perturbation, the orbit of a salite will experience
secular changes in the argument of perigee. Hence, the satellite mean altitude per latitude
will differ during the mission. This effect can be problematic, especially if we must
compare images by a camera onboard the satellite in different periods. The altitude
variation will change the resolution, leading to some problems when comparing the data.

We can avoid this problem if we compute an eccentricity and the argument of perigee that
yields theoretically:

    ∂e      ∂ω
    ── = 0, ── = 0
    ∂t      ∂t

This orbit is called **frozen**. Refer to **[1]** for more information.

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
    max_degree::Int = 53
)

    # Fetch EGM-96 gravity model if the user does not specify one.
    if isnothing(gravity_model)
        gm = GravityModels.load(IcgemFile, fetch_icgem_file(:EGM96))
    else
        gm = gravity_model
    end

    # If `max_degree` is negative, use maximum degree available.
    if max_degree <= 0
        max_degree = GravityModels.maximum_degree(grav_model)
    end

    # The maximum degree must be at least 3.
    max_degree = max(3, max_degree)

    # Obtain the maximum `p` given the maximum degree.
    p_max = floor(Int, (max_degree - 1) / 2)

    # If the coefficients are fully normalized, we must unormalize them because the theory
    # in [2].
    unnormalize = GravityModels.coefficient_norm(gm) == :full

    # We must using `BigInt` and `BigFloat` to compute for high degree.
    ab  = big(a)
    ib  = big(i)
    num = zero(BigFloat)
    den = zero(BigFloat)
    R_e = big(GravityModels.radius(gm))

    sin_i, cos_i = sincos(ib)

    # Auxiliar variables to compute (R_e / a)^2p.
    fact = R_e / ab
    prod = fact * fact # ............................................... (R_e / a)^2 (p = 0)

    # Compute the numerator and denominator of the frozen eccentricity as indicated in Eq.
    # 38 [3, p. 1297].
    for p in 1:p_max
        C_2p_0  = -(GravityModels.coefficients(gm, 2p,     0) |> first) |> big
        C_2p1_0 = -(GravityModels.coefficients(gm, 2p + 1, 0) |> first) |> big

        if unnormalize
            C_2p_0  *= √(4p + 1)
            C_2p1_0 *= √(4p + 3)
        end

        F_2p_0_p, ∂F_2p_0_p = _F_and_∂F_l0p(2p, p, ib)
        F_2p1_0_p, ~ = _F_and_∂F_l0p(2p + 1, p, ib)

        den  += prod * (∂F_2p_0_p * cos_i - p * (2p + 1) * F_2p_0_p * sin_i) * C_2p_0
        prod *= fact # .................................................. (R_e / a)^(2p + 1)
        num  += prod * p * F_2p1_0_p * C_2p1_0 * sin_i
        prod *= fact # ............................................. (R_e / a)^(2 * (p + 1))
    end

    # Compute the eccentricy and convert it back to `Float64`.
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
    _F_and_∂F_l0p(l::Integer, p::Integer, i::Number) -> BigFloat, BigFloat

Compute the inclination function `F_{l,0,p}(i)` and its derivative `∂F_{l,0,p} / ∂i`as
defined in **[1, p. 642]**.

!!! note
    Internally, we must use `BigFloat` to allow calculations on high degrees.

# References

- **[1]** Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications. 4th ed.
    Microcosm Press, Hawthorne, CA.
"""
function _F_and_∂F_l0p(l::Integer, p::Integer, i::Number)
    # We must using `BigInt` and `BigFloat` to compute for high degree.
    lb = big(l)
    pb = big(p)
    ib = big(i)

    F  = zero(BigFloat)
    ∂F = zero(BigFloat)
    kb = div(lb, 2)

    sin_i, cos_i = sincos(ib)

    for tb in big(0):min(kb, pb)
        k_t = factorial(2lb - 2tb) / (
            factorial(tb) *
            factorial(lb - tb) *
            factorial(pb - tb) *
            factorial(lb - pb - tb) *
            big(2)^(2lb - 2tb)
        )

        # This fact simulates the term (-1)^(p - t - k).
        fact = abs(pb - tb - kb) % 2 == 0 ? big(1) : big(-1)

        F  += k_t * sin_i^(lb - 2tb) * fact
        ∂F += k_t * (lb - 2tb) * sin_i^(lb - 2tb - 1) * fact
    end

    return F, ∂F * cos_i
end
