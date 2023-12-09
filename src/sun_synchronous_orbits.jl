## Description #############################################################################
#
# Functions to design Sun synchronous orbits.
#
## References ##############################################################################
#
# [1] Kozai, Y (1959). The Motion of a Close Earth Satellite. The Astronomical Journal,
#     v. 64, no. 1274, pp. 367 -- 377.
#
############################################################################################

export design_sun_sync_ground_repeating_orbit
export sun_sync_orbit_from_angular_velocity
export sun_sync_orbit_semi_major_axis
export sun_sync_orbit_inclination

"""
    design_sun_sync_ground_repeating_orbit(minimum_repetition::Int, maximum_repetition::Int; kwargs...) -> DataFrame

List all the Sun synchronous, ground repeating orbits in which their repetition period is in
the interval `[minimum_repetition, maximum_repetition]` days.

This function returns a `DataFrame` with the following columns:

- `semi_major_axis`: Orbit semi-major axis.
- `altitude`: Orbit altitude above the Equator `(a - R0)`.
- `inclination`: Orbit inclination.
- `period`: Orbital period.
- `rev_per_days`: If the keyword `pretify_rev_per_days` is `false`, this column contains
    `Tuple`s with the integer and rational parts of the number of revolutions per day.
    Otherwise, it contains a string with a prety representation of the number of revolutions
    per day.
- `adjacent_gt_distance`: Distance between two adjacent ground tracks at Equator.
- `adjacent_gt_angle`: Angle between two adjacent ground tracks at Equator measured from the
    satellite position.

!!! note
    The units of those columns depends on the keywords.

# Keywords

- `angle_unit::Symbol`: Unit for all the angles in the output `DataFrame`.  It can be `:deg`
    for degrees or `:rad` for radians. (**Default**: `:deg`)
- `distance_unit::Symbol`: The unit for all the distances in the output `DataFrame`. It can
    be `:m` for meters or `:km` for kilometers.  (**Default**: `:km`)
- `e::Number`: Orbit eccentricity. (**Default**: 0)
- `int_rev_per_day::Tuple`: `Tuple` with the integer parts of the number of revolutions per
    day to be analyzed. (**Default** = `(13, 14, 15, 16, 17)`)
- `pretity_rev_per_days::Bool`: If `true`, the column with the revolutions per day will be
    conveted to a string with a pretty representation of this information. (**Default**:
    `true`)
- `maximum_altitude::Union{Nothing, Number}`: Maximum altitude [m] of the orbits in the
    output `DataFrame`. If it is `nothing`, the algorithm will not apply a higher limit to
    the orbital altitude. (**Default** = `nothing`)
- `minimum_altitude::Union{Nothing, Number}`: Minimum altitude [m] of the orbits in the
    output `DataFrame`. If it is `nothing`, the algorithm will not apply a lower limit to
    the orbital altitude. (**Default** = `nothing`)
- `time_unit::Symbol`: Unit for all the time values in the output `DataFrame`.  It can be
    `:s` for seconds, `:m` for minutes, or `:h` for hours.  (**Default** = `:h`)
- `m0::Number`: Standard gravitational parameter for Earth [m³ / s²].
    (**Default** = `GM_EARTH`)
- `J2::Number`: J₂ perturbation term. (**Default** = EGM_2008_J2)
- `R0::Number`: Earth's equatorial radius [m]. (**Default** = EARTH_EQUATORIAL_RADIUS)
"""
function design_sun_sync_ground_repeating_orbit(
    minimum_repetition::Int,
    maximum_repetition::Int;
    angle_unit::Symbol = :deg,
    distance_unit::Symbol = :km,
    e::Number = 0,
    int_rev_per_day::Tuple = (13, 14, 15, 16, 17),
    pretify_rev_per_days::Bool = true,
    maximum_altitude::Union{Nothing, Number} = nothing,
    minimum_altitude::Union{Nothing, Number} = nothing,
    time_unit::Symbol = :m,
    # Constants.
    J2::Number = EGM_2008_J2,
    m0::Number = GM_EARTH,
    R0::Number = EARTH_EQUATORIAL_RADIUS
)
    R₀ = EARTH_EQUATORIAL_RADIUS

    # Check if the inputs are valid.
    minimum_repetition <= 0 && throw(
        ArgumentError("The minimum repetition must be greater than 0.")
    )

    maximum_repetition <= 0 && throw(
        ArgumentError("The maximum repetition must be greater than 0.")
    )

    maximum_repetition < minimum_repetition && throw(
        ArgumentError("The minimum repetition must be smaller or equal than the maximum repetition.")
    )

    !(0 <= e < 1) && throw(
        ArgumentError("The eccentricity must be within the interval [0, 1).")
    )

    # Create an empty `DataFrame` that will store the list of orbits.
    df = DataFrame(
        semi_major_axis      = Float64[],
        altitude             = Float64[],
        inclination          = Float64[],
        period               = Float64[],
        rev_per_days         = pretify_rev_per_days ? String[] : Tuple{Int, Rational}[],
        adjacent_gt_distance = Float64[],
        adjacent_gt_angle    = Float64[]
    )

    # Check the units for the values.
    dunit   = distance_unit === :km  ? 1.e-3     : 1.0
    angunit = angle_unit    === :deg ? 180.0 / π : 1.0
    tunit   = begin
        if time_unit == :h
            1.0 / 3600
        elseif time_unit == :m
            1.0 / 60
        else
            1.0
        end
    end

    # Loop through the possible repetition times.
    for den in minimum_repetition:maximum_repetition
        for num in 0:(den - 1)
            # Check if the fraction `num / den` is irreducible.
            gcd(num, den) != 1 && continue

            # Loop through the integer parts.
            for int in int_rev_per_day
                # Compute the number of revolutions per day of this orbit, and convert it to
                # angular velocity.
                num_rev_per_day = int + num / den
                n = num_rev_per_day * 2π / 86400

                # Find a Sun synchronous orbit with that angular velocity.
                a, i, converged = sun_sync_orbit_from_angular_velocity(
                    n,
                    e;
                    no_warnings = true,
                    J2 = J2,
                    m0 = m0,
                    R0 = R0
                )

                # If the algorithm has not converged or if the orbit is not valid, skip this
                # value.
                orbit_valid = a * (1 - e) > R₀
                (!converged || !orbit_valid) && continue

                # If we reach this point, add the orbit to the `DataFrame`.
                orb_angvel = orbital_angular_velocity(
                    a,
                    e,
                    i;
                    perturbation = :J2,
                    J2 = J2,
                    m0 = m0,
                    R0 = R0
                )

                orb_period = 2π / orb_angvel
                orb_cycle = num == 0 ? 1 : den

                h = a - R₀

                push!(df, (
                    a * dunit,
                    h * dunit,
                    i * angunit,
                    orb_period * tunit,
                    pretify_rev_per_days ?
                        _pretify_rev_per_days(int, num, den) :
                        (int, num // den),
                    ground_repeating_orbit_adjacent_track_distance(
                        orb_period,
                        i,
                        orb_cycle
                    ) * dunit,
                    ground_repeating_orbit_adjacent_track_angle(
                        h,
                        orb_period,
                        i,
                        orb_cycle
                    ) * angunit
                ))
            end
        end
    end

    # Filter the `DataFrame` with respect to the altitude.
    filter!(
        r -> begin
            if !isnothing(minimum_altitude) && (r[:altitude] < minimum_altitude * dunit)
                return false
            end

            if !isnothing(maximum_altitude) && (r[:altitude] > maximum_altitude * dunit)
                return false
            end

            return true
        end,
        df
    )

    # Sort the `DataFrame` by the semi-major axis.
    sort!(df, :semi_major_axis)

    return df
end

"""
    sun_sync_orbit_from_angular_velocity(angvel::T1, e::T2 = 0; kwargs...) where {T1 <: Number, T2 <: Number} -> T, T, Bool

Compute the Sun-synchronous orbit semi-major axis [m] and inclination [rad] given the
angular velocity `angvel` [rad / s] and the orbit eccentricity `e` [ ]. If the latter is
omitted, the orbit is considered circular, i.e., `e = 0`.

The algorithm here considers only the perturbation terms up to J₂.

!!! note
    Internally, this function uses the precision obtained by promoting `T1` and `T2` to a
    float-pointing number `T`.

# Keywords

- `max_iterations::Number`: Maximum number of iterations in the Newton-Raphson method.
    (**Default** = 3)
- `no_warnings::Bool`: If `true`, no warnings will be printed. (**Default** = `false`)
- `tolerance::Union{Nothing, NTuple{2, Number}}`: Residue tolerances to verify if the
    numerical method has converged. If it is `nothing`, `(√eps(T), √eps(T))` will be used,
    where `T` is the internal type for the computations. Notice that the residue function
    `f₁` unit is [deg / day], whereas the `f₂` unit is [deg / min]. (**Default** = 1e-18)
- `m0::Number`: Standard gravitational parameter for Earth [m³ / s²].
    (**Default** = `GM_EARTH`)
- `J2::Number`: J₂ perturbation term. (**Default** = EGM_2008_J2)
- `R0::Number`: Earth's equatorial radius [m]. (**Default** = EARTH_EQUATORIAL_RADIUS)

# Returns

- `T`: Semi-major axis [m].
- `T`: Inclination [rad].
- `Bool`: `true` if the Newton-Raphson algorithm converged, or `false` otherwise.

# Extended help

A Sun-synchronous orbit is defined as an orbit in which the precession of the right
ascension of the ascending node (RAAN) equals the Earth's orbit mean motion. In this case,
the orbit plane will have the same orientation to the Sun at the ascending node.

The RAAN time-derivative considering only the secular terms up to J₂ is [1, p. 372] is:

```
∂Ω      3                       n̄
── = - ─── R₀² . J₂ . cos(i) . ─── .
∂t      2                       p²
```

where:

```
         ┌                                                ┐
         │      3    R₀²                                  │
n̄ = n₀ . │ 1 + ─── . ─── . J₂ . √(1 - e²) . (2 - 3sin²(i))│.
         │      4     p²                                  │
         └                                                ┘
```

We can express the orbit angular velocity in terms of its nodal period, i.e., the period
it takes for the satellite to cross the ascending node two consecutive times:

```
          ∂M     ∂ω
angvel = ──── + ────,
          ∂t     ∂t

                  3    R₀²
angvel = n̄ + n̄ . ─── . ─── . J₂ . (4 - 5sin²(i)),
                  4     p²
```

where `n` is the perturbed mean motion due to the same consideration as presented for the
RAAN time-derivative.

Finally, this function finds the pair `(a, i)` that simultaneously solves the equations:

```
∂Ω
── (a, i) = EARTH_ORBIT_MEAN_MOTION,
∂t

 ∂M            ∂ω
──── (a, i) + ──── (a, i) = angvel,
 ∂t            ∂t
```

using the Newton-Raphson method with the presented equations.

## Examples

```julia-repl
julia> sun_sync_orbit_from_angular_velocity(0.06 |> deg2rad)
(7.130983932846816e6, 1.7175898375139984, true)

julia> sun_sync_orbit_from_angular_velocity(0.06 |> deg2rad, 0)
(7.130983932846816e6, 1.7175898375139984, true)

julia> sun_sync_orbit_from_angular_velocity(0.06 |> deg2rad, 0.1)
(7.13086251587883e6, 1.7146410689929386, true)
```

The user can verify some internal information of the solver by turning on the debugging
logs:

```julia-repl
julia> with_logger(ConsoleLogger(stderr, Logging.Debug)) do
           sun_sync_orbit_from_angular_velocity(0.06 |> deg2rad)
       end
┌ Debug: Iteration #1
│   Estimation :
│     a  = 7136.635455699327 km
│     i  = 81.57099271530629 °
│   Residues :
│     f₁ = 1.9706966881670205 ° / day
│     f₂ = 0.004266929859281898 ° / min
└ @ SatelliteAnalysis ~/.julia/dev/SatelliteAnalysis/src/sun_synchronous_orbits.jl:394
┌ Debug: Iteration #2
│   Estimation :
│     a  = 7128.856266265137 km
│     i  = 98.46515928332974 °
│   Residues :
│     f₁ = -0.0073785260175135425 ° / day
│     f₂ = -0.0016144146737784304 ° / min
└ @ SatelliteAnalysis ~/.julia/dev/SatelliteAnalysis/src/sun_synchronous_orbits.jl:394
┌ Debug: Iteration #3
│   Estimation :
│     a  = 7130.983594940013 km
│     i  = 98.41070350863473 °
│   Residues :
│     f₁ = -6.549620124363109e-6 ° / day
│     f₂ = -2.6066638092459016e-7 ° / min
└ @ SatelliteAnalysis ~/.julia/dev/SatelliteAnalysis/src/sun_synchronous_orbits.jl:394
┌ Debug: Iteration #4
│   Estimation :
│     a  = 7130.983932846698 km
│     i  = 98.41064862339992 °
│   Residues :
│     f₁ = 8.290634845309341e-11 ° / day
│     f₂ = -2.2648549702353193e-14 ° / min
└ @ SatelliteAnalysis ~/.julia/dev/SatelliteAnalysis/src/sun_synchronous_orbits.jl:394
(7.130983932846816e6, 1.7175898375139984, true)
```

## References

- **[1]** Kozai, Y (1959). The Motion of a Close Earth Satellite. The Astronomical Journal,
    v. 64, no. 1274, pp. 367 -- 377.
"""
function sun_sync_orbit_from_angular_velocity(
    angvel::T1,
    e::T2 = 0;
    max_iterations::Number = 30,
    no_warnings::Bool = false,
    tolerance::Union{Nothing, NTuple{2, Number}} = nothing,
    # Constants.
    J2::Number = EGM_2008_J2,
    m0::Number = GM_EARTH,
    R0::Number = EARTH_EQUATORIAL_RADIUS
) where {T1 <: Number, T2 <: Number}

    T = float(promote_type(T1, T2))

    # Constant to convert [rad / s] to [deg / day].
    rs_to_dd = T(86400 * 180 / π)

    # Constant to convert [rad / s] to [deg / min].
    rs_to_dm = T(60 * 180 / π)

    # Check if the arguments are valid.
    angvel <= 0 && throw(ArgumentError("The angular velocity must be greater than 0."))
    !(0 <= e < 1) && throw(ArgumentError("The eccentricity must be within the interval [0, 1)."))

    # Obtain the tolerance.
    tol = isnothing(tolerance) ? (√eps(T), √eps(T)) : (T(tolerance[1]), T(tolerance[2]))

    # Auxiliary variables.
    μ  = T(m0)
    R₀ = T(R0)
    J₂ = T(J2)
    β² = 1 - T(e)^2
    β  = √β²
    β³ = β² * β
    β⁴ = β² * β²

    # Auxiliary constant to compute the functions.
    k₅ = √(μ / R₀^3)
    k₁ = -(3 // 2) * J₂ * k₅ / β⁴
    k₂ = +(3 // 4) * J₂ / β³
    k₄ = -k₁ / 2
    k₃ = k₄ * β
    k₆ = ((3 // 4) * J₂ / β⁴)^2 * k₅ * β

    # We will change the units of the desired values to improve the numerical stability.
    k₁ *= rs_to_dd
    k₃ *= rs_to_dm
    k₄ *= rs_to_dm
    k₅ *= rs_to_dm
    k₆ *= rs_to_dm

    # The desired RAAN time-derivative is equal to the Earth's orbit mean motion
    # [deg / day].
    Ω̇_d = T(EARTH_ORBIT_MEAN_MOTION) * rs_to_dd

    # The desired angular velocity in [deg / min].
    ω_d = T(angvel) * rs_to_dm

    # == Solving for Zeros of f₁ and f₂ Using Newton-Raphson Method ========================
    #
    # The function `f₁` is the residue related to the time derivative of RAAN, and the
    # function `f₂` is the residue related to the angular velocity.

    # Initial guess based on the unperturbed model. Notice that we will estimate
    # `1 / √(a / R₀)` and `cos(i)`
    isqrt_ā = (ω_d / k₅)^(1 // 3)
    cos_i   = -Ω̇_d / isqrt_ā^7 / k₁

    # By setting the initial values of `f1` and `f2` to `10tol`, we assure that the loop
    # will be executed at least one time.
    f₁ = 10T(tol[1])
    f₂ = 10T(tol[2])

    # Loop
    it = 1
    converged = true

    while (abs(f₁) > tol[1]) || (abs(f₂) > tol[2])
        isqrt_ā² = isqrt_ā  * isqrt_ā
        isqrt_ā³ = isqrt_ā² * isqrt_ā
        isqrt_ā⁴ = isqrt_ā² * isqrt_ā²
        isqrt_ā⁶ = isqrt_ā³ * isqrt_ā³
        isqrt_ā⁷ = isqrt_ā⁴ * isqrt_ā³
        isqrt_ā⁸ = isqrt_ā⁴ * isqrt_ā⁴

        cos²_i = cos_i  * cos_i
        cos³_i = cos²_i * cos_i
        cos⁴_i = cos²_i * cos²_i

        c₁ = 3cos²_i - 1
        c₂ = 5cos²_i - 1
        c₃ = c₁ * c₂
        c₄ = k₃ * c₁ + k₄ * c₂

        # Compute the residue using the current estimates.
        f₁ = Ω̇_d - k₁ * cos_i * isqrt_ā⁷ * (1 + k₂ * c₁ * isqrt_ā⁴)
        f₂ = ω_d - (k₅ + c₄ * isqrt_ā⁴ + k₆ * c₃ * isqrt_ā⁸) * isqrt_ā³

        @debug """
        Iteration #$it
          Estimation :
            a  = $(R₀ / (isqrt_ā * isqrt_ā) / 1000) km
            i  = $(abs(cos_i) <= 1 ? acosd(cos_i) : "INVALID") °
          Residues :
            f₁ = $(f₁) ° / day
            f₂ = $(f₂) ° / min
        """

        # Compute the Jacobian.
        ∂f₁_∂isqrt_a = -k₁ * (7cos_i + 11k₂ * (3cos²_i - cos_i) * isqrt_ā⁴) * isqrt_ā⁶
        ∂f₁_∂cos_i   = -k₁ * (1 + k₂ * (9cos²_i - 1) * isqrt_ā⁴) * isqrt_ā⁷
        ∂f₂_∂isqrt_a = -(3k₅ + 7c₄ * isqrt_ā⁴ + 11k₆ * c₃ * isqrt_ā⁸) * isqrt_ā²
        ∂f₂_∂cos_i   = -((6k₃ + 10k₄) * cos_i + k₆ * (-16cos_i + 60cos³_i) * isqrt_ā⁴) * isqrt_ā⁷

        J = @SMatrix T[
            ∂f₁_∂isqrt_a ∂f₁_∂cos_i
            ∂f₂_∂isqrt_a ∂f₂_∂cos_i
        ]

        # Compute the new estimate using the Newton-Raphson method.
        isqrt_ā, cos_i = @SVector([isqrt_ā, cos_i]) - inv(J) * @SVector([f₁, f₂])

        # If the maximum number of iterations allowed has been reached, then
        # indicate that the solution did not converged and exit loop.
        if (it >= max_iterations)
            converged = false
            break
        end

        it += 1
    end

    # If `cos_i` absolute value is larger than 1, the solution does not have physical
    # meaning.
    (abs(cos_i) > 1) &&
        throw(ArgumentError(
            "It is not possible to find a Sun-synchronous orbit with the selected parameters (ang. vel = $(angvel) rad / s, e = $e)."
        ))

    a = R₀ / (isqrt_ā * isqrt_ā)
    i = acos(cos_i)

    # Check if the orbit is valid.
    ((a * (1 - e) < R₀) & !no_warnings) &&
        @warn("The orbit is not valid because the perigee is inside the Earth.")

    !converged && @warn("""
        The algorithm to compute the Sun-synchronous orbit has not converged!
        Residues :
          f₁ = $f₁ ° / day
          f₂ = $f₂ ° / min"""
    )

    # Return.
    return a, i, converged
end

"""
    sun_sync_orbit_semi_major_axis(i::T1, e::T2 = 0; kwargs...) where {T1 <: Number, T2 <: Number} -> T, Bool

Compute the semi-major axis [m] of the Sun-synchronous orbit with inclination `i` [rad] and
the eccentricity `e` [ ]. If the latter is omitted, the orbit is considered circular, i.e.,
`e = 0`.

The algorithm here considers only the perturbation terms up to J₂.

!!! note
    Internally, this function uses the precision obtained by promoting `T1` and `T2` to a
    float-pointing number `T`.

# Keywords

- `max_iterations::Number`: Maximum number of iterations in the Newton-Raphson method.
    (**Default** = 30)
- `tolerance::Union{Nothing, Number}`: Residue tolerance to verify if the numerical method
    has converged. If it is `nothing`, `√eps(T)` will be used, where `T` is the internal
    type for the computations. Notice that the residue unit is [deg / day].
    (**Default** = nothing)
- `m0::Number`: Standard gravitational parameter for Earth [m³ / s²].
    (**Default** = `GM_EARTH`)
- `J2::Number`: J₂ perturbation term. (**Default** = EGM_2008_J2)
- `R0::Number`: Earth's equatorial radius [m]. (**Default** = EARTH_EQUATORIAL_RADIUS)

# Returns

- `T`: Semi-major axis [m] of the Sun-synchronous orbit with inclination `i` and
    eccentricity `e`.
- `Bool`: `true` if the Newton-Raphson algorithm converged, or `false` otherwise.

# Extended help

A Sun-synchronous orbit is defined as an orbit in which the precession of the right
ascension of the ascending node (RAAN) equals the Earth's orbit mean motion. In this case,
the orbit plane will have the same orientation to the Sun at the ascending node.

The RAAN time-derivative considering only the secular terms up to J₂ is [1, p. 372] is:

```
∂Ω      3                       n̄
── = - ─── R₀² . J₂ . cos(i) . ─── .
∂t      2                       p²
```

where:

```
         ┌                                                ┐
         │      3    R₀²                                  │
n̄ = n₀ . │ 1 + ─── . ─── . J₂ . √(1 - e²) . (2 - 3sin²(i))│.
         │      4     p²                                  │
         └                                                ┘
```

Finally, this function solves the equation:

```
∂Ω
── (a) = EARTH_ORBIT_MEAN_MOTION
∂t
```

for `a` using the Newton-Raphson method with the presented equations.

## Examples

```julia-repl
julia> sun_sync_orbit_semi_major_axis(98.410 |> deg2rad)
(7.130827866508738e6, true)

julia> sun_sync_orbit_semi_major_axis(98.410 |> deg2rad, 0)
(7.130827866508738e6, true)

julia> sun_sync_orbit_semi_major_axis(98.410 |> deg2rad, 0.001111)
(7.1308328955274355e6, true)
```

The user can verify some internal information of the solver by turning on the debugging
logs:

```julia-repl
julia> with_logger(ConsoleLogger(stderr, Logging.Debug)) do
           sun_sync_orbit_semi_major_axis(98.41064163374567 |> deg2rad)
       end
┌ Debug: Iteration #1
│   Estimation : 7130.981820550704 km
│   Residue    : 0.0005989504045072862 ° / day
└ @ SatelliteAnalysis ~/.julia/dev/SatelliteAnalysis/src/sun_synchronous_orbits.jl:559
┌ Debug: Iteration #2
│   Estimation : 7130.982250931794 km
│   Residue    : -2.081337784770338e-7 ° / day
└ @ SatelliteAnalysis ~/.julia/dev/SatelliteAnalysis/src/sun_synchronous_orbits.jl:559
┌ Debug: Iteration #3
│   Estimation : 7130.982250931845 km
│   Residue    : -2.4312691616901194e-14 ° / day
└ @ SatelliteAnalysis ~/.julia/dev/SatelliteAnalysis/src/sun_synchronous_orbits.jl:559
(7.130982250931845e6, true)
```

## References

- **[1]** Kozai, Y (1959). The Motion of a Close Earth Satellite. The Astronomical Journal,
    v. 64, no. 1274, pp. 367 -- 377.
"""
function sun_sync_orbit_semi_major_axis(
    i::T1,
    e::T2 = 0;
    max_iterations::Number = 30,
    tolerance::Union{Nothing, Number} = nothing,
    # Constants.
    J2::Number = EGM_2008_J2,
    m0::Number = GM_EARTH,
    R0::Number = EARTH_EQUATORIAL_RADIUS
) where {T1 <: Number, T2 <: Number}

    T = float(promote_type(T1, T2))

    # Constant to convert [rad / s] to [deg / day].
    rs_to_dd = T(86400 * 180 / π)

    # Check the inputs.
    !(0 <= e < 1) && throw(ArgumentError("The eccentricity must be within the interval [0, 1)."))

    # Obtain the tolerance.
    tol = isnothing(tolerance) ? √(eps(T)) : T(tolerance)

    # The RAAN time-derivative considering only the terms up to J₂ is [1, p. 372]:
    #
    #   ∂Ω      3                       n̄
    #   ── = - ─── R₀² . J₂ . cos(i) . ─── .
    #   ∂t      2                       p²
    #
    # where:
    #
    #             ┌                                                ┐
    #             │      3    R₀²                                  │
    #    n̄ = n₀ . │ 1 + ─── . ─── . J₂ . √(1 - e²) . (2 - 3sin²(i))│.
    #             │      4     p²                                  │
    #             └                                                ┘

    # Auxiliary variables.
    R₀  = T(R0)
    kR₀ = √(R₀^3)
    μ   = T(m0)
    J₂  = T(J2)
    β²  = 1 - T(e)^2
    β   = √β²
    β³  = β² * β
    β⁴  = β² * β²

    sin_i, cos_i = sincos(T(i))
    k₁ = -(3 // 2) * J₂ * √μ * cos_i / (kR₀ * β⁴)
    k₂ = +(3 // 4) * J₂ * (2 - 3sin_i^2) / β³

    # If `k₁` is negative, Ω̇ will be also be negative. Hence, it is impossible to find a
    # Sun-synchronous orbit.
    k₁ < 0 &&
        throw(ArgumentError(
            "It is not possible to find a Sun-synchronous orbit with the selected parameters (i = $(rad2deg(i))°, e = $e)."
        ))

    # We will change the units of the desired values to improve the numerical stability.
    k₁ *= rs_to_dd

    # The desired RAAN time-derivative is equal to the Earth's orbit mean motion
    # [deg / day].
    Ω̇_d = T(EARTH_ORBIT_MEAN_MOTION) * rs_to_dd

    # The initial guess is computed by neglecting the mean motion secular perturbation in
    # the RAAN time-derivative. Notice that we will estimate `1 / √(a / R₀)`.
    isqrt_ā = (Ω̇_d / k₁)^(1 // 7)

    # By setting the initial values of `f₁` to `10tol`, we assure that the loop will be
    # executed at least one time.
    f₁ = 10T(tol)

    # Loop.
    it = 1
    converged = true

    while abs(f₁) > tol
        isqrt_ā³ = isqrt_ā  * isqrt_ā  * isqrt_ā
        isqrt_ā⁴ = isqrt_ā³ * isqrt_ā
        isqrt_ā⁶ = isqrt_ā³ * isqrt_ā³
        isqrt_ā⁷ = isqrt_ā⁶ * isqrt_ā

        # Compute the residue using the current estimate.
        f₁ = Ω̇_d - k₁ * (1 + k₂ * isqrt_ā⁴) * isqrt_ā⁷

        @debug """
        Iteration #$it
          Estimation : $(1 / (isqrt_ā * isqrt_ā) * R₀ / 1000) km
          Residue    : $(f₁) ° / day
        """

        # Compute the derivative.
        ∂f₁_∂isqrt_a = -k₁ * (7 + 11k₂ * isqrt_ā⁴) * isqrt_ā⁶

        # Update the estimate using the Newton-Raphson algorithm.
        isqrt_ā = isqrt_ā - f₁ / ∂f₁_∂isqrt_a

        # If the maximum number of iterations allowed has been reached, indicate that the
        # solution did not converged and exit loop.
        if (it >= max_iterations)
            @debug "The algorithm reached the maximum number of iterations without converging."
            converged = false
            break
        end

        it += 1
    end

    !converged && @warn("""
        The algorithm to compute the Sun-synchronous orbit semi-major axis has not converged!
        Residue: $f₁ ° / day """
    )

    # Compute the semi-major axis from the normalized semi-major axis.
    a = R₀ / (isqrt_ā * isqrt_ā)

    # Check if the orbit is valid.
    (a * (1 - e) < R₀) &&
        @warn("The orbit is not valid because the perigee is inside the Earth.")

    return a, converged
end

"""
    sun_sync_orbit_inclination(a::T1, e::T2 = 0; kwargs...) where {T1 <: Number, T2 <: Number} -> T, Bool

Compute the inclination [rad] of the Sun-synchronous orbit with semi-major axis `a` [m] and
the eccentricity `e` [ ]. If the latter is omitted, the orbit is considered circular, i.e.,
`e = 0`.

The algorithm here considers only the perturbation terms up to J₂.

!!! note
    Internally, this function uses the precision obtained by promoting `T1` and `T2` to a
    float-pointing number `T`.

# Keywords

- `max_iterations::Number`: Maximum number of iterations in the Newton-Raphson method.
    (**Default** = 30)
- `tolerance::Union{Nothing, Number}`: Residue tolerance to verify if the numerical method
    has converged. If it is `nothing`, `√eps(T)` will be used, where `T` is the internal
    type for the computations. Notice that the residue unit is [deg / day].
    (**Default** = nothing)
- `m0::Number`: Standard gravitational parameter for Earth [m³ / s²].
    (**Default** = `GM_EARTH`)
- `J2::Number`: J₂ perturbation term. (**Default** = EGM_2008_J2)
- `R0::Number`: Earth's equatorial radius [m]. (**Default** = EARTH_EQUATORIAL_RADIUS)

# Returns

- `T`: Inclination [rad] of the Sun-synchronous orbit with semi-major axis `a` and
    eccentricity `e`.
- `Bool`: `true` if the Newton-Raphson algorithm converged, or `false` otherwise.

# Extended help

A Sun-synchronous orbit is defined as an orbit in which the precession of the right
ascension of the ascending node (RAAN) equals the Earth's orbit mean motion. In this case,
the orbit plane will have the same orientation to the Sun at the ascending node.

The RAAN time-derivative considering only the secular terms up to J₂ is [1, p. 372] is:

```
∂Ω      3                       n̄
── = - ─── R₀² . J₂ . cos(i) . ─── .
∂t      2                       p²
```

where:

```
         ┌                                                ┐
         │      3    R₀²                                  │
n̄ = n₀ . │ 1 + ─── . ─── . J₂ . √(1 - e²) . (2 - 3sin²(i))│.
         │      4     p²                                  │
         └                                                ┘
```

Finally, this function solves the equation:

```
∂Ω
── (i) = EARTH_ORBIT_MEAN_MOTION
∂t
```

for `i` using the Newton-Raphson method with the presented equations.

## Examples

```julia-repl
julia> sun_sync_orbit_inclination(7130.982e3)
(1.7175896973066611, true)

julia> sun_sync_orbit_inclination(7130.982e3, 0)
(1.7175896973066611, true)

julia> sun_sync_orbit_inclination(7130.982e3, 0.001111)
(1.7175893324980402, true)
```

The user can verify some internal information of the solver by turning on the debugging
logs:

```julia
julia> with_logger(ConsoleLogger(stderr, Logging.Debug)) do
           sun_sync_orbit_inclination(7130.982e3)
       end
┌ Debug: Iteration #1
│   Estimation : 98.41064059121584 °
│   Residue    : 0.0005992085524891833 ° / day
└ @ SatelliteAnalysis ~/.julia/dev/SatelliteAnalysis/src/sun_synchronous_orbits.jl:686
┌ Debug: Iteration #2
│   Estimation : 98.41064059082426 °
│   Residue    : -4.556321986370904e-11 ° / day
└ @ SatelliteAnalysis ~/.julia/dev/SatelliteAnalysis/src/sun_synchronous_orbits.jl:686
(1.7175896973066611, true)
```

## References

- **[1]** Kozai, Y (1959). The Motion of a Close Earth Satellite. The Astronomical Journal,
    v. 64, no. 1274, pp. 367 -- 377.
"""
function sun_sync_orbit_inclination(
    a::T1,
    e::T2 = 0;
    max_iterations::Number = 30,
    tolerance::Union{Nothing, Number} = nothing,
    # Constants.
    J2::Number = EGM_2008_J2,
    m0::Number = GM_EARTH,
    R0::Number = EARTH_EQUATORIAL_RADIUS
) where {T1 <: Number, T2 <: Number}

    T  = float(promote_type(T1, T2))
    R₀ = T(EARTH_EQUATORIAL_RADIUS)
    rs_to_dd = T(86400 * 180 / π)

    # Check if the arguments are valid.
    !(0 <= e < 1) && throw(ArgumentError("The eccentricity must be within the interval [0, 1)."))
    (a * (1 - e) <= R₀) && throw(ArgumentError("The perigee must be larger than the Earth's radius."))

    # Obtain the tolerance.
    tol = isnothing(tolerance) ? √(eps(T)) : T(tolerance)

    # The RAAN time-derivative considering only the terms up to J₂ is [1, p. 372]:
    #
    #   ∂Ω      3                       n̄
    #   ── = - ─── R₀² . J₂ . cos(i) . ─── .
    #   ∂t      2                       p²
    #
    # where:
    #
    #             ┌                                                ┐
    #             │      3    R₀²                                  │
    #    n̄ = n₀ . │ 1 + ─── . ─── . J₂ . √(1 - e²) . (2 - 3sin²(i))│.
    #             │      4     p²                                  │
    #             └                                                ┘

    # Auxiliary variables.
    μ  = T(m0)
    J₂ = T(J2)
    β² = 1 - T(e)^2
    β  = √β²
    n₀ = √(μ / a^3)
    p  = a * β²
    k₁ = -(3 // 2) * J₂ * (R₀ / p)^2 * n₀
    k₂ = +(3 // 4) * J₂ * (R₀ / p)^2 * β

    # We will change the units of the desired values to improve the numerical stability.
    k₁ *= rs_to_dd

    # We can reduce the numerical errors by switching the unit from [rad / s] to
    # [deg / day].
    A = (3 * k₁ * k₂)
    B = k₁ * (1 - k₂)

    # The desired RAAN time-derivative is equal to the Earth's orbit mean motion
    # [deg / day].
    Ω̇_d = T(EARTH_ORBIT_MEAN_MOTION) * rs_to_dd

    # The initial guess is computed by neglecting the mean motion secular perturbation in
    # the RAAN time-derivative.  Notice that we will estimate `cos(i)`.
    cos_i = Ω̇_d / k₁

    # By setting the initial values of `f₁` to `10tol`, we assure that the loop will be
    # executed at least one time.
    f₁ = 10T(tol)

    # We will solve for `cos(i)` using the Newton-Raphson method.
    it = 1
    converged = true

    while abs(f₁) > tol
        cos²_i = cos_i  * cos_i
        cos³_i = cos²_i * cos_i

        # Compute the residue using the current estimate.
        f₁ = Ω̇_d - A * cos³_i - B * cos_i

        @debug """
        Iteration #$it
          Estimation : $(abs(cos_i) <= 1 ? acosd(cos_i) : "INVALID") °
          Residue    : $(f₁) ° / day
        """

        # Compute the derivative.
        ∂f₁_∂cos_i = -3 * A * cos²_i - B

        # Update the estimate using the Newton-Raphson algorithm.
        cos_i = cos_i - f₁ / ∂f₁_∂cos_i

        # If the maximum number of iterations allowed has been reached, indicate that the
        # solution did not converged and exit loop.
        if (it >= max_iterations)
            @debug "The algorithm reached the maximum number of iterations without converging."
            converged = false
            break
        end

        it += 1
    end

    !converged && @warn("""
        The algorithm to compute the Sun-synchronous orbit inclination has not converged!
        Residue: $f₁ ° / day """
    )

    # If `cos_i` absolute value is larger than 1, the solution does not have physical
    # meaning.
    (abs(cos_i) > 1) &&
        throw(ArgumentError(
            "It is not possible to find a Sun-synchronous orbit with the selected parameters (a = $(a / 1000) km, e = $e)."
        ))

    # Recover the inclination.
    i = acos(cos_i)

    return i, converged
end

############################################################################################
#                                    Private Functions                                     #
############################################################################################

function _pretify_rev_per_days(i::Int, num::Int, den::Int)
    if num == 0
        return string(i)
    else
        return string(i) * " + " * pretty_number(String, num // den)
    end
end
