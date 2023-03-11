# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==============================================================================
#
#   Functions to design Sun synchronous orbits.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export design_sun_sync_ground_repeating_orbit

"""
    design_sun_sync_ground_repeating_orbit(minimum_repetition::Int, maximum_repetition::Int; kwargs...)

List all the Sun synchronous, ground repeating orbits in which their repetition
period is in the interval `[minimum_repetition, maximum_repetition]` days.

This function returns a `DataFrame` with the following columns:

- `semi_major_axis`: The orbit semi-major axis.
- `altitude`: The orbit altitude above the Equator `(a - R0)`.
- `inclination`: The orbit inclination.
- `period`: The orbital period.
- `rev_per_days`: A `Tuple` with the integer and rational parts of the number of
    revolutions per day.

!!! note
    The units of those columns depends on the keywords.

# Keywords

- `R0::Number`: Earth's equatorial radius for the analsyis [m].
    (**Default** = `R0`)
- `angle_unit::Symbol`: Unit for all the angles in the output `DataFrame`.
    It can be `:deg` for degrees or `:rad` for radians. (**Default**: `:deg`)
- `distance_unit::Symbol`: The unit for all the distances in the output
    `DataFrame`. It can be `:m` for meters or `:km` for kilometers.
    (**Default**: `:km`)
- `e::Number`: Orbit eccentricity. (**Default**: 0)
- `int_rev_per_day::Tuple`: `Tuple` with the integer parts of the number of
    revolutions per day to be analyzed. (**Default** = `(13, 14, 15, 16, 17)`)
- `maximum_altitude::Union{Nothing, Number}`: Maximum altitude [m] of the orbits
    in the output `DataFrame`. If it is `nothing`, the algorithm will not apply
    a higher limit to the orbital altitude. (**Default** = `nothing`)
- `minimum_altitude::Union{Nothing, Number}`: Minimum altitude [m] of the orbits
    in the output `DataFrame`. If it is `nothing`, the algorithm will not apply
    a lower limit to the orbital altitude. (**Default** = `nothing`)
- `time_unit::Symbol`: Unit for all the time values in the output `DataFrame`.
    It can be `:s` for seconds, `:m` for minutes, or `:h` for hours.
    (**Default** = `:h`)
"""
function design_sun_sync_ground_repeating_orbit(
    minimum_repetition::Int,
    maximum_repetition::Int;
    R0::Number = R0,
    angle_unit::Symbol = :deg,
    distance_unit::Symbol = :km,
    e::Number = 0,
    int_rev_per_day::Tuple = (13, 14, 15, 16, 17),
    maximum_altitude::Union{Nothing, Number} = nothing,
    minimum_altitude::Union{Nothing, Number} = nothing,
    time_unit::Symbol = :m
)
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
        semi_major_axis = Float64[],
        altitude = Float64[],
        inclination = Float64[],
        period = Float64[],
        rev_per_days = Tuple{Int, Rational}[]
    )

    # Check the units for the values.
    dunit   = distance_unit === :km  ? 1.e-3     : 1.0
    angunit = angle_unit    === :deg ? 180.0 / π : 1.0
    tunit   = begin
        if time_unit === :h
            1.0 / 3600
        elseif time_unit === :m
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
                # Compute the number of revolutions per day of this orbit, and
                # convert it to angular velocity.
                num_rev_per_day = int + num / den
                n = num_rev_per_day * 2π / 86400

                # Find a Sun synchronous orbit with that angular velocity.
                a, i, converged = sun_sync_orbit_from_ang_vel(n, e)

                # If the algorithm has not converged or if the orbit is not
                # valid, skip this value.
                (!converged || !@check_orbit(a, e)) && continue

                # If we reach this point, add the orbit to the `DataFrame`.
                orbit_period = period(a, e, i, :J2)

                push!(df, (
                    a * dunit,
                    (a - R0) * dunit,
                    i * angunit,
                    orbit_period * tunit,
                    (int, num // den)
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
