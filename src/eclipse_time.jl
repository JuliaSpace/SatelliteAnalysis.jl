## Description #############################################################################
#
# Compute the satellite eclipse time summary.
#
## References ##############################################################################
#
# [1] Longo, C. R. O., Rickman, S. L (1995). Method for the Calculation of Spacecraft Umbra
#     and Penumbra Shadow Terminator Points. NASA Technical Paper 3547.
#
############################################################################################

export eclipse_time_summary

"""
    eclipse_time_summary(orbp::OrbitPropagator; kwargs...) -> DataFrame

Compute the eclipse time summary for the orbit propagator `orbp`. The summary is computed as
the total time the object stays in the sunlight, penumbra, and umbra regions per orbit at
each day.

# Keywords

- `num_days::Integer`: Number of days in which the analysis will be performed.
    (**Default** = 365)
- `step::Union{Nothing, Number}`: The step [s] in which the propagation will occur. Notice
    that this function has a crossing estimation to accurately estimate the transition
    between the regions, including those entirely inside one step, such as a penumbra
    passage between the sunlight and the umbra. However, if this step is very large, we may
    miss a region if the lighting condition is the same in two consecutive instants. If it
    is `nothing`, it will be selected as the time in which the mean anomaly advances 0.5°.
    (**Default** = `nothing`)
- `time_unit::Symbol`: Select the unit in which the results will be generated. The possible
    values are:
    - `:s` for seconds (**Default**);
    - `:min` for minutes; or
    - `:h` for hours.

# Returns

- `DataFrame`: The function returns a `DataFrame` with four columns:
    - `date`: Date of the analysis [UTC] encoded using `Date`.
    - `sunlight`: Total sunlight time per orbit at each day [`time_unit`].
    - `penumbra`: Total penumbra time per orbit at each day [`time_unit`].
    - `umbra`: Total umbra time per orbit at each day [`time_unit`].
    The unit of each column is stored in the `DataFrame` using metadata.

# Extended Help

## Throws

- `ArgumentError`: If `num_days` is lower than 1, if `step` is not positive or not lower
    than the orbital period, or if `time_unit` is not `:s`, `:min`, or `:h`.

## Examples

```julia-repl
julia> using SatelliteAnalysis

julia> jd₀ = date_to_jd(2021, 1, 1, 0, 0, 0)

julia> orb = KeplerianElements(
           jd₀,
           7130.982e3,
           0.001111,
           98.405 |> deg2rad,
           ltdn_to_raan(10.5, jd₀),
           90     |> deg2rad,
           0
       )

julia> orbp = Propagators.init(Val(:J2), orb)

julia> df = eclipse_time_summary(orbp; num_days = 5)
5×4 DataFrame
 Row │ date        sunlight  penumbra  umbra
     │ Date        Float64   Float64   Float64
─────┼─────────────────────────────────────────
   1 │ 2021-01-01   3972.63   20.4117  2006.96
   2 │ 2021-01-02   3973.85   20.4376  2005.71
   3 │ 2021-01-03   3974.77   20.4575  2004.77
   4 │ 2021-01-04   3975.74   20.4758  2003.79
   5 │ 2021-01-05   3976.94   20.5022  2002.55

julia> df = eclipse_time_summary(orbp; num_days = 5, time_unit = :min)
5×4 DataFrame
 Row │ date        sunlight  penumbra  umbra
     │ Date        Float64   Float64   Float64
─────┼─────────────────────────────────────────
   1 │ 2021-01-01   66.2105  0.340195  33.4493
   2 │ 2021-01-02   66.2308  0.340627  33.4285
   3 │ 2021-01-03   66.2461  0.340958  33.4129
   4 │ 2021-01-04   66.2623  0.341263  33.3964
   5 │ 2021-01-05   66.2824  0.341704  33.3759

julia> colmetadata(df)
Dict{Symbol, Dict{String, Symbol}} with 3 entries:
  :penumbra => Dict("Unit"=>:min)
  :sunlight => Dict("Unit"=>:min)
  :umbra    => Dict("Unit"=>:min)
```
"""
function eclipse_time_summary(
    orbp::OrbitPropagator;
    num_days::Integer = 365,
    step::Union{Nothing, Number} = nothing,
    time_unit::Symbol = :s,
)
    num_days < 1 && throw(ArgumentError("The number of days must be greater than 0."))

    # Factor to convert the time from seconds to the selected unit. Notice that this
    # function also validates the input.
    time_factor = _time_unit_factor(time_unit)

    jd₀ = Propagators.epoch(orbp)
    dt₀ = julian2datetime(jd₀)

    # TODO: Improve how the orbit period is computed.
    # We must obtain the mean elements to compute the orbit period. Maybe there is a better
    # way to do this.
    mean_elements = Propagators.mean_elements(orbp)

    # If the propagator cannot return the mean elements, we will compute the orbit period by
    # converting the osculating elements to Keplerian elements.
    if isnothing(mean_elements)
        r_i, v_i = Propagators.propagate!(orbp, 0)
        mean_elements = rv_to_kepler(r_i, v_i, jd₀)
    end

    # We need the orbit period because we will propagate one orbit per day.
    orb_period = orbital_period(mean_elements)

    # Check the propagation step we need to use. If the user did not specify the step, we
    # select the time in which the mean anomaly advances 0.5°.
    step′ = isnothing(step) ? orb_period * (one(orb_period) / 2) / 360 : step

    (0 < step′ < orb_period) || throw(
        ArgumentError(
            "The step must be positive and lower than the orbital period ($orb_period s).",
        ),
    )

    time_type = promote_type(typeof(orb_period), typeof(step′))
    Δt₀ = time_type(step′)

    # Vector of the days in which the eclipse time will be computed.
    days = 0:1:(num_days - 1)

    # Pre-allocate the output variables.
    date          = Vector{DateTime}(undef, num_days)
    sunlight_time = zeros(time_type, num_days)
    penumbra_time = zeros(time_type, num_days)
    umbra_time    = zeros(time_type, num_days)

    # == Loop ==============================================================================

    @inbounds for d in days
        date[d + 1] = dt₀ + Day(d)

        # TODO: Should we transform between MOD => TOD/TEME?
        # Get the Sun position represented in the inertial reference frame.
        s_i = sun_position_mod(jd₀ + d)

        # Initial state.
        old_state = _get_lighting_condition(orbp, 0, d, s_i)

        # Compute the eclipse time during one orbit.
        Δt  = Δt₀
        t_k = Δt

        while true
            new_state = _get_lighting_condition(orbp, t_k, d, s_i)

            # Check if the state has changed.
            if new_state != old_state
                # Refine to find the edges. Notice that we can have more than one edge
                # inside the step. For example, the satellite can cross the entire penumbra
                # region when going from the sunlight to the umbra. Hence, after finding an
                # edge, we must verify the state just after it and keep searching until we
                # reach the state at the end of the step.
                t_e     = t_k - Δt
                state_e = old_state

                while state_e != new_state
                    t_kc = find_crossing(
                        _lighting_condition_crossing,
                        t_e,
                        t_k,
                        true,
                        false,
                        orbp,
                        d,
                        s_i,
                        state_e,
                    )

                    _accumulate(
                        t_kc - t_e, state_e, d + 1, sunlight_time, penumbra_time, umbra_time
                    )

                    t_e     = t_kc
                    state_e = _get_lighting_condition(orbp, t_kc, d, s_i)
                end

                # The remaining time in the step must be added to the new state.
                _accumulate(
                    t_k - t_e, new_state, d + 1, sunlight_time, penumbra_time, umbra_time
                )

                # If not, just add the time step to the current state.
            else
                _accumulate(Δt, new_state, d + 1, sunlight_time, penumbra_time, umbra_time)
            end

            old_state = new_state

            abs(t_k - orb_period) < 1e-3 && break

            # Make sure that the last interval will have the exact size so that the end of
            # the analysis is the end of the orbit.
            if (t_k + Δt > orb_period)
                Δt = orb_period - t_k
            end

            t_k += Δt
        end
    end

    # Convert to the right units.
    sunlight_time .*= time_factor
    penumbra_time .*= time_factor
    umbra_time    .*= time_factor

    # Create and returns the DataFrame.
    df = DataFrame(;
        date = Date.(date),
        sunlight = sunlight_time,
        penumbra = penumbra_time,
        umbra = umbra_time,
    )

    # Add metadata to the DataFrame. The style `:note` makes the metadata propagate through
    # DataFrame transformations.
    metadata!(
        df, "Description", "Eclipse time PER ORBIT computed at each day."; style = :note
    )

    colmetadata!(df, :sunlight, "Unit", time_unit; style = :note)
    colmetadata!(df, :penumbra, "Unit", time_unit; style = :note)
    colmetadata!(df, :umbra, "Unit", time_unit; style = :note)

    return df
end

############################################################################################
#                                    Private Functions                                     #
############################################################################################

# Accumulate the time in a specific lighting state.
function _accumulate(
    Δts::Number,
    state::Symbol,
    ind::Integer,
    sunlight_time::AbstractArray,
    penumbra_time::AbstractArray,
    umbra_time::AbstractArray,
)
    @inbounds if state == :sunlight
        sunlight_time[ind] += Δts
    elseif state == :penumbra
        penumbra_time[ind] += Δts
    elseif state == :umbra
        umbra_time[ind] += Δts
    end

    return nothing
end

# Without the `@noinline` we get **a lot** of allocations when calling the function
# `find_crossing`.
@noinline function _get_lighting_condition(
    orbp::OrbitPropagator, t::Number, d::Number, s_i::AbstractVector
)
    r_i, ~ = Propagators.propagate!(orbp, 86400d + t)
    return lighting_condition(r_i, s_i)
end

# Function used in `find_crossing` to precisely obtain the instant in which the lighting
# condition changed.
@noinline function _lighting_condition_crossing(
    t::Number, orbp::OrbitPropagator, d::Number, s_i::AbstractVector, old_state::Symbol
)
    return _get_lighting_condition(orbp, t, d, s_i) == old_state
end
