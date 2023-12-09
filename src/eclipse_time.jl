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

- `num_days::Number`: Number of days in which the analysis will be performed.
    (**Default** = 365)
- `step::Number`: The step in which the propagation will occur. Notice that this function
    has a crossing estimation to accurately estimate the transition between the regions.
    However, if this step is very large, we may miss some small regions. If it is negative,
    it will be selected as the time in which the mean anomaly advances 0.5°.
- `unit::Symbol`: Select the unit in which the results will be generated. The possible
    values are:
    - `:s` for seconds (**Default**);
    - `:m` for minutes; or
    - `:h` for hours.

# Returns

- `DataFrame`: The function returns a `DataFrame` with three columns:
    - `sunlight`: Total sunlight time per orbit at each day [`unit`].
    - `penumbra`: Total penumbra time per orbit at each day [`unit`].
    - `umbra`: Total umbra time per orbit at each day [`unit`].
    The unit of each column is stored in the `DataFrame` using metadata.

# Extended Help

## Examples

```
julia> orb = KeplerianElements(
           date_to_jd(2021, 1, 1, 0, 0, 0),
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

julia> df = eclipse_time_summary(orbp; num_days = 5, unit = :m)
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
  :penumbra => Dict("Unit"=>:m)
  :sunlight => Dict("Unit"=>:m)
  :umbra    => Dict("Unit"=>:m)
```
"""
function eclipse_time_summary(
    orbp::OrbitPropagator;
    num_days::Number = 365,
    step::Number = -1,
    unit::Symbol = :s
)
    jd₀ = Propagators.epoch(orbp)
    dt₀ = julian2datetime(jd₀)

    # TODO: Improve how the orbit period is computed.
    # We must obtain the mean elements to compute the orbit period. Maybe there is a better
    # way to do this.
    mean_elements = Propagators.mean_elements(orbp)

    # We need the orbit period because we will propagate one orbit per day.
    orb_period = orbital_period(mean_elements)

    # Check the propagation step we need to use.
    Δt₀ = step < 0 ? orb_period * 0.5 / 360 : Float64(step)

    # Vector of the days in which the eclipse time will be computed.
    days = 0:1:num_days-1

    # Pre-allocate the output variables.
    date          = Vector{DateTime}(undef, num_days)
    sunlight_time = zeros(num_days)
    penumbra_time = zeros(num_days)
    umbra_time    = zeros(num_days)

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
                # Refine to find the edge.
                t_k₀ = t_k - Δt
                t_k₁ = t_k
                t_kc = find_crossing(
                    _lighting_condition_crossing,
                    t_k₀,
                    t_k₁,
                    true,
                    false,
                    orbp,
                    d,
                    s_i,
                    old_state
                )

                # Times to be added in the old and new states.
                Δts₀ = t_kc - t_k₀
                Δts₁ = t_k₁ - t_kc

                _accumulate(
                    Δts₀,
                    old_state,
                    d + 1,
                    sunlight_time,
                    penumbra_time,
                    umbra_time
                )

                _accumulate(
                    Δts₁,
                    new_state,
                    d + 1,
                    sunlight_time,
                    penumbra_time,
                    umbra_time
                )

            # If not, just add the time step to the current state.
            else
                _accumulate(
                    Δt,
                    new_state,
                    d + 1,
                    sunlight_time,
                    penumbra_time,
                    umbra_time
                )
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
    if unit == :h
        sunlight_time ./= 3600
        penumbra_time ./= 3600
        umbra_time    ./= 3600
    elseif unit == :m
        sunlight_time ./= 60
        penumbra_time ./= 60
        umbra_time    ./= 60
    else
        # If the symbol is not known, we must use seconds.
        unit = :s
    end

    # Create and returns the DataFrame.
    df = DataFrame(
        date = Date.(date),
        sunlight = sunlight_time,
        penumbra = penumbra_time,
        umbra = umbra_time
    )

    # Add metadata to the DataFrame.
    metadata!(df, "Description", "Eclipse time PER ORBIT computed at each day.")
    colmetadata!(df, :sunlight, "Unit", unit)
    colmetadata!(df, :penumbra, "Unit", unit)
    colmetadata!(df, :umbra,    "Unit", unit)

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
    umbra_time::AbstractArray
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
    orbp::OrbitPropagator,
    t::Number,
    d::Number,
    s_i::AbstractVector
)
    r_i, ~ = Propagators.propagate!(orbp, 86400d + t)
    return lighting_condition(r_i, s_i)
end

# Function used in `find_crossing` to precisely obtain the instant in which the lightning
# condition changed.
@noinline function _lighting_condition_crossing(
    t::Number,
    orbp::OrbitPropagator,
    d::Number,
    s_i::AbstractVector,
    old_state::Symbol
)
    return _get_lighting_condition(orbp, t, d, s_i) == old_state
end
