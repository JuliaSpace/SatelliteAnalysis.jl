# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==============================================================================
#
#   Compute the satellite eclipse time summary.
#
# References
# ==============================================================================
#
# [1] Longo, C. R. O., Rickman, S. L (1995). Method for the Calculation of
#     Spacecraft Umbra and Penumbra Shadow Terminator Points. NASA Technical
#     Paper 3547.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export eclipse_time_summary

"""
    eclipse_time_summary([io::IO,] orbp::OrbitPropagator; kwargs...)

Compute the eclipse time summary for the orbit propagator `orbp`. The summary is
computed as the total time the object stays in the sunlight, penumbra, and umbra
regions per day.

# Returns

If the `io` parameter **is not** present, then the function returns three
vectors:

- The total sunlight time per day [s].
- The total penumbra time per day [s].
- The total umbra time per day [s].

If the `io` parameter **is** present, then the function prints to `io` the
results in table format.

# Keywords

- `num_days::Number`: Number of days in which the analysis will be performed.
    (**Default** = 365)
- `step::Number`: The step in which the propagation will occur. Notice that this
    function has a crossing estimation to accurately estimate the transition
    between the regions. However, if this step is very large, we may miss some
    small regions. If it is negative, then it will be selected as the time in
    which the mean anomaly advances 0.5°.

If the parameter `io` is passed, then the following additional keywords are
available:

- `unit::Symbol`: Select the unit in which the results will be printed. The
    possible values are:
    - `:s` for seconds (**Default**);
    - `:m` for minutes; or
    - `:h` for hours.
- `use_pager::Bool`: If `true`, then the result will be printed using
    **TerminalPager.jl**.
"""
function eclipse_time_summary(
    orbp::OrbitPropagator;
    num_days::Number = 365,
    step::Number = -1
)
    jd₀ = get_epoch(orbp)

    # TODO: Improve how the orbit period is computed.
    # We must obtain the mean elements to compute the orbit period. Maybe there
    # is a better way to do this.
    mean_elements = get_mean_elements(orbp)

    # We need the orbit period because we will propagate one orbit per day.
    orbit_period = period(mean_elements)

    # Check the propagation step we need to use.
    Δt₀ = step < 0 ? orbit_period * 0.5 / 360 : Float64(step)

    # Vector of the days in which the eclipse time will be computed.
    days = 0:1:num_days-1

    # Preallocate the output variables.
    sunlight_time = zeros(num_days)
    penumbra_time = zeros(num_days)
    umbra_time    = zeros(num_days)

    # Loop
    # ==========================================================================

    @inbounds for d in days
        # TODO: Should we transform between MOD => TOD/TEME?
        # Get the Sun position represented in the inertial reference frame.
        s_i = sun_position_i(jd₀ + d)

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
                    _fc,
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

            abs(t_k - orbit_period) < 1e-3 && break

            # Make sure that the last interval will have the exact size so that
            # the end of the analysis is the end of the orbit.
            if (t_k + Δt > orbit_period)
                Δt = orbit_period - t_k
            end

            t_k += Δt
        end
    end

    return sunlight_time, penumbra_time, umbra_time
end

function eclipse_time_summary(
    io::IO,
    orbp::OrbitPropagator;
    num_days::Number = 365,
    step::Number = -1,
    unit::Symbol = :s,
    use_pager::Bool = false
)
    # Compute the summary.
    sunlight_time, penumbra_time, umbra_time = eclipse_time_summary(
        orbp;
        num_days,
        step
    )

    # Get the date vector.
    vjd = get_epoch(orbp) .+ collect(0:1:num_days-1)
    vdate = string.(jd_to_date.(DateTime, vjd))

    # Obtain the conversion given the unit selection.
    if unit == :h
        unit_label = "h"
        unit_factor = 3600
    elseif unit == :m
        unit_label = "min"
        unit_factor = 60
    else
        unit_label = "s"
        unit_factor = 1
    end

    sunlight_time ./= unit_factor
    penumbra_time ./= unit_factor
    umbra_time ./= unit_factor

    str = pretty_table(
        String,
        vcat(
            hcat(vdate, sunlight_time, penumbra_time, umbra_time),
            hcat(
                AnsiTextCell(string(crayon"red bold") * "Maximum"),
                maximum(sunlight_time),
                maximum(penumbra_time),
                maximum(umbra_time)
            ),
            hcat(
                AnsiTextCell(string(crayon"yellow bold") * "Mean"),
                mean(sunlight_time),
                mean(penumbra_time),
                mean(umbra_time)
            ),
            hcat(
                AnsiTextCell(string(crayon"blue bold") * "Minimum"),
                minimum(sunlight_time),
                minimum(penumbra_time),
                minimum(umbra_time)
            ),
        );
        color = get(stdout, :color, false),
        crop = use_pager ? :none : :horizontal,
        display_size = use_pager ? (-1, -1) : displaysize(stdout),
        header = [
            "Day",
            "Sunlight time [$unit_label]",
            "Penumbra time [$unit_label]",
            "Umbra time [$unit_label]"
        ],
        hlines = [:header, num_days + 1],
        alignment_anchor_regex = Dict(
            2 => [r"\."],
            3 => [r"\."],
            4 => [r"\."]
        ),
        vlines = [1]
    )

    # Check if we need to use the pager.
    if use_pager
        pager(
            str;
            freeze_columns = max(length.(vdate)..., 6) + 3,
            freeze_rows = 2
        )
    else
        println(str)
    end

    return nothing
end


#                              Private functions
# ==============================================================================

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

# Without the `@noinline` we get **a lot** of allocations when calling the
# function `find_crossing`.
@noinline function _get_lighting_condition(
    orbp::OrbitPropagator,
    t::Number,
    d::Number,
    s_i::AbstractVector
)
    r_i, ~ = propagate!(orbp, 86400d + t)
    return lighting_condition(r_i, s_i)
end

@noinline function _fc(
    t::Number,
    orbp::OrbitPropagator,
    d::Number,
    s_i::AbstractVector,
    old_state::Symbol
)
    return _get_lighting_condition(orbp, t, d, s_i) == old_state
end
