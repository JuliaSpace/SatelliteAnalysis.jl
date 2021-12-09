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
