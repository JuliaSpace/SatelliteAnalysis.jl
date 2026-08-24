## Description #############################################################################
#
# Progress interface shown during the decay analysis when the keyword `verbose` is enabled.
#
############################################################################################

# Inner width of the progress panel, i.e., the number of columns between the borders.
const _PROGRESS_PANEL_INNER_WIDTH = 52

# Width of the progress bar [columns].
const _PROGRESS_BAR_WIDTH = 38

# Number of terminal lines occupied by the progress panel.
const _PROGRESS_PANEL_LINES = 6

# Partial block characters used to smooth the progress bar, ordered by fill fraction.
const _PROGRESS_PARTIAL_BLOCKS = ("", "▏", "▎", "▍", "▌", "▋", "▊", "▉")

# Minimum wall time [s] between two panel redraws not caused by a progress advance.
const _PROGRESS_REDRAW_PERIOD = 0.02

"""
    mutable struct DecayProgress

Store the state of the progress interface shown during the decay analysis.

# Fields

- `io::IO`: Output stream of the interface.
- `ansi::Bool`: If `true`, the interface renders a live panel using ANSI escape sequences.
    Otherwise, it prints plain progress lines at every 10%.
- `tf::Float64`: Maximum propagation time [s] after the orbit epoch.
- `perigee₀::Float64`: Initial mean perigee altitude [m].
- `terminate_altitude::Float64`: Mean perigee altitude [m] that terminates the analysis.
- `start_wall::Float64`: Wall clock [s] at the beginning of the analysis.
- `last_draw_wall::Float64`: Wall clock [s] of the last panel redraw, used for throttling.
- `last_percent::Int`: Last progress percentage rendered, used to force a redraw whenever
    the progress advances.
- `drawn::Bool`: If `true`, the panel is currently drawn in the terminal.
"""
mutable struct DecayProgress
    io::IO
    ansi::Bool
    tf::Float64
    perigee₀::Float64
    terminate_altitude::Float64
    start_wall::Float64
    last_draw_wall::Float64
    last_percent::Int
    drawn::Bool
end

"""
    DecayProgress(io::IO, tf::Number, perigee₀::Number, terminate_altitude::Number; kwargs...) -> DecayProgress

Create the progress interface state writing to `io` for an analysis with the maximum
propagation time `tf` [s], the initial mean perigee altitude `perigee₀` [m], and the
terminate altitude `terminate_altitude` [m].

# Keywords

- `ansi::Bool`: If `true`, the interface renders a live ANSI panel. Otherwise, it prints
    plain progress lines.
    (**Default**: `io isa Base.TTY`)
"""
function DecayProgress(
    io::IO,
    tf::Number,
    perigee₀::Number,
    terminate_altitude::Number;
    ansi::Bool = io isa Base.TTY
)
    return DecayProgress(
        io,
        ansi,
        Float64(tf),
        Float64(perigee₀),
        Float64(terminate_altitude),
        0.0,
        0.0,
        -1,
        false
    )
end

############################################################################################
#                                     Public Interface                                     #
############################################################################################

"""
    struct DecayProgressCondition

Condition of the progress callback of the decay analysis.

# Fields

- `progress::Union{Nothing, DecayProgress}`: Progress interface state, or `nothing` when
    the interface is disabled.
"""
struct DecayProgressCondition
    progress::Union{Nothing, DecayProgress}
end

"""
    (c::DecayProgressCondition)(u, t, integrator) -> Bool

Return whether the progress interface of the condition `c` is enabled.
"""
function (c::DecayProgressCondition)(u, t, integrator)
    return !isnothing(c.progress)
end

"""
    struct DecayProgressAffect

Affect of the progress callback of the decay analysis, updating the progress interface at
every accepted step of the numerical integration.

# Fields

- `progress::Union{Nothing, DecayProgress}`: Progress interface state, or `nothing` when
    the interface is disabled.
"""
struct DecayProgressAffect
    progress::Union{Nothing, DecayProgress}
end

"""
    (a::DecayProgressAffect)(integrator) -> Nothing

Update the progress interface of the affect `a` with the current state of `integrator`.
"""
function (a::DecayProgressAffect)(integrator)
    progress = a.progress
    isnothing(progress) && return nothing

    ae, ee, _, _, _, _ = _equinoctial_to_classical(integrator.u)

    perigee = ae * (1 - ee) - EARTH_EQUATORIAL_RADIUS
    apogee  = ae * (1 + ee) - EARTH_EQUATORIAL_RADIUS

    _update_decay_progress!(progress, integrator.t, perigee, apogee)

    # The state was not modified, avoiding an unnecessary function re-evaluation.
    u_modified!(integrator, false)

    return nothing
end

"""
    _decay_progress_callback(progress::Union{Nothing, DecayProgress}) -> DiscreteCallback

Create the callback that updates the progress interface `progress` at every accepted step
of the numerical integration, or a disabled callback if `progress` is `nothing`. The
callback type does not depend on whether the interface is enabled, so the solver
specialization is shared by the verbose and silent paths. The callback does not modify
the state, does not save additional points, and hence does not change the analysis
result.
"""
function _decay_progress_callback(progress::Union{Nothing, DecayProgress})
    return DiscreteCallback(
        DecayProgressCondition(progress),
        DecayProgressAffect(progress);
        save_positions = (false, false)
    )
end

"""
    _start_decay_progress!(progress::DecayProgress, apogee₀::Number) -> Nothing

Start the progress interface `progress`, recording the wall clock and rendering the
initial state with the initial mean apogee altitude `apogee₀` [m]. In ANSI mode, the
terminal cursor is hidden until the interface finishes.
"""
function _start_decay_progress!(progress::DecayProgress, apogee₀::Number)
    progress.start_wall = time()

    progress.ansi && print(progress.io, "\e[?25l")

    _update_decay_progress!(progress, 0.0, progress.perigee₀, Float64(apogee₀))

    return nothing
end

"""
    _update_decay_progress!(progress::DecayProgress, t::Number, perigee::Number, apogee::Number) -> Nothing

Update the progress interface `progress` at the model time `t` [s] with the current mean
perigee and apogee altitudes [m]. In ANSI mode, the panel is redrawn whenever the progress
advances by at least one percentage point, and otherwise at most 50 times per second of
wall time. In the plain fallback, a line is printed when the progress crosses a multiple
of 5%.
"""
function _update_decay_progress!(
    progress::DecayProgress,
    t::Number,
    perigee::Number,
    apogee::Number
)
    fraction = _decay_progress_fraction(progress, t, perigee)

    if progress.ansi
        now = time()
        pct = floor(Int, fraction * 100)

        redraw = !progress.drawn ||
            (fraction >= 1) ||
            (pct > progress.last_percent) ||
            (now - progress.last_draw_wall >= _PROGRESS_REDRAW_PERIOD)

        redraw || return nothing

        _draw_decay_progress_panel(
            progress,
            fraction,
            t,
            perigee,
            apogee,
            now - progress.start_wall
        )

        progress.last_draw_wall = now
        progress.last_percent   = pct
        progress.drawn          = true
    else
        pct = 5 * floor(Int, fraction * 20)
        (pct <= progress.last_percent) && return nothing

        progress.last_percent = pct

        println(
            progress.io,
            "Decay analysis: $(pct)% (perigee ",
            _format_progress_altitude(perigee),
            ", model time ",
            _format_progress_span(t),
            ")"
        )
    end

    return nothing
end

"""
    _finish_decay_progress!(progress::DecayProgress, t_end::Number, perigee_end::Number, apogee_end::Number, reentered::Bool) -> Nothing

Finish the progress interface `progress` at the final model time `t_end` [s] with the
final mean perigee and apogee altitudes [m]. In ANSI mode, the panel is redrawn with the
final state and left on the screen. A summary line is then printed with the outcome: a
reentry after `t_end` if `reentered` is `true`, or no reentry within the maximum
propagation time otherwise.
"""
function _finish_decay_progress!(
    progress::DecayProgress,
    t_end::Number,
    perigee_end::Number,
    apogee_end::Number,
    reentered::Bool
)
    wall = time() - progress.start_wall

    msg = reentered ?
        "Decay analysis: reentry after " * _format_progress_span(t_end) :
        "Decay analysis: no reentry within " * _format_progress_span(progress.tf)

    msg *= " (wall time: " * _format_progress_span(wall) * ")"

    if progress.ansi
        # Redraw the panel with the final state, leaving it on the screen, and print the
        # summary below it. Notice that the progress fraction is 1 at either termination
        # condition by construction.
        fraction = _decay_progress_fraction(progress, t_end, perigee_end)

        _draw_decay_progress_panel(progress, fraction, t_end, perigee_end, apogee_end, wall)
        progress.drawn = true

        println(progress.io, "\e[32m✓\e[0m ", msg)
    else
        println(progress.io, msg)
    end

    return nothing
end

"""
    _cleanup_decay_progress!(progress::DecayProgress) -> Nothing

Restore the terminal cursor hidden by the progress interface `progress`. This function
must be called even if the numerical integration throws.
"""
function _cleanup_decay_progress!(progress::DecayProgress)
    progress.ansi && print(progress.io, "\e[?25h")
    return nothing
end

############################################################################################
#                                    Private Functions                                     #
############################################################################################

"""
    _decay_progress_fraction(progress::DecayProgress, t::Number, perigee::Number) -> Float64

Compute the progress fraction at the model time `t` [s] with the current mean perigee
altitude `perigee` [m]. The fraction is the maximum between the time fraction and the
perigee descent fraction, clamped to [0, 1]. Hence, it reaches 100% at either termination
condition of the analysis: the perigee reaching the terminate altitude or the time
reaching the maximum propagation time.
"""
function _decay_progress_fraction(progress::DecayProgress, t::Number, perigee::Number)
    time_fraction = t / progress.tf

    Δh = progress.perigee₀ - progress.terminate_altitude
    altitude_fraction = Δh > 0 ? (progress.perigee₀ - perigee) / Δh : 1.0

    return clamp(max(time_fraction, altitude_fraction), 0.0, 1.0)
end

"""
    _draw_decay_progress_panel(progress::DecayProgress, fraction::Number, t::Number, perigee::Number, apogee::Number, wall::Number) -> Nothing

Draw the ANSI progress panel of `progress` with the progress `fraction`, the model time
`t` [s], the mean perigee and apogee altitudes [m], and the elapsed wall time `wall` [s].
On redraws, the panel is rewritten in place.
"""
function _draw_decay_progress_panel(
    progress::DecayProgress,
    fraction::Number,
    t::Number,
    perigee::Number,
    apogee::Number,
    wall::Number
)
    io    = progress.io
    inner = _PROGRESS_PANEL_INNER_WIDTH

    # Move the cursor to the panel origin on redraws.
    progress.drawn && print(io, "\e[$(_PROGRESS_PANEL_LINES)A")

    # -- Title Border ----------------------------------------------------------------------

    print(io, "\e[2K\e[90m╭─ \e[0m\e[1;36mDECAY ANALYSIS\e[0m\e[90m ")
    print(io, "─"^(inner - 17), "╮\e[0m\n")

    # -- Progress Bar ----------------------------------------------------------------------

    pct = lpad(string(round(fraction * 100; digits = 1)) * "%", 6)
    _print_panel_row(
        io,
        " " * _render_progress_bar(fraction) * "  \e[1m" * pct * "\e[0m",
        1 + _PROGRESS_BAR_WIDTH + 2 + 6
    )

    # -- Separator and Data Rows -----------------------------------------------------------

    _print_panel_row(io, "", 0)

    perigee_str = _format_progress_altitude(perigee)
    apogee_str  = _format_progress_altitude(apogee)
    model_str   = _format_progress_span(t)
    wall_str    = _format_progress_span(wall)

    _print_panel_row(
        io,
        " \e[90m" * rpad("Perigee", 12) * "\e[0m" * rpad(perigee_str, 13) * "\e[90m" *
            rpad("Apogee", 12) * "\e[0m" * apogee_str,
        1 + 12 + 13 + 12 + length(apogee_str)
    )

    _print_panel_row(
        io,
        " \e[90m" * rpad("Model time", 12) * "\e[0m" * rpad(model_str, 13) * "\e[90m" *
            rpad("Elapsed", 12) * "\e[0m" * wall_str,
        1 + 12 + 13 + 12 + length(wall_str)
    )

    # -- Bottom Border ---------------------------------------------------------------------

    print(io, "\e[2K\e[90m╰", "─"^inner, "╯\e[0m\n")

    return nothing
end

"""
    _print_panel_row(io::IO, content::String, visible_length::Int) -> Nothing

Print one row of the progress panel to `io` with the borders around `content`, padding it
to the panel inner width. Since `content` may carry ANSI escape sequences, its visible
length must be provided in `visible_length`.
"""
function _print_panel_row(io::IO, content::String, visible_length::Int)
    padding = max(_PROGRESS_PANEL_INNER_WIDTH - visible_length, 0)
    print(io, "\e[2K\e[90m│\e[0m", content, " "^padding, "\e[90m│\e[0m\n")
    return nothing
end

"""
    _render_progress_bar(fraction::Number) -> String

Render the progress bar for the progress `fraction`, using partial block characters to
smooth the leading edge. The returned string carries ANSI color sequences.
"""
function _render_progress_bar(fraction::Number)
    cells   = clamp(fraction, 0, 1) * _PROGRESS_BAR_WIDTH
    full    = floor(Int, cells)
    partial = _PROGRESS_PARTIAL_BLOCKS[floor(Int, (cells - full) * 8) + 1]
    empty   = _PROGRESS_BAR_WIDTH - full - (isempty(partial) ? 0 : 1)

    return "\e[36m" * "█"^full * partial * "\e[90m" * "░"^empty * "\e[0m"
end

"""
    _format_progress_altitude(h::Number) -> String

Format the altitude `h` [m] in kilometers with one decimal digit.
"""
function _format_progress_altitude(h::Number)
    return string(round(h / 1000; digits = 1)) * " km"
end

"""
    _format_progress_span(seconds::Number) -> String

Format the time span `seconds` [s] using seconds, minutes, hours, days, or Julian years,
selecting the unit that best fits its magnitude.
"""
function _format_progress_span(seconds::Number)
    s = Float64(seconds)

    if s < 120
        return string(round(s; digits = 1)) * " s"
    elseif s < 2 * 3600
        return string(round(s / 60; digits = 1)) * " min"
    elseif s < 2 * 86400
        return string(round(s / 3600; digits = 1)) * " h"
    elseif s < 2 * 365.25 * 86400
        return string(round(s / 86400; digits = 1)) * " days"
    else
        return string(round(s / (365.25 * 86400); digits = 2)) * " years"
    end
end
