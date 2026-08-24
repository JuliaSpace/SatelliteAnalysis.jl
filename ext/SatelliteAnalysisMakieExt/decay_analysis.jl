## Description #############################################################################
#
# Function to plot the decay analysis.
#
############################################################################################

############################################################################################
#                                     Public Functions                                     #
############################################################################################

function SatelliteAnalysis.plot_decay_analysis(
    df::DataFrame;
    fontscale::Real = 1,
    mission_name::Union{Nothing, String}        = nothing,
    mono_ticklabels::Bool = false,
    panel_width::Union{Nothing, Int} = nothing,
    satellite_mass::Union{Nothing, Number}      = nothing,
    satellite_mean_area::Union{Nothing, Number} = nothing,
    show_assumptions::Bool = true,
    show_f107::Bool = false,
    show_reentry_callout::Bool = true,
    show_reentry_date::Bool = false,
    subtitle::Union{Nothing, String, Symbol} = :auto,
    terminate_altitude::Union{Nothing, Number}  = nothing,
    theme::Symbol = :light,
    title::String = "Orbital Decay Analysis",
    xlims::Union{Nothing, Tuple} = nothing,
    ylims::Union{Nothing, Tuple} = nothing,
    size = (1280, 720),
    kwargs...
)
    # == Input Validation ==================================================================

    required_columns = (:time, :date, :apogee_altitude, :perigee_altitude)

    for c in (show_f107 ? (required_columns..., :f107) : required_columns)
        hasproperty(df, c) || throw(
            ArgumentError(
                "The input `DataFrame` must have the column `$c`. It should be obtained " *
                "using the function `decay_analysis`."
            )
        )
    end

    isempty(df) && throw(ArgumentError("The input `DataFrame` is empty."))

    (subtitle isa Symbol && subtitle != :auto) && throw(
        ArgumentError("The keyword `subtitle` must be a `String`, `nothing`, or `:auto`.")
    )

    # Build the theme first since it also validates the variant in `theme`.
    sa_theme = makie_theme(theme; fontscale = fontscale, mono_ticklabels = mono_ticklabels)

    # Width of the column with the information panel and the legend, scaling with the
    # figure width by default.
    panel_width′ = isnothing(panel_width) ?
        clamp(round(Int, 0.22 * size[1]), 240, 360) :
        panel_width

    # == Information Panel Values ==========================================================

    # The keywords have precedence over the `DataFrame` metadata.
    mass = isnothing(satellite_mass) ?
        metadata(df, "Satellite Mass", nothing) :
        satellite_mass

    area = isnothing(satellite_mean_area) ?
        metadata(df, "Satellite Mean Area", nothing) :
        satellite_mean_area

    term = isnothing(terminate_altitude) ?
        metadata(df, "Terminate Altitude", nothing) :
        terminate_altitude

    # Assumption values, resolved only from the `DataFrame` metadata.
    atm_model_name = metadata(df, "Atmospheric Model", nothing)
    C_d            = metadata(df, "Drag Coefficient",  nothing)
    C_r            = metadata(df, "SRP Coefficient",   nothing)
    f107_source    = metadata(df, "F10.7 Source",      nothing)

    # == Units =============================================================================

    time_unit     = colmetadata(df, :time,            "Unit", :y)
    distance_unit = colmetadata(df, :apogee_altitude, "Unit", :km)

    # Human-readable axis labels, falling back to the raw symbol if it is unknown.
    time_unit_label     = get(_TIME_UNIT_LABELS,     time_unit,     string(time_unit))
    distance_unit_label = get(_DISTANCE_UNIT_LABELS, distance_unit, string(distance_unit))

    # == Reentry Detection =================================================================

    reentered = false

    if !isnothing(term)
        # The terminate altitude is stored in meters, whereas the altitude columns use
        # `distance_unit`, which is either meters or kilometers.
        term_alt = distance_unit == :m ? term : term / 1000

        # The numerical integration stops very close to the terminate altitude, but not
        # exactly at it. Hence, we need a small tolerance to detect the reentry.
        reentered = last(df.perigee_altitude) <= term_alt * (1 + 1e-6)
    end

    # == Colors ============================================================================

    dark = theme == :dark

    accent_color    = dark ? MAGENTA_DARK        : MAGENTA_LIGHT
    border_color    = dark ? BORDER_DARK         : BORDER_LIGHT
    card_color      = dark ? NAVY_CARD           : SURFACE_CARD
    f107_color      = dark ? CATEGORICAL_DARK[3] : CATEGORICAL_LIGHT[3]
    subline_color   = dark ? TEXT_TERTIARY_DARK  : TEXT_TERTIARY_LIGHT
    threshold_color = dark ? CATEGORICAL_DARK[6] : CATEGORICAL_LIGHT[6]
    title_color     = dark ? CYAN_DARK           : CYAN_LIGHT

    # == Assemble the Information Panel Cards =============================================

    # Each card has a title, its value lines, and a flag selecting the compact rendering,
    # in which all the lines use a small font size.
    cards = Tuple{String, Vector{String}, Bool}[]

    !isnothing(mass) &&
        push!(cards, ("SATELLITE MASS", ["$(_format_number(mass)) kg"], false))
    !isnothing(area) &&
        push!(cards, ("MEAN AREA", ["$(_format_number(area)) m²"], false))

    if show_assumptions && !isnothing(mass) && !isnothing(area) && !isnothing(C_d)
        # Drag-facing ballistic coefficient convention, as used by STELA.
        bc = C_d * area / mass
        push!(cards, ("BALLISTIC COEFF.", ["$(_format_number(bc)) m²/kg"], false))
    end

    if !isnothing(term)
        if reentered
            value_lines = [_format_duration(last(df.date) - first(df.date))]

            if show_reentry_date
                reentry_epoch = Dates.format(last(df.date), dateformat"yyyy-mm-dd HH:MM")
                push!(value_lines, reentry_epoch * " UTC")
            end

            push!(cards, ("TIME TO REENTER", value_lines, false))
        else
            push!(cards, ("TIME TO REENTER", ["No reentry"], false))
        end
    end

    if show_assumptions
        assumption_lines = String[]

        !isnothing(atm_model_name) &&
            push!(assumption_lines, "Atm. model: " * string(atm_model_name))
        !isnothing(f107_source) && push!(assumption_lines, "F10.7: " * string(f107_source))

        if !isnothing(C_d) && !isnothing(C_r)
            push!(
                assumption_lines,
                "C_d $(_format_number(C_d)) · C_r $(_format_number(C_r))"
            )
        elseif !isnothing(C_d)
            push!(assumption_lines, "C_d $(_format_number(C_d))")
        elseif !isnothing(C_r)
            push!(assumption_lines, "C_r $(_format_number(C_r))")
        end

        !isempty(assumption_lines) && push!(cards, ("ASSUMPTIONS", assumption_lines, true))
    end

    # == Plot ==============================================================================

    # Every object must be created inside `with_theme` because Makie resolves the theme
    # attributes at object-creation time.
    return with_theme(sa_theme) do
        fig = Figure(; size = size, kwargs...)

        legend_plots  = Any[]
        legend_labels = String[]

        # == F10.7 Twin Y Axis =============================================================

        # The twin axis must be created before the main axis so that the F10.7 line is
        # rendered below all the main axis elements. Hence, it keeps the default axis
        # background, whereas the main axis background is made transparent.
        ax_f107   = nothing
        f107_line = nothing

        if show_f107
            ax_f107 = Axis(
                fig[1, 1];
                yaxisposition     = :right,
                ygridvisible      = false,
                ylabel            = "F10.7 [sfu]",
                yminorgridvisible = false,
            )

            hidexdecorations!(ax_f107)

            f107_line = lines!(
                ax_f107,
                df.time,
                df.f107;
                color     = (f107_color, 0.5),
                linewidth = 1.5,
            )
        end

        # == Main Axis =====================================================================

        main_axis_background = show_f107 ? (; backgroundcolor = :transparent) : (;)

        # The automatic subtitle shows the analysis timespan.
        subtitle_str = subtitle isa Symbol ?
            Dates.format(first(df.date), dateformat"yyyy-mm-dd") * " → " *
                Dates.format(last(df.date), dateformat"yyyy-mm-dd") * " UTC" :
            subtitle

        subtitle_attrs = isnothing(subtitle_str) ? (;) : (;
            subtitle      = subtitle_str,
            subtitlecolor = subline_color,
            subtitlesize  = 16.0 * fontscale,
        )

        ax = Axis(
            fig[1, 1];
            title  = title,
            xlabel = "Time [$time_unit_label]",
            ylabel = "Altitude [$distance_unit_label]",
            subtitle_attrs...,
            main_axis_background...
        )

        push!(legend_plots, lines!(ax, df.time, df.apogee_altitude))
        push!(legend_labels, "Apogee Altitude")

        push!(legend_plots, lines!(ax, df.time, df.perigee_altitude))
        push!(legend_labels, "Perigee Altitude")

        # == Terminate Altitude Threshold ==================================================

        if !isnothing(term)
            threshold_line = hlines!(
                ax,
                [term_alt];
                color     = (threshold_color, 0.8),
                linestyle = :dash,
                linewidth = 1.5,
            )

            # Render the threshold below the altitude lines.
            translate!(threshold_line, 0, 0, -1)

            push!(legend_plots, threshold_line)
            push!(legend_labels, "Reentry Altitude")
        end

        # == Reentry Marker and Callout ====================================================

        if reentered
            reentry_marker = scatter!(
                ax,
                [last(df.time)],
                [last(df.perigee_altitude)];
                color      = accent_color,
                markersize = 14,
            )
            translate!(reentry_marker, 0, 0, 10)

            push!(legend_plots, reentry_marker)
            push!(legend_labels, "Reentry")

            if show_reentry_callout
                callout_text = "Reentry: " *
                    Dates.format(last(df.date), dateformat"yyyy-mm-dd")

                # The pixel offset anchors the annotation next to the marker regardless of
                # the axis limits.
                callout = text!(
                    ax,
                    [Point2f(last(df.time), last(df.perigee_altitude))];
                    text     = callout_text,
                    align    = (:right, :bottom),
                    color    = accent_color,
                    font     = :bold,
                    fontsize = 15 * fontscale,
                    offset   = (-12, 12),
                )
                translate!(callout, 0, 0, 11)
            end
        end

        # == Axis Limits ===================================================================

        !isnothing(xlims) && xlims!(ax, xlims...)
        !isnothing(ylims) && ylims!(ax, ylims...)

        if show_f107
            push!(legend_plots, f107_line)
            push!(legend_labels, "F10.7")

            linkxaxes!(ax, ax_f107)

            # Align the twin axis ticks with the main axis grid using canonical values,
            # keeping them aligned if the main axis limits or ticks change.
            f107_min, f107_max = extrema(df.f107)

            onany(ax.finallimits, ax.yaxis.tickvalues) do main_limits, main_tickvalues
                _align_twin_yticks!(ax_f107, main_limits, main_tickvalues, f107_min, f107_max)
            end

            _align_twin_yticks!(
                ax_f107,
                ax.finallimits[],
                ax.yaxis.tickvalues[],
                f107_min,
                f107_max
            )
        end

        # == Side Column: Information Panel and Legend =====================================

        # The information panel and the legend share a single column layout, so they can
        # never overlap: the cards stack from the top, the legend card sits at the bottom,
        # and a flexible spacer row absorbs the remaining vertical space.
        side = GridLayout(fig[1, 2])

        for (k, (card_title, value_lines, compact)) in enumerate(cards)
            _add_stat_card!(
                side,
                k,
                card_title,
                value_lines;
                border_color  = border_color,
                card_color    = card_color,
                compact       = compact,
                fontscale     = fontscale,
                subline_color = subline_color,
                title_color   = title_color,
            )
        end

        # The legend is placed at the bottom of the side column, inside a card matching
        # the information panel style.
        legend_row = length(cards) + 2

        Box(
            side[legend_row, 1];
            color        = card_color,
            cornerradius = 8,
            strokecolor  = border_color,
            strokewidth  = 1,
        )

        legend_content = GridLayout(side[legend_row, 1]; alignmode = Outside(12))

        Label(
            legend_content[1, 1],
            "LEGEND";
            color     = title_color,
            font      = :bold,
            fontsize  = 14 * fontscale,
            halign    = :left,
            tellwidth = false,
        )

        Legend(
            legend_content[2, 1],
            legend_plots,
            legend_labels;
            backgroundcolor = :transparent,
            framevisible    = false,
            halign          = :left,
            labelsize       = 15 * fontscale,
            padding         = (0.0f0, 0.0f0, 0.0f0, 0.0f0),
            rowgap          = 2,
            tellheight      = true,
            tellwidth       = false,
        )

        rowgap!(legend_content, 6)

        # The spacer row between the cards and the legend absorbs the leftover space.
        rowsize!(side, legend_row - 1, Auto(false))
        rowgap!(side, 10)

        colsize!(fig.layout, 2, Fixed(panel_width′))

        # == Mission Name ==================================================================

        if !isnothing(mission_name)
            Label(
                fig[0, 1],
                uppercase(mission_name);
                color     = title_color,
                font      = :bold,
                fontsize  = 16 * fontscale,
                tellwidth = false,
            )

            # Tighten the gap between the mission name and the plot title.
            rowgap!(fig.layout, 1, 4)
        end

        return fig, ax
    end
end

############################################################################################
#                                    Private Functions                                     #
############################################################################################

"""
    _add_stat_card!(
        panel::GridLayout,
        row::Int,
        title::String,
        value_lines::Vector{String};
        kwargs...
    ) -> Nothing

Add to `panel` at `row` a card with the statistic `title` and its `value_lines`, modifying
the layout of the figure that owns `panel`. The first element of `value_lines` is rendered
as the card value, whereas the other elements are rendered as smaller complementary lines.

# Keywords

- `border_color::Colorant`: Color of the card border.
- `card_color::Colorant`: Color of the card background.
- `compact::Bool`: If `true`, all the value lines are rendered with a small font size,
    yielding a denser card suited for multi-line textual content.
    (**Default**: `false`)
- `fontscale::Real`: Scale applied to the font sizes of the card.
    (**Default**: 1)
- `subline_color::Colorant`: Color of the card complementary lines.
- `title_color::Colorant`: Color of the card title.
"""
function _add_stat_card!(
    panel::GridLayout,
    row::Int,
    title::String,
    value_lines::Vector{String};
    border_color::Colorant,
    card_color::Colorant,
    compact::Bool = false,
    fontscale::Real = 1,
    subline_color::Colorant,
    title_color::Colorant,
)
    Box(
        panel[row, 1];
        color        = card_color,
        cornerradius = 8,
        strokecolor  = border_color,
        strokewidth  = 1,
    )

    card = GridLayout(panel[row, 1]; alignmode = Outside(12))

    Label(
        card[1, 1],
        title;
        color     = title_color,
        font      = :bold,
        fontsize  = 14 * fontscale,
        halign    = :left,
        tellwidth = false,
    )

    Label(
        card[2, 1],
        first(value_lines);
        font      = compact ? :regular : :bold,
        fontsize  = (compact ? 15 : 21) * fontscale,
        halign    = :left,
        tellwidth = false,
    )

    for k in 2:length(value_lines)
        Label(
            card[k + 1, 1],
            value_lines[k];
            color     = subline_color,
            fontsize  = 15 * fontscale,
            halign    = :left,
            tellwidth = false,
        )
    end

    rowgap!(card, 4)

    return nothing
end

"""
    _align_twin_yticks!(
        ax_twin::Axis,
        main_limits,
        main_tickvalues::Vector,
        data_min::Number,
        data_max::Number
    ) -> Nothing

Align the y-axis ticks of the twin axis `ax_twin` with the y-axis grid of the main axis
using canonical tick values, modifying the y-axis limits and the attribute `yticks` of
`ax_twin`. The main axis state is provided by its final limits `main_limits` and its current
tick values `main_tickvalues`, whereas `data_min` and `data_max` are the extrema of the data
plotted in the twin axis.

The algorithm selects the smallest canonical tick step (1, 2, 2.5, or 5 times a power of
ten) and the twin axis limits such that the data fits within the limits and each twin axis
tick, a multiple of the tick step, shares the screen position of a main axis tick.
"""
function _align_twin_yticks!(
    ax_twin::Axis,
    main_limits,
    main_tickvalues::Vector,
    data_min::Number,
    data_max::Number,
)
    ylo = minimum(main_limits)[2]
    yhi = maximum(main_limits)[2]
    Δ   = yhi - ylo

    ((Δ > 0) && (length(main_tickvalues) >= 2)) || return nothing

    t₁     = first(main_tickvalues)
    s_main = main_tickvalues[2] - t₁

    (s_main > 0) || return nothing

    # Data span, falling back to a small value if the data is constant.
    D = data_max - data_min
    (D <= 0) && (D = max(0.1 * abs(data_max), 1.0))

    # Number of main axis tick steps that span the full main axis limits.
    m = Δ / s_main

    # Find the smallest canonical tick step such that the data fits within the resulting
    # twin axis limits.
    s_twin = _canonical_step(D / m)
    u₁ = lo = hi = 0.0

    while true
        # `u₁` is the twin axis tick at the first main axis tick: the largest multiple of
        # the tick step that keeps the lower limit below the data minimum.
        u₁ = s_twin * floor((data_min + (t₁ - ylo) * s_twin / s_main) / s_twin)
        lo = u₁ - (t₁ - ylo) * s_twin / s_main
        hi = lo + m * s_twin

        (hi >= data_max) && break

        s_twin = _canonical_step(1.01 * s_twin)
    end

    ylims!(ax_twin, lo, hi)

    tickvalues = [u₁ + k * s_twin for k in 0:(length(main_tickvalues) - 1)]
    ax_twin.yticks = (tickvalues, _format_number.(tickvalues))

    return nothing
end

"""
    _canonical_step(x::Number) -> Float64

Return the smallest canonical step, defined as 1, 2, 2.5, or 5 times a power of ten, that
is equal to or greater than `x`. The input must be positive.
"""
function _canonical_step(x::Number)
    base = exp10(floor(log10(x)))

    for mult in (1.0, 2.0, 2.5, 5.0, 10.0)
        (mult * base >= x) && return mult * base
    end

    return 10.0 * base
end

"""
    _format_duration(Δ::Dates.Period) -> String

Format the duration `Δ` using hours, days, or Julian years, selecting the unit that best
fits its magnitude.
"""
function _format_duration(Δ::Dates.Period)
    days = Dates.value(Dates.Millisecond(Δ)) / 86_400_000

    if days < 1
        return string(round(24 * days; digits = 1)) * " h"
    elseif days < 365.25
        return string(round(days; digits = 1)) * " days"
    else
        return string(round(days / 365.25; digits = 2)) * " years"
    end
end

"""
    _format_number(v::Number) -> String

Format the number `v` with four significant digits, omitting the decimal part if the
rounded value is an integer and grouping the digits of large integers with thin spaces.
"""
function _format_number(v::Number)
    rounded = round(v; sigdigits = 4)
    str     = isinteger(rounded) ? string(Int(rounded)) : string(rounded)

    return _group_digits(str)
end

"""
    _group_digits(str::AbstractString) -> String

Insert a thin space between each group of three digits in the number in `str`, counting
from the right, when it is an integer with more than four digits.
"""
function _group_digits(str::AbstractString)
    neg  = startswith(str, '-')
    body = neg ? str[2:end] : str

    # Grouping only applies to plain integers with more than four digits.
    (occursin('.', body) || occursin('e', body) || length(body) <= 4) && return String(str)

    n       = length(body)
    grouped = sprint() do io
        for (k, c) in pairs(body)
            print(io, c)
            r = n - k
            (r > 0) && (r % 3 == 0) && print(io, ' ')
        end
    end

    return (neg ? "-" : "") * grouped
end
