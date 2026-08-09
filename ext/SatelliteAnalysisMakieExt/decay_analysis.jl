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
    satellite_mass::Union{Nothing, Number}      = nothing,
    satellite_mean_area::Union{Nothing, Number} = nothing,
    terminate_altitude::Union{Nothing, Number}  = nothing,
    theme::Symbol = :light,
    size = (1280, 720),
    kwargs...
)
    # == Input Validation ==================================================================

    for c in (:time, :date, :apogee_altitude, :perigee_altitude)
        hasproperty(df, c) || throw(
            ArgumentError(
                "The input `DataFrame` must have the column `$c`. It should be obtained " *
                "using the function `decay_analysis`."
            )
        )
    end

    isempty(df) && throw(ArgumentError("The input `DataFrame` is empty."))

    # Build the theme first since it also validates the variant in `theme`.
    sa_theme = makie_theme(theme)

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

    # == Units =============================================================================

    time_unit     = colmetadata(df, :time,            "Unit", :y)
    distance_unit = colmetadata(df, :apogee_altitude, "Unit", :km)

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
    secondary_color = dark ? TEXT_SECONDARY_DARK : TEXT_SECONDARY_LIGHT

    # == Assemble the Information Panel Cards =============================================

    cards = Tuple{String, Vector{String}}[]

    !isnothing(mass) && push!(cards, ("SATELLITE MASS", ["$(_format_number(mass)) kg"]))
    !isnothing(area) && push!(cards, ("MEAN AREA", ["$(_format_number(area)) m²"]))

    if !isnothing(term)
        if reentered
            reentry_epoch = Dates.format(last(df.date), dateformat"yyyy-mm-dd HH:MM")
            push!(
                cards,
                (
                    "TIME TO REENTER",
                    [
                        _format_duration(last(df.date) - first(df.date)),
                        reentry_epoch * " UTC"
                    ]
                )
            )
        else
            push!(cards, ("TIME TO REENTER", ["No reentry"]))
        end
    end

    # == Plot ==============================================================================

    # Every object must be created inside `with_theme` because Makie resolves the theme
    # attributes at object-creation time.
    return with_theme(sa_theme) do
        fig = Figure(; size = size, kwargs...)

        ax = Axis(
            fig[1, 1];
            title  = "Orbital Decay Analysis",
            xlabel = "Time [$time_unit]",
            ylabel = "Altitude [$distance_unit]",
        )

        lines!(ax, df.time, df.apogee_altitude;  label = "Apogee Altitude")
        lines!(ax, df.time, df.perigee_altitude; label = "Perigee Altitude")

        if reentered
            reentry_marker = scatter!(
                ax,
                [last(df.time)],
                [last(df.perigee_altitude)];
                color      = accent_color,
                label      = "Reentry",
                markersize = 14,
            )
            translate!(reentry_marker, 0, 0, 10)
        end

        axislegend(ax; position = :rt)

        if !isempty(cards)
            panel = GridLayout(fig[1, 2]; tellheight = false, valign = :top)
            colsize!(fig.layout, 2, Fixed(280))

            for (k, (title, value_lines)) in enumerate(cards)
                _add_stat_card!(
                    panel,
                    k,
                    title,
                    value_lines;
                    border_color    = border_color,
                    card_color      = card_color,
                    secondary_color = secondary_color,
                )
            end

            (length(cards) > 1) && rowgap!(panel, 15)
        end

        return fig, ax
    end
end

############################################################################################
#                                    Private Functions                                     #
############################################################################################

"""
    _add_stat_card!(panel::GridLayout, row::Int, title::String, value_lines::Vector{String}; kwargs...) -> Nothing

Add to `panel` at `row` a card with the statistic `title` and its `value_lines`, modifying
the layout of the figure that owns `panel`. The first element of `value_lines` is rendered
as the card value, whereas the other elements are rendered as smaller complementary lines.

# Keywords

- `border_color::Colorant`: Color of the card border.
- `card_color::Colorant`: Color of the card background.
- `secondary_color::Colorant`: Color of the card title and complementary lines.
"""
function _add_stat_card!(
    panel::GridLayout,
    row::Int,
    title::String,
    value_lines::Vector{String};
    border_color::Colorant,
    card_color::Colorant,
    secondary_color::Colorant,
)
    Box(
        panel[row, 1];
        color        = card_color,
        cornerradius = 8,
        strokecolor  = border_color,
        strokewidth  = 1,
    )

    card = GridLayout(panel[row, 1]; alignmode = Outside(16))

    Label(
        card[1, 1],
        title;
        color     = secondary_color,
        font      = :bold,
        fontsize  = 14,
        halign    = :left,
        tellwidth = false,
    )

    Label(
        card[2, 1],
        first(value_lines);
        font      = :bold,
        fontsize  = 24,
        halign    = :left,
        tellwidth = false,
    )

    for k in 2:length(value_lines)
        Label(
            card[k + 1, 1],
            value_lines[k];
            color     = secondary_color,
            fontsize  = 15,
            halign    = :left,
            tellwidth = false,
        )
    end

    rowgap!(card, 6)

    return nothing
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
rounded value is an integer.
"""
function _format_number(v::Number)
    rounded = round(v; sigdigits = 4)
    return isinteger(rounded) ? string(Int(rounded)) : string(rounded)
end
