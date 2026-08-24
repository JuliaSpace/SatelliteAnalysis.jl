## Description #############################################################################
#
# Function to plot the decay analysis.
#
############################################################################################

export plot_decay_analysis

"""
    plot_decay_analysis(df::DataFrame; kwargs...) -> Figure, Axis

Plot the decay analysis in `df`, computed using the function [`decay_analysis`](@ref). It
returns the objects `Figure` and `Axis` used to plot the data. For more information, please
refer to **Makie.jl** documentation.

The figure shows the evolution of the mean apogee and perigee altitudes and an information
panel with the satellite mass, the satellite mean area, the estimated time to reenter, and
a card with the analysis assumptions (atmospheric model, F10.7 source, and drag and SRP
coefficients).
A dashed line marks the terminate altitude used to declare the reentry, and the estimated
reentry date is annotated next to the reentry marker.
The legend is placed outside the plot, at the bottom of the column that contains the
information panel, inside a card that matches the information panel style.
The time to reenter is only shown if the analysis detected a reentry, _i.e._, if the mean
perigee altitude reached the terminate altitude. The satellite mass [kg], the satellite
mean area [m²], and the terminate altitude [m] are obtained from the `DataFrame` metadata
written by [`decay_analysis`](@ref) (`Satellite Mass`, `Satellite Mean Area`, and
`Terminate Altitude`), but they can be overridden using keywords. The assumptions are
resolved from the metadata `Atmospheric Model`, `Drag Coefficient`, `SRP Coefficient`, and
`F10.7 Source`. Information whose value cannot be resolved, because the metadata is absent
and the related keyword was not passed, is omitted from the panel.

!!! warning

    This function **only works** after loading the package **Makie.jl**. Furthermore, the
    user must also load one Makie.jl backend (CairoMakie.jl or GLMakie.jl, for example) to
    see the result.

# Keywords

- `fontscale::Real`: Factor to uniformly scale every font size of the figure, useful when
    rendering at a size other than the default.
    (**Default**: 1)
- `mission_name::Union{Nothing, String}`: Mission name rendered in uppercase above the plot
    title. If it is `nothing`, no mission name is added to the figure.
    (**Default**: `nothing`)
- `mono_ticklabels::Bool`: If `true`, the tick labels are rendered using a monospaced
    font.
    (**Default**: `false`)
- `panel_width::Union{Nothing, Int}`: Width [px] of the column with the information panel
    and the legend. If it is `nothing`, the width scales with the figure width.
    (**Default**: `nothing`)
- `satellite_mass::Union{Nothing, Number}`: Satellite mass [kg] shown in the information
    panel. If it is `nothing`, the value is obtained from the metadata `Satellite Mass` of
    `df`.
    (**Default**: `nothing`)
- `satellite_mean_area::Union{Nothing, Number}`: Satellite mean area [m²] shown in the
    information panel. If it is `nothing`, the value is obtained from the metadata
    `Satellite Mean Area` of `df`.
    (**Default**: `nothing`)
- `show_assumptions::Bool`: If `true`, the information panel shows a card with the
    analysis assumptions resolved from the `DataFrame` metadata written by
    [`decay_analysis`](@ref).
    (**Default**: `true`)
- `show_f107::Bool`: If `true`, the 10.7 cm solar flux index in the column `f107` of `df`
    is plotted [sfu] using a twin y-axis placed at the right side of the figure. The line
    is rendered with transparency below the other plot elements, and the twin y-axis ticks
    use canonical values aligned with the grid of the main axis.
    (**Default**: `false`)
- `show_reentry_callout::Bool`: If `true` and the analysis detected a reentry, the
    estimated reentry date [UTC] is annotated in the plot next to the reentry marker.
    (**Default**: `true`)
- `show_reentry_date::Bool`: If `true`, the estimated reentry date [UTC] is shown in the
    information panel below the time to reenter. Otherwise, only the timespan is shown.
    (**Default**: `false`)
- `subtitle::Union{Nothing, String, Symbol}`: Subtitle rendered below the plot title. If it
    is `:auto`, the subtitle shows the analysis timespan [UTC]. If it is `nothing`, no
    subtitle is added to the figure. Any other `Symbol` raises an `ArgumentError`.
    (**Default**: `:auto`)
- `terminate_altitude::Union{Nothing, Number}`: Mean perigee altitude [m] that terminates
    the decay analysis, used to detect if a reentry happened. If it is `nothing`, the value
    is obtained from the metadata `Terminate Altitude` of `df`.
    (**Default**: `nothing`)
- `theme::Symbol`: Theme variant used to style the figure, applied locally through the
    function `SatelliteAnalysis.makie_theme`. It can be `:light` or `:dark`.
    (**Default**: `:light`)
- `title::String`: Title of the plot.
    (**Default**: `"Orbital Decay Analysis"`)
- `xlims::Union{Nothing, Tuple}`: Limits of the x-axis of the main plot. If it is
    `nothing`, the limits are computed automatically.
    (**Default**: `nothing`)
- `ylims::Union{Nothing, Tuple}`: Limits of the y-axis of the main plot. If it is
    `nothing`, the limits are computed automatically.
    (**Default**: `nothing`)
- `size::Tuple`: Size of the figure.
    (**Default**: `(1280, 720)`)

All other `kwargs...` are passed to the function `Figure`.

To export the figure in high resolution for reports, use
`save("plot.png", fig; px_per_unit = 2)`.

# Extended help

## Throws

- `ArgumentError`: If `df` does not have the columns `time`, `date`, `apogee_altitude`, and
    `perigee_altitude`, or if the keyword `show_f107` is `true` and `df` does not have the
    column `f107`.
- `ArgumentError`: If the theme variant in `theme` is not `:dark` or `:light`.
- `ArgumentError`: If the keyword `subtitle` is a `Symbol` other than `:auto`.

## Examples

```julia-repl
julia> using SatelliteAnalysis, OrdinaryDiffEqAdamsBashforthMoulton, CairoMakie

julia> jd₀ = date_to_jd(2024, 1, 1);

julia> orb = KeplerianElements(
           jd₀,
           EARTH_EQUATORIAL_RADIUS + 300e3,
           0.001,
           98.0 |> deg2rad,
           ltdn_to_raan(10.5, jd₀),
           90 |> deg2rad,
           0
       );

julia> df = decay_analysis(
           orb;
           satellite_mass = 100.0,
           satellite_mean_area = 1.0,
           F107 = 140
       );

julia> fig, ax = plot_decay_analysis(df);

julia> fig
```
"""
function plot_decay_analysis(::Any; kwargs...)
    return error("Wrong input or the package Makie.jl is not loaded.")
end
