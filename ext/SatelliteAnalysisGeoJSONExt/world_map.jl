## Description #############################################################################
#
# Function to create a plot with the world map.
#
############################################################################################

function SatelliteAnalysis.plot_world_map(;
    theme::Union{Nothing, Symbol, Makie.Theme} = :light,
    size = (1450, 800),
    kwargs...
)
    # Every object must be created inside this function because Makie resolves the theme
    # attributes at object-creation time.
    return SatelliteAnalysis._with_plot_theme(theme) do
        return _create_world_map(_theme_variant(theme); size = size, kwargs...)
    end
end

############################################################################################
#                                    Private Functions                                     #
############################################################################################

"""
    _theme_variant(theme::Union{Nothing, Symbol, Makie.Theme}) -> Symbol

Return the variant of the SatelliteAnalysis.jl theme used to select the colors of the
elements that are not styled by the Makie theme, such as the country polygons. If `theme` is
a `Symbol`, it is the variant. Otherwise, the variant is `:light`.
"""
_theme_variant(theme::Symbol) = theme
_theme_variant(::Union{Nothing, Makie.Theme}) = :light

"""
    _create_world_map(variant::Symbol; kwargs...) -> Figure, Axis

Create a figure with the world map using the country colors of the theme `variant` (`:light`
or `:dark`). This function does not apply any theme. Hence, the public functions must call
it inside `with_theme`, which allows them to build the theme only once.

# Keywords

- `size::Tuple`: Figure size.
    (**Default**: `(1450, 800)`)

All the other keywords are passed to the function `Figure`.
"""
function _create_world_map(variant::Symbol; size = (1450, 800), kwargs...)
    fig = Figure(; size = size, kwargs...)

    ax = Axis(
        fig[1, 1];
        aspect = 2,
        title  = "World Map",
        xlabel = "Longitude [°]",
        ylabel = "Latitude [°]",
    )

    xlims!(ax, -180, +180)
    ylims!(ax, -90, +90)
    ax.xticks = -180:20:+180
    ax.yticks = -90:15:90

    _draw_world_map!(ax, variant)

    return fig, ax
end

"""
    _draw_world_map!(ax::Axis, variant::Symbol) -> Poly

Draw the country polygons in the axis `ax` using the colors of the theme `variant` (`:light`
or `:dark`), returning the created plot. The polygons are obtained from the file fetched by
the function [`fetch_country_polygons`](@ref). Hence, if this file does not exist, the
algorithm tries to download it.
"""
function _draw_world_map!(ax::Axis, variant::Symbol)
    # Get the GeoJSON file with the countries.
    countries_filename = fetch_country_polygons(; force_download = false)

    # Load the polygons of the countries. Notice that **GeoJSON.jl** provides the recipes
    # to plot its geometries when **Makie.jl** is loaded.
    country_polys = GeoJSON.read(countries_filename)

    return poly!(
        ax,
        country_polys.geometry;
        color       = _COUNTRY_FILL_COLOR[variant],
        strokecolor = _COUNTRY_STROKE_COLOR[variant],
        strokewidth = 1,
    )
end
