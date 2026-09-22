## Description #############################################################################
#
# Function to create a plot with the world map.
#
############################################################################################

export plot_world_map, plot_world_map!

"""
    plot_world_map(; kwargs...) -> Figure, Axis

Create a **Makie.jl** `Figure` and `Axis` with the world map. The figure is styled with the
theme obtained from the function `SatelliteAnalysis.makie_theme`, selected by the keyword
`theme`. All other `kwargs...` are passed to the function `Figure`. For more information,
please refer to **Makie.jl** documentation.

!!! note

    This function plots the countries' borders in the created figure using the file with the
    country polygons fetched with the function [`fetch_country_polygons`](@ref). Hence, if
    this file does not exist, the algorithm tries to download it.

!!! warning

    This function **only works** after loading the package **GeoJSON.jl** and one
    **Makie.jl** backend (**CairoMakie.jl** or **GLMakie.jl**, for example).

# Keywords

- `theme::Union{Nothing, Symbol, Makie.Theme}`: Theme used to style the figure, which is
    applied locally. If it is a `Symbol`, it selects the variant of the theme created by the
    function `SatelliteAnalysis.makie_theme`, which can be `:light` or `:dark`. If it is a
    `Makie.Theme`, this theme is applied. If it is `nothing`, no theme is applied, and the
    figure uses the current Makie theme. In the last two cases, the elements that are not
    styled by the theme use the colors of the variant `:light`.
    (**Default**: `:light`)
- `size::Tuple`: Size of the figure.
    (**Default**: `(1450, 800)`)

# Extended help

## Throws

- `ArgumentError`: If the theme variant in `theme` is not `:dark` or `:light`.
"""
function plot_world_map(args...; kwargs...)
    return _extension_error(
        "plot_world_map",
        "Makie.jl and GeoJSON.jl",
        "CairoMakie, GeoJSON",
        "plot_world_map(; kwargs...)",
    )
end

"""
    plot_world_map!(ax::Axis; kwargs...) -> Poly

Plot the country polygons of the world map in the **Makie.jl** axis `ax`, in which the
X-axis is the longitude [°] and the Y-axis is the latitude [°], returning the created plot.
This function does not change the axis limits, ticks, or labels.

!!! note

    This function plots the countries' borders using the file with the country polygons
    fetched with the function [`fetch_country_polygons`](@ref). Hence, if this file does not
    exist, the algorithm tries to download it.

!!! note

    Since this function draws into an existing axis, it does not apply the theme provided by
    the function `SatelliteAnalysis.makie_theme`. The plot inherits the styling of the
    figure that owns `ax`.

!!! warning

    This function **only works** after loading the package **GeoJSON.jl** and one
    **Makie.jl** backend (**CairoMakie.jl** or **GLMakie.jl**, for example).

# Keywords

- `theme::Union{Nothing, Symbol, Makie.Theme}`: Select the colors of the country polygons,
    which are not styled by the Makie theme. If it is `:dark`, we use the colors of the dark
    variant of the theme provided by `SatelliteAnalysis.makie_theme`. Otherwise, we use the
    colors of the light variant.
    (**Default**: `:light`)

All other `kwargs...` are passed to the function `poly!`, overriding the attributes selected
by this function (`color`, `strokecolor`, and `strokewidth`).
"""
function plot_world_map!(args...; kwargs...)
    return _extension_error(
        "plot_world_map!",
        "Makie.jl and GeoJSON.jl",
        "CairoMakie, GeoJSON",
        "plot_world_map!(ax::Axis; kwargs...)",
    )
end
