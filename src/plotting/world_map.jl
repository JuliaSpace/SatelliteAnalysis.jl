## Description #############################################################################
#
# Function to create a plot with the world map.
#
############################################################################################

export plot_world_map

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

    This function **only works** after loading the package **GeoMakie.jl**. Furthermore, the
    user must also load one Makie.jl backend (CairoMakie.jl or GLMakie.jl, for example) to
    see the result.

# Keywords

- `theme::Symbol`: Theme variant used to style the figure, applied locally through the
    function `SatelliteAnalysis.makie_theme`. It can be `:light` or `:dark`.
    (**Default**: `:light`)
- `size::Tuple`: Size of the figure.
    (**Default**: `(1450, 800)`)

# Extended help

## Throws

- `ArgumentError`: If the theme variant in `theme` is not `:dark` or `:light`.
"""
function plot_world_map(args...; kwargs...)
    return error("Wrong input or the package GeoMakie.jl is not loaded.")
end
