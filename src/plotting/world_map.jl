## Description #############################################################################
#
# Function to create a plot with the world map.
#
############################################################################################

export plot_world_map

"""
    plot_world_map(; kwargs...) -> Figure, Axis

Create a **Makie.jl** `Figure` and `Axis` with the world map. All the `kwargs...` are passed
to the function `Figure`. For more information, please, refer to **Makie.jl** documentation.

!!! note

    This function plots the countries' borders in the created figure using the file with the
    country polygons fetched with the function [`fetch_country_polygons`](@ref). Hence, if
    this files does not exist, the algorithm tries to download it.
"""
function plot_world_map(args...)
    error("Wrong input or the package GeoMakie.jl is not loaded.")
end
