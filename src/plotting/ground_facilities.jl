## Description #############################################################################
#
# Functions to plot results related to the ground facilities.
#
############################################################################################

export plot_ground_facility_visibility_circles, plot_ground_facility_visibility_circles!

"""
    plot_ground_facility_visibility_circles(vgf_vc::Vector{Vector{NTuple{2, Number}}}; kwargs...) -> Figure, Axis

Plot the ground facility visibility circles in the vector `vgf_vc`, where each element
is computed using the function [`ground_facility_visibility_circle`](@ref). It returns
the objects `Figure` and `Axis` used to plot the data. For more information, please, refer
to **Makie.jl** documentation.

!!! note

    This function plots the countries' borders in the created figure using the file with the
    country polygons fetched with the function [`fetch_country_polygons`](@ref). Hence, if
    this files does not exist, the algorithm tries to download it.

!!! warning

    This function **only works** after loading the package **GeoMakie.jl**. Furthermore, the
    user must also load one Makie.jl back end (CairoMakie.jl or GLMakie.jl, for example) to
    see the result.

# Keywords

- `ground_facility_names::Union{Nothing, Vector{String}}`: The user can provide a vector of
    `String`s with the length of `vgf_vc` to be plotted with the visibility circles. If this
    parameter is `nothing`, no ground facility name is added to the figure.
    (**Default** = `nothing`)

All other `kwargs...` are passed to the function [`plot_world_map`](@ref).

# Extended Help

## Examples

```julia
julia> using SatelliteAnalysis, GeoMakie, GLMakie

julia> gfv1 = ground_facility_visibility_circle((0, 0, 0), EARTH_EQUATORIAL_RADIUS + 700e3);

julia> gfv2 = ground_facility_visibility_circle((-40 |> deg2rad, -60 |> deg2rad, 0), EARTH_EQUATORIAL_RADIUS + 700e3);

julia> fig, ax = plot_ground_facility_visibility_circles(
           [gfv1, gfv2];
           ground_facility_names = ["GF 1", "GF 2"]
       )
(Scene (1600px, 800px):
  0 Plots
  1 Child Scene:
    └ Scene (1600px, 800px), Axis (7 plots))

julia> fig
```
"""
function plot_ground_facility_visibility_circles(args...; kwargs...)
    error("Wrong input or the package GeoMakie.jl is not loaded.")
    return nothing
end

"""
    plot_ground_facility_visibility_circles!(ax::Axis, vgf_vc::Vector{Vector{NTuple{2, Number}}}; kwargs...) -> Nothing

Plot in the **Makie.jl** axis `ax` the ground facility visibility circles in the vector
`vgf_vc`, where each element is computed using the function
[`ground_facility_visibility_circle`](@ref).

!!! warning

    This function **only works** after loading the package **GeoMakie.jl**. Furthermore, the
    user must also load one Makie.jl back end (CairoMakie.jl or GLMakie.jl, for example) to
    see the result.

# Keywords

- `ground_facility_names::Union{Nothing, Vector{String}}`: The user can provide a vector of
    `String`s with the length of `vgf_vc` to be plotted with the visibility circles. If this
    parameter is `nothing`, no ground facility name is added to the figure.
    (**Default** = `nothing`)

# Extended Help

## Examples

```julia
julia> using SatelliteAnalysis, GeoMakie, GLMakie

julia> gfv1 = ground_facility_visibility_circle((0, 0, 0), EARTH_EQUATORIAL_RADIUS + 700e3);

julia> gfv2 = ground_facility_visibility_circle((-40 |> deg2rad, -60 |> deg2rad, 0), EARTH_EQUATORIAL_RADIUS + 700e3);

julia> fig = Figure(size = (1000, 1000))

julia> ax = Axis(fig[1, 1])

julia> plot_ground_facility_visibility_circles!(
           ax,
           [gfv1, gfv2];
           ground_facility_names = ["GF 1", "GF 2"]
       )
```
"""
function plot_ground_facility_visibility_circles!(args...)
    error("Wrong input or the package GeoMakie.jl is not loaded.")
    return nothing
end
