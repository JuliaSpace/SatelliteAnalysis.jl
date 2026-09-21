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
the objects `Figure` and `Axis` used to plot the data. For more information, please refer
to **Makie.jl** documentation.

!!! note

    This function plots the countries' borders in the created figure using the file with the
    country polygons fetched with the function [`fetch_country_polygons`](@ref). Hence, if
    this file does not exist, the algorithm tries to download it.

!!! warning

    This function **only works** after loading the package **GeoJSON.jl** and one
    **Makie.jl** backend (**CairoMakie.jl** or **GLMakie.jl**, for example).

# Keywords

- `ground_facilities::Union{Nothing, AbstractVector{<:Tuple}}`: Vector with the WGS84
    position of each ground facility `(latitude [rad], longitude [rad], altitude [m])`, as
    used to compute the visibility circles, which selects the position of the ground
    facility markers. If it is `nothing`, the positions are estimated using the visibility
    circles.
    (**Default** = `nothing`)
- `ground_facility_names::Union{Nothing, Vector{String}}`: The user can provide a vector of
    `String`s with the length of `vgf_vc` to be plotted with the visibility circles. If this
    parameter is `nothing`, no ground facility name is added to the figure.
    (**Default** = `nothing`)
- `theme::Union{Nothing, Symbol, Makie.Theme}`: Theme used to style the figure, which is
    applied locally. If it is a `Symbol`, it selects the variant of the theme created by the
    function `SatelliteAnalysis.makie_theme`, which can be `:light` or `:dark`. If it is a
    `Makie.Theme`, this theme is applied. If it is `nothing`, no theme is applied, and the
    figure uses the current Makie theme. In the last two cases, the elements that are not
    styled by the theme use the colors of the variant `:light`.
    (**Default**: `:light`)

All other `kwargs...` are passed to the function [`plot_world_map`](@ref).

# Extended help

## Throws

- `ArgumentError`: If the theme variant in `theme` is not `:dark` or `:light`.

## Examples

```julia
julia> using SatelliteAnalysis, GeoJSON, GLMakie

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
    return _extension_error(
        "plot_ground_facility_visibility_circles",
        "Makie.jl and GeoJSON.jl",
        "CairoMakie, GeoJSON",
        "plot_ground_facility_visibility_circles(vgf_vc::Vector{Vector{NTuple{2, T}}}; kwargs...) where {T <: Number}",
    )
end

"""
    plot_ground_facility_visibility_circles!(ax::Axis, vgf_vc::Vector{Vector{NTuple{2, T}}}; kwargs...) where {T <: Number} -> Vector{Lines}

Plot in the **Makie.jl** axis `ax` the ground facility visibility circles in the vector
`vgf_vc`, where each element is computed using the function
[`ground_facility_visibility_circle`](@ref). It returns a vector with the plot of each
visibility circle, which can be used, for example, to build a legend.

!!! note

    Since this function draws into an existing axis, it does not apply the theme provided by
    the function `SatelliteAnalysis.makie_theme`. The plot inherits the styling of the
    figure that owns `ax`.

!!! warning

    This function **only works** after loading the package **GeoJSON.jl** and one
    **Makie.jl** backend (**CairoMakie.jl** or **GLMakie.jl**, for example).

# Keywords

- `ground_facilities::Union{Nothing, AbstractVector{<:Tuple}}`: Vector with the WGS84
    position of each ground facility `(latitude [rad], longitude [rad], altitude [m])`, as
    used to compute the visibility circles, which selects the position of the ground
    facility markers. If it is `nothing`, the positions are estimated using the visibility
    circles.
    (**Default** = `nothing`)
- `ground_facility_names::Union{Nothing, Vector{String}}`: The user can provide a vector of
    `String`s with the length of `vgf_vc` to be plotted with the visibility circles. If this
    parameter is `nothing`, no ground facility name is added to the figure.
    (**Default** = `nothing`)

All other `kwargs...` are passed to the function `lines!` that plots each visibility circle,
allowing the selection of attributes such as `linestyle` and `linewidth`.

# Extended Help

## Examples

```julia
julia> using SatelliteAnalysis, GeoJSON, GLMakie

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
function plot_ground_facility_visibility_circles!(args...; kwargs...)
    return _extension_error(
        "plot_ground_facility_visibility_circles!",
        "Makie.jl and GeoJSON.jl",
        "CairoMakie, GeoJSON",
        "plot_ground_facility_visibility_circles!(ax::Axis, vgf_vc::Vector{Vector{NTuple{2, T}}}; kwargs...) where {T <: Number}",
    )
end
