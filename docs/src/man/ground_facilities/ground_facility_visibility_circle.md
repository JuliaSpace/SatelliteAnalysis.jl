# Ground Facility Visibility Circle

```@meta
CurrentModule = SatelliteAnalysis
```

```@setup ground_facility_visibility_circle
using SatelliteAnalysis
```

We can use the function:

```julia
ground_facility_visibility_circle(gf_wgs84::Tuple, satellite_position_norm::Number; kwargs...) -> Vector{NTuple{2, Float64}}
```

to compute the ground facility visibility circle from the position `gf_wgs84` (WGS84) to a
satellite in which its distance from the Earth's center is `satellite_position_norm` [m]. It
returns a vector of `NTuple{2, Float64}` where the first element is the latitude [rad] and
the second is the longitude [rad] of each point in the visibility circle.

The ground facility is specified using a tuple with its WGS84 position:

```text
(latitude [rad], longitude [rad], altitude [m])
```

The following keywords are available:

- `azimuth_step::Number`: The step in the azimuth used to compute the visibility circle.
    (**Default**: `0.1 |> deg2rad`)
- `minimum_elevation::Number`: Minimum elevation angle for communication between the
    satellite and the ground facility [rad].
    (**Default**: `10 |> deg2rad`)

!!! note

    If we want to verify if a satellite has line-of-sight to a ground facility, see the
    function [`is_ground_facility_visible`](@ref).

## Examples

We can obtain the visibility circle between the Amazonia-1 satellite and INPE's ground
station at Cuiabá, MT, Brazil, using:

```@setup ground_facility_visibility_circle_example
using SatelliteAnalysis, GeoMakie, CairoMakie

gf = ground_facility_visibility_circle(
    (-(15 + 33 / 60) |> deg2rad, -(56 + 04 / 60) |> deg2rad, 0),
    7130.982e3
)

countries_filename = fetch_country_polygons(; force_download = false)

country_polys = GeoMakie.GeoJSON.read(countries_filename)

fig = Figure(; size = (800, 800))

ax = Axis(
    fig[1, 1],
    aspect         = 1,
    title          = "Ground Facility Visibility",
    titlegap       = 16,
    titlesize      = 30,
    xlabel         = "Longitude [°]",
    xlabelsize     = 30,
    xticklabelsize = 26,
    ylabel         = "Latitude [°]",
    ylabelsize     = 30,
    yticklabelsize = 26,
)

xlims!(ax, -95, -15)
ylims!(ax, -60, +20)
ax.xticks = -95:10:-15
ax.yticks = -60:10:20

poly!(
    ax,
    GeoMakie.to_multipoly(country_polys.geometry);
    color       = :white,
    strokecolor = :black,
    strokewidth = 1
)

gf_lat = first.(gf)
gf_lon = last.(gf)

vc = lines!(ax, gf_lon .|> rad2deg, gf_lat .|> rad2deg; linewidth = 2)

dot = scatter!(ax, -(56 + 04 / 60), -(15 + 33 / 60))
translate!(dot, 0, 0, 10)

label = text!(
    ax,
    "Cuiabá";
    fontsize = 26,
    position = (-(56 + 04 / 60), -(15 + 33 / 60))
)
translate!(label, 0, 0, 10)

save("gf_visibility_circle_01.png", fig)
```

```@repl ground_facility_visibility_circle
ground_facility_visibility_circle(
    (-(15 + 33 / 60) |> deg2rad, -(56 + 04 / 60) |> deg2rad, 0),
    7130.982e3
)
```

If we plot the result using Makie, we obtain:

![Cuiabá ground facility visibility circle](./gf_visibility_circle_01.png)

## Plotting

If the user loads the package [GeoMakie.jl](https://github.com/MakieOrg/GeoMakie.jl)
together with a [Makie.jl](https://docs.makie.org/stable/) back end, an extension is loaded
and adds the possibility to plot the ground facility visibility circle. In this case, the
following functions are available:

```julia
plot_ground_facility_visibility_circles(vgf_vc::Vector{Vector{NTuple{2, Number}}}; kwargs...) -> Figure, Axis
```

It plots the ground facility visibility circles in the vector `vgf_vc`, where each element
is computed using the function [`ground_facility_visibility_circle`](@ref). It returns the
objects `Figure` and `Axis` used to plot the data. For more information, please, refer to
[Makie.jl](https://docs.makie.org/stable/) documentation.

The following keywords are available:

- `ground_facility_names::Union{Nothing, Vector{String}}`: The user can provide a vector of
    `String`s with the length of `vgf_vc` to be plotted with the visibility circles. If this
    parameter is `nothing`, no ground facility name is added to the figure.
    (**Default** = `nothing`)

All other `kwargs...` are passed to the function [`plot_world_map`](@ref).

!!! note

    This function plots the countries' borders in the created figure using the file with the
    country polygons fetched with the function [`fetch_country_polygons`](@ref). Hence, if
    this files does not exist, the algorithm tries to download it.

```julia
plot_ground_facility_visibility_circles!(ax::Axis, vgf_vc::Vector{Vector{NTuple{2, Number}}}; kwargs...) -> Nothing
```

It plots in the [Makie.jl](https://docs.makie.org/stable/) axis `ax` the ground facility
visibility circles in the vector `vgf_vc`, where each element is computed using the function
[`ground_facility_visibility_circle`](@ref).

The following keywords are available:

- `ground_facility_names::Union{Nothing, Vector{String}}`: The user can provide a vector of
    `String`s with the length of `vgf_vc` to be plotted with the visibility circles. If this
    parameter is `nothing`, no ground facility name is added to the figure.
    (**Default** = `nothing`)

The user can use this function to plot the ground facility visibility circle on top of an
existing figure.

### Example

The code:

```@repl ground_facility_visibility_circle
using GeoMakie, CairoMakie

gf1_vc = ground_facility_visibility_circle(
    (-(15 + 33 / 60) |> deg2rad, -(56 + 04 / 60) |> deg2rad, 0),
    7130.982e3
);

gf2_vc = ground_facility_visibility_circle(
    (-22.6763 |> deg2rad, -44.9973 |> deg2rad, 0),
    7130.982e3
);

gf3_vc = ground_facility_visibility_circle(
    (+78.228 |> deg2rad, +15.399 |> deg2rad, 0),
    7130.982e3
);

fig, ax = plot_ground_facility_visibility_circles(
    [gf1_vc, gf2_vc, gf3_vc];
    ground_facility_names = ["Cuiabá", "Cachoeira Paulista", "Svalbard"]
);

save("gf_visibility_circle_02.png", fig)
```

produces the following figure:

![Ground facility visibility circles](./gf_visibility_circle_02.png)
