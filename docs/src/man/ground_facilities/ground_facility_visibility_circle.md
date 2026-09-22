# Ground Facility Visibility Circle

```@meta
CurrentModule = SatelliteAnalysis
```

```@setup ground_facility_visibility_circle
using SatelliteAnalysis
```

We can compute the visibility circle of a ground facility using the function
`ground_facility_visibility_circle`:

```@docs; canonical = false
ground_facility_visibility_circle
```

!!! note

    If we want to verify if a satellite has line-of-sight to a ground facility, see the
    function [`is_ground_facility_visible`](@ref).

## Examples

We can obtain the visibility circle between the Amazonia-1 satellite and INPE's ground
station at Cuiabá, MT, Brazil, using:

```@setup ground_facility_visibility_circle_example
using SatelliteAnalysis, GeoJSON, CairoMakie

gf = ground_facility_visibility_circle(
    (-(15 + 33 / 60) |> deg2rad, -(56 + 04 / 60) |> deg2rad, 0),
    7130.982e3
)

countries_filename = fetch_country_polygons(; force_download = false)

country_polys = GeoJSON.read(countries_filename)

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
    country_polys.geometry;
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

If the user loads the package [GeoJSON.jl](https://github.com/JuliaGeo/GeoJSON.jl)
together with a [Makie.jl](https://docs.makie.org/stable/) backend, an extension is loaded
and adds the possibility to plot the ground facility visibility circle. In this case, the
following functions are available:

```@docs; canonical = false
plot_ground_facility_visibility_circles
```

```@docs; canonical = false
plot_ground_facility_visibility_circles!
```

### Example

The code:

```@repl ground_facility_visibility_circle
using GeoJSON, CairoMakie

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
