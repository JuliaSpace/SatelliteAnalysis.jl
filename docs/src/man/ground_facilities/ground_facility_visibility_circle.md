Ground Facility Visibility Circle
=================================

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

```
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

```@repl ground_facility_visibility_circle
ground_facility_visibility_circle(
    (-(15 + 33 / 60) |> deg2rad, -(56 + 04 / 60) |> deg2rad, 0),
    7130.982e3
)
```

If we plot the result using Makie, we obtain:

```@raw html
<div align="center">
  <img src="../../../assets/ground_station_visibility_circle.png" alt="Ground Station Visibility Circle" width="100%"/>
</div>
```

## Plotting

If the user loads the package [GeoMakie.jl](https://github.com/MakieOrg/GeoMakie.jl)
together with a [Makie.jl](https://docs.makie.org/stable/) back end, an extension is loaded
and adds the possibility to plot the ground facility visibility circle. In this case, the
following function is available:

```julia
plot_ground_facility_visibility_circles(vgf_vc::Vector{Vector{NTuple{2, T}}}; kwargs...) where T <: Number -> Figure, Axis
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

All other `kwargs...` are passed to the function `Figure`.

### Example

The code:

```julia-repl ground_track_plotting
julia> using GeoMakie, GLMakie

julia> gf1_vc = ground_facility_visibility_circle(
           (-(15 + 33 / 60) |> deg2rad, -(56 + 04 / 60) |> deg2rad, 0),
           7130.982e3
       );

julia> gf2_vc = ground_facility_visibility_circle(
           (-22.6763 |> deg2rad, -44.9973 |> deg2rad, 0),
           7130.982e3
       );

julia> gf3_vc = ground_facility_visibility_circle(
           (+78.228 |> deg2rad, +15.399 |> deg2rad, 0),
           7130.982e3
       );

julia> fig, ax = plot_ground_facility_visibility_circles(
           [gf1_vc, gf2_vc, gf3_vc];
           ground_facility_names = ["Cuiabá", "Cachoeira Paulista", "Svalbard"]
       );

julia> fig
```

produces the following figure if **GLMakie.jl** is loaded:

```@raw html
<div align="center">
  <img src="../../../assets/ground_facility_plotting_extension.png" alt="Ground Facility Visibility Circle Plot" width="100%"/>
</div>
```
