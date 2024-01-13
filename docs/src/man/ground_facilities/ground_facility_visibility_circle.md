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


