# Ground Repeating Orbits

```@meta
CurrentModule = SatelliteAnalysis
```

```@setup ground_repeating
using SatelliteAnalysis
```

A ground repeating orbit is any orbit in which the number of revolutions per day is a
rational number:

```math
R_d = I + \frac{N}{D}\ ,
```

where ``I, N, D \in \mathbb{N}``. Hence, this type of orbit repeats its ground track after
``D`` (solar) days, which is called the orbit cycle. Those orbits are largely used in remote
sensing missions because the satellite revisits the same regions periodically.

Since the ground track repeats, the distance between two adjacent ground tracks at the
Equator is constant. This information is paramount to design the payload of a remote sensing
mission: if the swath is smaller than this distance, there will be gaps in the coverage.

## Adjacent Track Distance

We can compute the distance between two adjacent ground tracks at the Equator using the
function `ground_repeating_orbit_adjacent_track_distance`:

```@docs; canonical = false
ground_repeating_orbit_adjacent_track_distance
```

## Adjacent Track Angle

The angle between two adjacent ground tracks at the Equator measured from the satellite
position is also useful to design the payload because it defines the required field of view
of a camera, for example. We can compute it using the function
`ground_repeating_orbit_adjacent_track_angle`:

```@docs; canonical = false
ground_repeating_orbit_adjacent_track_angle
```

## Examples

The Sun-synchronous orbit of the Amazonia-1 mission repeats its ground track after 5 days
with 14 + 2/5 revolutions per day. The distance [km] between two adjacent ground tracks at
the Equator is:

```@repl ground_repeating
a = 7130.982e3
e = 0.001111
i = 98.410 |> deg2rad

ground_repeating_orbit_adjacent_track_distance(a, e, i, 5) / 1000
```

Hence, the swath of the camera must be larger than this value to avoid gaps in the coverage.
The angle [°] between two adjacent ground tracks measured from the satellite is:

```@repl ground_repeating
ground_repeating_orbit_adjacent_track_angle(a, e, i, 5) |> rad2deg
```

!!! note

    The function [`design_sun_sync_ground_repeating_orbit`](@ref) lists all the
    Sun-synchronous, ground-repeating orbits in a range of orbit cycles together with those
    two quantities. See [Sun-Synchronous Orbits](@ref) for more information.
