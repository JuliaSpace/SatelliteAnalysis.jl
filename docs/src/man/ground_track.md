# Ground Track

```@meta
CurrentModule = SatelliteAnalysis
```

```@setup ground_track
using SatelliteAnalysis
```

We can obtain the ground track of a satellite using the function `ground_track`:

```@docs; canonical = false
ground_track
```

## Ground Track Inclination

Many analyses require the ground track inclination. For example, if we are designing a
remote sensing mission with an optical payload, we must know how much two images overlap.
This information can only be computed using spherical trigonometry and the ground track
inclination instead of the orbital one.

The ground track inclination ``i_{gt}``, shown in the following figure, is a composition of
the orbit inclination, the Earth's angular speed, and the RAAN time derivative. We can
compute it using:

```math
i_{gt} = \tan^{-1}\left(\frac{
    \omega_s \sin i
}{
    \omega_s \cos i - \omega_e + \dot{\Omega}
}\right)\ ,
```

where ``i`` is the orbital inclination, ``\omega_s`` is the satellite angular speed at
Equator, ``\omega_e`` is the Earth angular speed, and ``\dot{\Omega}`` is the RAAN time derivative.

```@raw html
<div align="center">
  <img src="../../assets/ground_track_inclination.png" alt="Ground Track Inclination" width="100%"/>
</div>
```

Formally, we should use the satellite instantaneous angular speed at the Equator in
``\omega_s``. However, given the perturbations caused by the Earth's gravitational
potential, this speed is not simple to compute. The calculation would require implementing
an orbit propagator. Thus, we simplify it by assuming that the orbit eccentricity is small.
This assumption is reasonable given the missions that would benefit from the computation of
the ground track inclination. In this case, we approximate ``\omega_s`` as the mean
satellite angular speed.

We can compute it using the function `ground_track_inclination`:

```@docs; canonical = false
ground_track_inclination
```

## Examples

We will compute the descending ground tracks of the Amazonia-1 mission for five days. The
first thing we need to do is define the orbit:

```@repl ground_track
jd₀ = date_to_jd(2021, 1, 1)

orb = KeplerianElements(
    jd₀,
    7130.982e3,
    0.001111,
    98.405 |> deg2rad,
    ltdn_to_raan(10.5, jd₀),
    π / 2,
    0
)
```

The next step is to define the desired propagator:

```@repl ground_track
orbp = Propagators.init(Val(:J2), orb)
```

Now, we can use the function `ground_track` to obtain the satellite ground track considering
only the descending passages:

```@repl ground_track
gt = ground_track(orbp; duration = 5 * 86400, track_types = :descending)
```

Finally, we can extract the latitude and longitude of each point in the ground track using:

```@repl ground_track
gt_lat = first.(gt)
gt_lon = last.(gt)
```

If we use **Makie.jl** to plot, we obtain:

```@setup ground_track
using GeoJSON, CairoMakie

fig, ax = plot_ground_track(gt)

ax.title = "Amazonia-1 Descending Ground Tracks"

save("amz1_descending_ground_tracks.png", fig)
```

![Amazonia-1 descending ground track](./amz1_descending_ground_tracks.png)

Finally, the ground track inclination is:

```@repl ground_track
ground_track_inclination(orb) |> rad2deg
```

## Plotting

If the user loads the package [GeoJSON.jl](https://github.com/JuliaGeo/GeoJSON.jl)
together with a [Makie.jl](https://docs.makie.org/stable/) backend, an extension is loaded
and adds the possibility to plot the ground track. In this case, the following functions are
available:

```@docs; canonical = false
plot_ground_track
```

```@docs; canonical = false
plot_ground_track!
```

### Example

The code:

```@repl ground_track
using GeoJSON, CairoMakie

jd₀ = date_to_jd(2021, 1, 1)

orb = KeplerianElements(
    jd₀,
    7130.982e3,
    0.001111,
    98.405 |> deg2rad,
    ltdn_to_raan(10.5, jd₀),
    π / 2,
    0
)

orbp = Propagators.init(Val(:J2), orb)

gt = ground_track(orbp; duration = 5 * 86400, track_types = :ascending)

fig, ax = plot_ground_track(gt)

save("ground_track.png", fig)
```

produces the following figure:

![Ground track](./ground_track.png)
