# Ground Track

```@meta
CurrentModule = SatelliteAnalysis
```

```@setup ground_track
using SatelliteAnalysis
```

We can obtain the ground track of a satellite using the function:

```julia
ground_track(orbp::OrbitPropagator; kwargs...) -> Vector{NTuple{2, Float64}}
```

It computes the satellite ground track using the orbit propagator `orbp`. It returns a
vector of `NTuple{2, Float64}` where the first element is the latitude [rad] and the second
is the longitude [rad] of each point in the ground track.

The following keywords are available:

- `add_nans::Bool`: If `true`, we add `NaN` if there is a discontinuity in the ground track
    to improve plotting.
    (**Default**: true)
- `duration::Number`: Duration of the analysis.
    (**Default**: 86400)
- `initial_time::Number`: Initial time regarding the orbit propagator `orbp` epoch [s].
    (**Default**: 0)
- `f_eci_to_ecef::Function`: Function to convert the orbit propagator position represented
    in the Earth-centered inertial (ECI) reference frame to the Earth-centered, Earth-fixed
    (ECEF) reference frame. The signature must be

    `f_eci_to_ecef(r_i::AbstractVector, jd::Number) -> AbstractVector`

    and it must return the position vector `r_i` represented in the ECEF at the instant `jd`
    [Julian Day]. By default, we use TEME as the ECI and PEF as the ECEF.
    (**Default**: `_ground_track_default_eci_to_ecef`)
- `step::Union{Nothing, Number}`: Step for the computation. If `nothing`, we will roughly
    compute the step to approximate 1° in the mean anomaly.
    (**Default**: `nothing`)
- `track_types::Symbol`: A symbol describing what kind of track types we must add to the
    output vector. It can be `:ascending` for only ascending passages, `:descending` for
    only descending passages, or `:all` for both.
    (**Default**: `:all`)

## Ground Track Inclination

Many analyses require the ground track inclination. For example, if we are designing a
remote sensing mission with an optical payload, we must know how much two images overlap.
This information can only be computed using spherical trigonometry and the ground track
inclination instead of the orbital one.

The ground track inclination ``i_{gt}``, shown in the following figure, is a composition of the orbit inclination, the
Earth's angular speed, and the RAAN time derivative. We can compute it using:

```math
i_{gt} = \tan^{-1}\left(\frac{
    \omega_s \sin i 
}{
    \omega_s \cos i - \omega_e + \dot{\Omega}
}\right)\ ,
```

where ``i`` is the orbital inclination, ``\omega_s`` is the satellite angular speed at
Equator, ``\omega_e`` is the Earth angular speed, and ``\Omega`` is the RAAN.

```@raw html
<div align="center">
  <img src="../../assets/ground_track_inclination.png" alt="Ground Track Inclination" width="100%"/>
</div>
```

Formally, we should use the satellite instantaneous angular speed at the Equator in
``\omega_s``. However, given the perturbations caused by the Earth's gravitational
potential, this speed is not simple to compute. The calculation would required to implement
an orbit propagator. Thus, we simplify it by assuming that the orbit eccentricity is small.
This assumption is reasonable given the missions that would benefit from the computation of
the ground track inclination. In this case, we approximate ``\omega_s`` as the mean
satellite angular speed.

The function:

```julia
ground_track_inclination(a::Number, e::Number, i::Number; kwargs...) -> T
ground_track_inclination(orb::Orbit{Tepoch, T}); kwargs...) where {Tepoch <: Number, T <: Number} -> T
```

computes the ground track inclination at the Equator [rad] in an orbit with semi-major axis
`a` [m], eccentricity `e` [ ], and inclination `i` [rad]. The orbit can also be specified by
`orb` (see `Orbit`).

!!! note
    The output type `T` in the first signature is obtained by promoting the inputs to a
    float type.

The following keywords are available:

- `perturbation::Symbol`: Symbol to select the perturbation terms that will be used. It can
    be `:J0`, `:J2`, or `:J4`.
    (**Default**: `:J2`)
- `m0::Number`: Standard gravitational parameter for Earth [m³ / s²].
    (**Default**: `GM_EARTH`)
- `J2::Number`: J₂ perturbation term.
    (**Default**: `EGM_2008_J2`)
- `J4::Number`: J₄ perturbation term.
    (**Default**: `EGM_2008_J4`)
- `R0::Number`: Earth's equatorial radius [m].
    (**Default**: `EARTH_EQUATORIAL_RADIUS`)
- `we::Number`: Earth's angular speed [rad / s].
    (**Default**: `EARTH_ANGULAR_SPEED`)

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

If we use **GeoMakie.jl** to plot, we obtain:

```@setup ground_track
using GeoMakie, CairoMakie

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

If the user loads the package [GeoMakie.jl](https://github.com/MakieOrg/GeoMakie.jl)
together with a [Makie.jl](https://docs.makie.org/stable/) back end, an extension is loaded
and adds the possibility to plot the ground track. In this case, the following function is
available:

```julia
plot_ground_track(gt::Vector{NTuple{2, T}}; kwargs...) where T<:Number -> Figure, Axis
```

It plots the ground track `gt` computed using the function [`ground_track`](@ref). It
returns the objects `Figure` and `Axis` used to plot the data. For more information, please,
refer to [Makie.jl](https://docs.makie.org/stable/) documentation.

All `kwargs...` are passed to the function `Figure`.

### Example

The code:

```@repl ground_track
using GeoMakie, CairoMakie

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
