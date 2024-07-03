## Description #############################################################################
#
# Functions to plot the ground track.
#
############################################################################################

export plot_ground_track, plot_ground_track!

"""
    plot_ground_track(gt::Vector{NTuple{2, Number}}; kwargs...) -> Figure, Axis

Plot the ground track `gt` computed using the function [`ground_track`](@ref). It returns
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

All `kwargs...` are passed to the function [`plot_world_map`](@ref).

# Extended Help

## Examples

```julia-repl
julia> using SatelliteAnalysis, GeoMakie, GLMakie

julia> jd₀ = date_to_jd(2021, 1, 1)
2.4592155e6

julia> orb = KeplerianElements(
           jd₀,
           7130.982e3,
           0.001111,
           98.405 |> deg2rad,
           ltdn_to_raan(10.5, jd₀),
           π / 2,
           0
       )
KeplerianElements{Float64, Float64}:
           Epoch :    2.45922e6 (2021-01-01T00:00:00)
 Semi-major axis : 7130.98     km
    Eccentricity :    0.001111
     Inclination :   98.405    °
            RAAN :   78.4021   °
 Arg. of Perigee :   90.0      °
    True Anomaly :    0.0      °

julia> orbp = Propagators.init(Val(:J2), orb)
OrbitPropagatorJ2{Float64, Float64}:
   Propagator name : J2 Orbit Propagator
  Propagator epoch : 2021-01-01T00:00:00
  Last propagation : 2021-01-01T00:00:00

julia> gt = ground_track(orbp; track_types = :descending, duration = 5 * 86400);

julia> fig, ax = plot_ground_track(gt; size = (2000, 1000))
(Scene (2000px, 1000px):
  0 Plots
  1 Child Scene:
    └ Scene (2000px, 1000px), Axis (2 plots))

julia> fig
```
"""
function plot_ground_track(::Any; kwargs...)
    error("Wrong input or the package GeoMakie.jl is not loaded.")
    return nothing
end

"""
    plot_ground_track!(ax:Axis, gt::Vector{NTuple{2, Number}}) -> Nothing

Plot in the **Makie.jl** axis `ax` the ground track `gt` computed using the function
[`ground_track`](@ref).

!!! warning

    This function **only works** after loading the package **GeoMakie.jl**. Furthermore, the
    user must also load one Makie.jl back end (CairoMakie.jl or GLMakie.jl, for example) to
    see the result.

# Extended Help

## Examples

```julia-repl
julia> using SatelliteAnalysis, GeoMakie, GLMakie

julia> jd₀ = date_to_jd(2021, 1, 1)
2.4592155e6

julia> orb = KeplerianElements(
           jd₀,
           7130.982e3,
           0.001111,
           98.405 |> deg2rad,
           ltdn_to_raan(10.5, jd₀),
           π / 2,
           0
       )
KeplerianElements{Float64, Float64}:
           Epoch :    2.45922e6 (2021-01-01T00:00:00)
 Semi-major axis : 7130.98     km
    Eccentricity :    0.001111
     Inclination :   98.405    °
            RAAN :   78.4021   °
 Arg. of Perigee :   90.0      °
    True Anomaly :    0.0      °

julia> orbp = Propagators.init(Val(:J2), orb)
OrbitPropagatorJ2{Float64, Float64}:
   Propagator name : J2 Orbit Propagator
  Propagator epoch : 2021-01-01T00:00:00
  Last propagation : 2021-01-01T00:00:00

julia> gt = ground_track(orbp; track_types = :descending, duration = 5 * 86400);


julia> fig = Figure(size = (1000, 1000))

julia> ax = Axis(fig[1, 1])

julia> plot_ground_track!(ax, gt)
```
"""
function plot_ground_track!(args...)
    error("Wrong input or the package GeoMakie.jl is not loaded.")
    return nothing
end
