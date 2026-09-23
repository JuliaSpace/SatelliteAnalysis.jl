# SatelliteAnalysis.jl

```@meta
CurrentModule = SatelliteAnalysis
```

This package contains several functions to perform analysis related to satellites. Those
functions were split from the package
[SatelliteToolbox.jl](https://github.com/JuliaSpace/SatelliteToolbox.jl).

## Installation

This package can be installed using:

```julia-repl
julia> using Pkg
julia> Pkg.add("SatelliteAnalysis")
```

## Features

- **Beta angle**: `beta_angle`.
- **Eclipse time**: `eclipse_time_summary` and `lighting_condition`.
- **Ground facilities**: `ground_facility_accesses`, `ground_facility_gaps`,
  `ground_facility_visibility_circle`, and `is_ground_facility_visible`.
- **Ground track**: `ground_track` and `ground_track_inclination`.
- **Orbit design**: `frozen_orbit`, `design_sun_sync_ground_repeating_orbit`,
  `sun_sync_orbit_from_angular_velocity`, `sun_sync_orbit_inclination`,
  `sun_sync_orbit_semi_major_axis`, `ground_repeating_orbit_adjacent_track_angle`, and
  `ground_repeating_orbit_adjacent_track_distance`.
- **Orbital decay**: `decay_analysis`.
- **Plotting**: `plot_world_map`, `plot_ground_track`,
  `plot_ground_facility_visibility_circles`, `plot_decay_analysis`,
  `fetch_country_polygons`, and a Makie theme.

## Optional Features

Some features are provided by package extensions, which are loaded automatically when the
corresponding packages are loaded:

| Load                                    | Features                                                                                            |
|:----------------------------------------|:----------------------------------------------------------------------------------------------------|
| A Makie backend (e.g., **CairoMakie**)  | `makie_theme`, `makie_palette`, and `plot_decay_analysis`                                            |
| A Makie backend and **GeoJSON**         | `plot_world_map`, `plot_ground_track`, and `plot_ground_facility_visibility_circles` (and their `!` variants) |
| **OrdinaryDiffEqAdamsBashforthMoulton** | `decay_analysis` and the macros that select its atmospheric model                                   |
