# World Map

```@meta
CurrentModule = SatelliteAnalysis
```

```@setup world_map
using SatelliteAnalysis
```

**SatelliteAnalysis.jl** has a built-in algorithm to plot the World Map provided that the
user loaded **GeoMakie.jl** and one of the **Makie.jl** backends. This empty plot can be
used to add analysis on top of it. We can create it using the function:

```julia
plot_world_map(; kwargs...) -> Figure, Axis
```

It returns a **Makie.jl** `Figure` and `Axis` with the World map. The figure is styled with
the theme obtained from the function `SatelliteAnalysis.makie_theme`, selected by the
keyword `theme`, which can be `:light` (default) or `:dark`. All other `kwargs...` are
passed to the function `Figure`.

!!! note

    This function plots the countries' borders in the created figure using the file with the
    country polygons fetched with the function [`fetch_country_polygons`](@ref). Hence, if
    this file does not exist, the algorithm tries to download it.

```@repl world_map
using GeoMakie, CairoMakie

fig, ax = plot_world_map()
fig
```

```@setup world_map
using GeoMakie, CairoMakie

fig, ax = plot_world_map()
save("world_map.png", fig)
```

![World Map](./world_map.png)

We can also create the World map using the dark theme variant:

```@setup world_map
fig_dark, ax_dark = plot_world_map(; theme = :dark)
save("world_map_dark.png", fig_dark)
```

```julia
fig, ax = plot_world_map(; theme = :dark)
```

![World Map (Dark)](./world_map_dark.png)
