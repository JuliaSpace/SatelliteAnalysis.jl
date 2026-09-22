# World Map

```@meta
CurrentModule = SatelliteAnalysis
```

```@setup world_map
using SatelliteAnalysis
```

**SatelliteAnalysis.jl** has a built-in algorithm to plot the World Map provided that the
user loaded **GeoJSON.jl** and one of the **Makie.jl** backends. This empty plot can be
used to add analysis on top of it. We can create it using the function `plot_world_map`:

```@docs; canonical = false
plot_world_map
```

```@docs; canonical = false
plot_world_map!
```

```@repl world_map
using GeoJSON, CairoMakie

fig, ax = plot_world_map();
```

```@setup world_map
using GeoJSON, CairoMakie

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
