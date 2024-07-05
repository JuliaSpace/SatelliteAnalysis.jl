# World Map

```@meta
CurrentModule = SatelliteAnalysis
```

```@setup world_map
using SatelliteAnalysis
```

**SatelliteAnalysis.jl** has a built-in algorithm to plot the World Map provided that the
user loaded **GeoMakie.jl** and one of the **Makie.jl**'s back end. This empty plot can be
used to add analysis on top of it. We can create it using the function:

```julia
plot_world_map(; kwargs...) -> Figure, Axis
```

It returns a **Makie.jl** `Figure` and `Axis` with the World map. All `kwargs...` are passed
to the function `Figure`.

!!! note

    This function plots the countries' borders in the created figure using the file with the
    country polygons fetched with the function [`fetch_country_polygons`](@ref). Hence, if
    this files does not exist, the algorithm tries to download it.

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
