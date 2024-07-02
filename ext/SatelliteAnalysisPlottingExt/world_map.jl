## Description #############################################################################
#
# Function to create a plot with the world map.
#
############################################################################################

function SatelliteAnalysis.plot_world_map(; size = (1450, 800), kwargs...)
    # Plot the ground trace.
    fig = Figure(; size = size, kwargs...)

    ax = Axis(
        fig[1, 1],
        aspect         = 2,
        title          = "World Map",
        titlegap       = 16,
        titlesize      = _TITLE_SIZE,
        xlabel         = "Longitude [°]",
        xlabelsize     = _LABEL_SIZE,
        xticklabelsize = _TICK_LABEL_SIZE,
        ylabel         = "Latitude [°]",
        ylabelsize     = _LABEL_SIZE,
        yticklabelsize = _TICK_LABEL_SIZE,
    )

    xlims!(ax, -180, +180)
    ylims!(ax, -90,  +90)
    ax.xticks = -180:20:+180
    ax.yticks = -90:20:90

    # Get the GeoJSON file with the countries.
    countries_filename = fetch_country_polygons(; force_download = false)

    # Load the polygons of the countries.
    country_polys = GeoMakie.GeoJSON.read(countries_filename)

    poly!(
        ax,
        GeoMakie.to_multipoly(country_polys.geometry);
        color       = :white,
        strokecolor = :black,
        strokewidth = 1
    )

    return fig, ax
end
