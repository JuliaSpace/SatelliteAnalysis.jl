## Description #############################################################################
#
# Function to create a plot with the world map.
#
############################################################################################

function SatelliteAnalysis.plot_world_map(;
    theme::Symbol = :light,
    size = (1450, 800),
    kwargs...
)
    # Build the theme first since it also validates the variant in `theme`.
    sa_theme = SatelliteAnalysis.makie_theme(theme)

    # Every object must be created inside `with_theme` because Makie resolves the theme
    # attributes at object-creation time.
    return with_theme(sa_theme) do
        fig = Figure(; size = size, kwargs...)

        ax = Axis(
            fig[1, 1];
            aspect = 2,
            title  = "World Map",
            xlabel = "Longitude [°]",
            ylabel = "Latitude [°]",
        )

        xlims!(ax, -180, +180)
        ylims!(ax, -90, +90)
        ax.xticks = -180:20:+180
        ax.yticks = -90:15:90

        # Get the GeoJSON file with the countries.
        countries_filename = fetch_country_polygons(; force_download = false)

        # Load the polygons of the countries.
        country_polys = GeoMakie.GeoJSON.read(countries_filename)

        poly!(
            ax,
            GeoMakie.to_multipoly(country_polys.geometry);
            color       = _COUNTRY_FILL_COLOR[theme],
            strokecolor = _COUNTRY_STROKE_COLOR[theme],
            strokewidth = 1,
        )

        return fig, ax
    end
end
