## Description #############################################################################
#
# Plot ground track.
#
############################################################################################

function SatelliteAnalysis.plot_ground_track(
    gt::Vector{NTuple{2, T}};
    kwargs...
) where T<:Number
    gt_lat = first.(gt) .|> rad2deg
    gt_lon = last.(gt)  .|> rad2deg

    # Get the GeoJSON file with the countries.
    countries_filename = fetch_country_polygons(; force_download = false)

    # Load the polygons of the countries.
    country_polys = GeoMakie.GeoJSON.read(read(countries_filename))

    # Plot the ground trace.
    fig = Figure(; size = (1450, 800), kwargs...)

    ax = Axis(
        fig[1, 1],
        aspect = 2,
        xlabel = "Longitude [°]",
        xlabelsize = _LABEL_SIZE,
        xticklabelsize = _TICK_LABEL_SIZE,
        ylabel = "Latitude [°]",
        ylabelsize = _LABEL_SIZE,
        yticklabelsize = _TICK_LABEL_SIZE,
        title = "Ground Track",
        titlegap = 16,
        titlesize = _TITLE_SIZE,
    )

    xlims!(ax, -180, +180)
    ylims!(ax, -90,  +90)
    ax.xticks = -180:20:+180
    ax.yticks = -90:20:90

    poly!(
        ax,
        country_polys;
        color = :white,
        strokecolor = :black,
        strokewidth = 1
    )

    lines!(ax, gt_lon, gt_lat; linewidth = 2)

    return fig, ax
end
