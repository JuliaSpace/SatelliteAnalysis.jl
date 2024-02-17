## Description #############################################################################
#
# Plot results related to the ground facilities.
#
############################################################################################

function SatelliteAnalysis.plot_ground_facility_visibility_circles(
    vgf_vc::Vector{Vector{NTuple{2, T}}};
    ground_facility_names::Union{Nothing, Vector{String}} = nothing,
    kwargs...
) where T <: Number
    # Check inputs.
    if !isnothing(ground_facility_names) && (length(vgf_vc) != length(ground_facility_names))
        throw(ArgumentError(
            "The number of elements in `vgf_vc` and `ground_facility_names` must be equal."
        ))
    end

    # Get the GeoJSON file with the countries.
    countries_filename = fetch_country_polygons(; force_download = false)

    # Load the polygons of the countries.
    country_polys = GeoMakie.GeoJSON.read(read(countries_filename))

    # Plot the ground trace.
    fig = Figure(; size = (1600, 800), kwargs...)

    ax = Axis(
        fig[1, 1],
        aspect = 2,
        xlabel = "Longitude [°]",
        xlabelsize = _LABEL_SIZE,
        xticklabelsize = _TICK_LABEL_SIZE,
        ylabel = "Latitude [°]",
        ylabelsize = _LABEL_SIZE,
        yticklabelsize = _TICK_LABEL_SIZE,
        title = "Ground Facility Visibility Circles",
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

    # Plot the visibility circles.
    for k in 1:length(vgf_vc)
        gf_vc  = vgf_vc[k - 1 + begin]
        gf_lat = first.(gf_vc) .|> rad2deg
        gf_lon = last.(gf_vc)  .|> rad2deg

        lines!(ax, gf_lon, gf_lat; linewidth = 2)

        gf_lat_mean = sum(gf_lat) ./ length(gf_lat)
        gf_lon_mean = sum(gf_lon) ./ length(gf_lon)

        scatter!(ax, gf_lon_mean, gf_lat_mean)

        !isnothing(ground_facility_names) && text!(
            ax,
            ground_facility_names[k];
            position = (gf_lon_mean, gf_lat_mean),
            fontsize = _TICK_LABEL_SIZE,
        )
    end

    return fig, ax
end
