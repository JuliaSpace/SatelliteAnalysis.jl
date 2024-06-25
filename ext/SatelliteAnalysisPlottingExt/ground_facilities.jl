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
    country_polys = GeoMakie.GeoJSON.read(countries_filename)

    # Plot the ground trace.
    fig = Figure(; size = (1450, 800), kwargs...)

    ax = Axis(
        fig[1, 1],
        aspect         = 2,
        title          = "Ground Facility Visibility Circles",
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

    poly!(
        ax,
        GeoMakie.to_multipoly(country_polys.geometry);
        color       = :white,
        strokecolor = :black,
        strokewidth = 1
    )

    # Plot the visibility circles.
    for k in 1:length(vgf_vc)
        gf_vc  = vgf_vc[k - 1 + begin]
        gf_lat = first.(gf_vc)
        gf_lon = last.(gf_vc)

        vc = lines!(ax, gf_lon .|> rad2deg, gf_lat .|> rad2deg; linewidth = 2)

        # We need to compute the vectors in the ECEF reference frame to obtain the ground
        # station position, which is computed by averaging them.
        vr_ecef = geodetic_to_ecef.(gf_lat, gf_lon, 0)
        gf_ecef = sum(vr_ecef) / length(gf_vc)
        gf_lat, gf_lon, ~ = ecef_to_geodetic(gf_ecef)

        dot = scatter!(ax, gf_lon |> rad2deg, gf_lat |> rad2deg; color = vc.color)
        translate!(dot, 0, 0, 10)

        if !isnothing(ground_facility_names)
            label = text!(
                ax,
                ground_facility_names[k];
                color    = vc.color,
                fontsize = _TICK_LABEL_SIZE,
                position = (gf_lon |> rad2deg, gf_lat |> rad2deg),
            )
            translate!(label, 0, 0, 10)
        end
    end

    return fig, ax
end
