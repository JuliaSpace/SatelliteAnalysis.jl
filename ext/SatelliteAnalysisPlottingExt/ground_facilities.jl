## Description #############################################################################
#
# Plot results related to the ground facilities.
#
############################################################################################

function SatelliteAnalysis.plot_ground_facility_visibility_circles!(
    ax::Axis,
    vgf_vc::Vector{Vector{NTuple{2, T}}};
    ground_facility_names::Union{Nothing, Vector{String}} = nothing
) where T<:Number
    # Check inputs.
    if !isnothing(ground_facility_names) && (length(vgf_vc) != length(ground_facility_names))
        throw(ArgumentError(
            "The number of elements in `vgf_vc` and `ground_facility_names` must be equal."
        ))
    end

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

    return nothing
end

function SatelliteAnalysis.plot_ground_facility_visibility_circles(
    vgf_vc::Vector{Vector{NTuple{2, T}}};
    ground_facility_names::Union{Nothing, Vector{String}} = nothing,
    kwargs...
) where T <: Number

    fig, ax = plot_world_map(; kwargs...)

    ax.title = "Ground Facility Visibility Circles"

    plot_ground_facility_visibility_circles!(
        ax,
        vgf_vc;
        ground_facility_names = ground_facility_names
    )

    return fig, ax
end
