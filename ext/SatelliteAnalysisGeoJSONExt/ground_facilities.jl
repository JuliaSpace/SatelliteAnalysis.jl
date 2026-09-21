## Description #############################################################################
#
# Plot results related to the ground facilities.
#
############################################################################################

function SatelliteAnalysis.plot_ground_facility_visibility_circles!(
    ax::Axis,
    vgf_vc::Vector{Vector{NTuple{2, T}}};
    ground_facility_names::Union{Nothing, Vector{String}} = nothing,
    kwargs...,
) where {T <: Number}
    # Check inputs.
    if !isnothing(ground_facility_names) &&
        (length(vgf_vc) != length(ground_facility_names))
        throw(
            ArgumentError(
                "The number of elements in `vgf_vc` and `ground_facility_names` must be equal.",
            ),
        )
    end

    # Vector with the plots of the visibility circles, which is returned to the user.
    plots = Lines[]

    # Plot the visibility circles.
    for (k, gf_vc) in enumerate(vgf_vc)
        gf_lat = first.(gf_vc)
        gf_lon = last.(gf_vc)

        vc = lines!(ax, gf_lon .|> rad2deg, gf_lat .|> rad2deg; kwargs...)
        push!(plots, vc)

        # We need to compute the vectors in the ECEF reference frame to obtain the ground
        # station position, which is computed by averaging them. Notice that we must neglect
        # the NaNs that indicate the discontinuities in the visibility circle.
        valid_points = Iterators.filter(p -> !isnan(first(p)), gf_vc)

        gf_ecef =
            sum(p -> geodetic_to_ecef(first(p), last(p), 0), valid_points) /
            count(p -> !isnan(first(p)), gf_vc)

        center_lat, center_lon, ~ = ecef_to_geodetic(gf_ecef)

        dot = scatter!(ax, center_lon |> rad2deg, center_lat |> rad2deg; color = vc.color)
        translate!(dot, 0, 0, 10)

        if !isnothing(ground_facility_names)
            label = text!(
                ax,
                ground_facility_names[k];
                color    = vc.color,
                position = (center_lon |> rad2deg, center_lat |> rad2deg),
            )
            translate!(label, 0, 0, 10)
        end
    end

    return plots
end

function SatelliteAnalysis.plot_ground_facility_visibility_circles(
    vgf_vc::Vector{Vector{NTuple{2, T}}};
    ground_facility_names::Union{Nothing, Vector{String}} = nothing,
    theme::Union{Nothing, Symbol, Makie.Theme} = :light,
    kwargs...,
) where {T <: Number}
    # Wrap the entire body in this function so that the plot calls also resolve their
    # attributes, such as the color cycle, using the selected theme.
    return SatelliteAnalysis._with_plot_theme(theme) do
        fig, ax = _create_world_map(_theme_variant(theme); kwargs...)

        ax.title = "Ground Facility Visibility Circles"

        plot_ground_facility_visibility_circles!(
            ax, vgf_vc; ground_facility_names = ground_facility_names
        )

        return fig, ax
    end
end
