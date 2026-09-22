## Description #############################################################################
#
# Plot ground track.
#
############################################################################################

function SatelliteAnalysis.plot_ground_track!(
    ax::Axis, gt::Vector{NTuple{2, T}}; kwargs...
) where {T <: Number}
    gt_lat = first.(gt) .|> rad2deg
    gt_lon = last.(gt) .|> rad2deg

    return lines!(ax, gt_lon, gt_lat; kwargs...)
end

function SatelliteAnalysis.plot_ground_track(
    gt::Vector{NTuple{2, T}}; theme::Union{Nothing, Symbol, Makie.Theme} = :light, kwargs...
) where {T <: Number}
    # Wrap the entire body in this function so that the plot calls also resolve their
    # attributes, such as the color cycle, using the selected theme.
    return SatelliteAnalysis._with_plot_theme(theme) do
        fig, ax = _create_world_map(_theme_variant(theme); kwargs...)
        ax.title = "Ground Track"

        plot_ground_track!(ax, gt)

        return fig, ax
    end
end
