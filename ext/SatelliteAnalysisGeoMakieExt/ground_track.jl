## Description #############################################################################
#
# Plot ground track.
#
############################################################################################

function SatelliteAnalysis.plot_ground_track!(
    ax::Axis, gt::Vector{NTuple{2, T}}
) where {T <: Number}
    gt_lat = first.(gt) .|> rad2deg
    gt_lon = last.(gt) .|> rad2deg

    lines!(ax, gt_lon, gt_lat; linewidth = 2)

    return nothing
end

function SatelliteAnalysis.plot_ground_track(
    gt::Vector{NTuple{2, T}};
    theme::Symbol = :light,
    kwargs...
) where {T <: Number}
    # Wrap the entire body in `with_theme` so that the plot calls also resolve their
    # attributes, such as the color cycle, using the SatelliteAnalysis.jl theme.
    return with_theme(SatelliteAnalysis.makie_theme(theme)) do
        fig, ax = plot_world_map(; theme = theme, kwargs...)
        ax.title = "Ground Track"

        plot_ground_track!(ax, gt)

        return fig, ax
    end
end
