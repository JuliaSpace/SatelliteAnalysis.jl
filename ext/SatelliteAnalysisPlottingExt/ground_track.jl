## Description #############################################################################
#
# Plot ground track.
#
############################################################################################

function SatelliteAnalysis.plot_ground_track(
    gt::Vector{NTuple{2, T}};
    kwargs...
) where T<:Number

    fig, ax = plot_world_map(; kwargs...)
    ax.title = "Ground Track"

    gt_lat = first.(gt) .|> rad2deg
    gt_lon = last.(gt)  .|> rad2deg

    lines!(ax, gt_lon, gt_lat; linewidth = 2)

    return fig, ax
end
