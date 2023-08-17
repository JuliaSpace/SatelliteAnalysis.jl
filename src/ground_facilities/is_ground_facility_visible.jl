# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#   Functions to verify is a ground facility is visible given a satellite position.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export is_ground_facility_visible

"""
    is_ground_facility_visible(sat_r_e::AbstractVector, gs_r_e::AbstractVector, θ::Number) -> Bool

Check if the satellite with position vector `sat_r_e` (ECEF) is inside the visibility circle
of a ground facility with position vector `gs_r_e` (ECEF) and a minimum elevation angle of
`θ` [rad].

Notice that `sat_r_e` and `gs_r_e` must be represented in the same ECEF frame, and must have
the same unit.

# Returns

- `Bool`: `true` if the satellite is inside the visibility circle, or `false` otherwise.
"""
function is_ground_facility_visible(
    sat_r_e::AbstractVector,
    gs_r_e::AbstractVector,
    θ::Number
)
    # Check if the satellite is within the visibility circle of the facility.
    Δr_e = sat_r_e - gs_r_e
    cos_β = dot(Δr_e / norm(Δr_e), gs_r_e / norm(gs_r_e))

    return cos_β > cos(π / 2 - θ)
end

"""
    is_ground_facility_visible(sat_r_e::AbstractVector, gs_lat::Number, gs_lon::Number, gs_h::Number, θ::Number) -> Bool

Check if the satellite with position vector `sat_r_e` (ECEF) is inside the visibility circle
of a ground facility with latitude `gs_lat` [rad], longitude `gs_lon` [rad], altitude `gs_h`
(WGS-84), and a minimum elevation angle of `θ` [rad].

Notice that `sat_r_e` and `gs_h` must have the same units.

# Returns

- `Bool`: `true` if the satellite is inside the visibility circle, or `false` otherwise.
"""
function is_ground_facility_visible(
    sat_r_e::AbstractVector,
    gs_lat::Number,
    gs_lon::Number,
    gs_h::Number,
    θ::Number
)
    gs_r_e = geodetic_to_ecef(gs_lat, gs_lon, gs_h)
    return is_ground_facility_visible(sat_r_e, gs_r_e, θ)
end
