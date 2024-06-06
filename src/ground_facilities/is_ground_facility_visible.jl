## Description #############################################################################
#
# Functions to verify is a ground facility is visible given a satellite position.
#
############################################################################################

export is_ground_facility_visible

"""
    is_ground_facility_visible(sat_r_e::AbstractVector, gf_lat::Number, gf_lon::Number, gf_h::Number, θ::Number) -> Bool
    is_ground_facility_visible(sat_r_e::AbstractVector, gf_r_e::AbstractVector, gf_rot_ecef_enu::SMatrix{3,3}, θ::Number) -> Bool

Check if the satellite with position vector `sat_r_e` (ECEF) is inside the visibility circle
of a ground facility with latitude `gf_lat` [rad], longitude `gf_lon` [rad], altitude `gf_h`
(WGS-84), and a minimum elevation angle of `θ` [rad].

The second method is an alternative to the main function allowing to pre-calculate the rotation matrix for each of the ground facility location. This avoids computing the ECEF to Local Coordinates Plane at every propagation time step.

Notice that `sat_r_e` and `gf_h` must have the same units.

# Returns

- `Bool`: `true` if the satellite is inside the visibility circle, or `false` otherwise.
"""
function is_ground_facility_visible(
    sat_r_e::AbstractVector,
    gf_lat::Number,
    gf_lon::Number,
    gf_h::Number,
    θ::Number
)
    r_ned = ecef_to_ned(sat_r_e, gf_lat, gf_lon, gf_h; translate = true)
    return is_ground_facility_visible(r_ned, θ)
end

function is_ground_facility_visible(
    r_ned::AbstractVector,
    minimum_elevation::Number
)
    # Check if the satellite is within the minimum elevation supported by the facility.
    # Using the NED vector of the satellite, wrt to the current ground facility, it is sufficient 
    # to check the angle θ between `r_ned` and the local vertical (-Z axis), then: el = π/2 - θ.
    (;x,y,z) = r_ned
    r = norm(r_ned)
    θ = acos(-z/r)
    
    return (π/2 - θ) > minimum_elevation
end