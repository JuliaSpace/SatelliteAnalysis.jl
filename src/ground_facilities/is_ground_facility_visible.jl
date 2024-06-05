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
    gf_r_e = geodetic_to_ecef(gf_lat, gf_lon, gf_h)

    # Check if the satellite is within the minimum elevation supported by the facility.
    Δr_e = sat_r_e - gf_r_e 
    Δr_e_ned = ecef_to_ned(Δr_e, gf_lat, gf_lon, gf_h)
    Δr_e_enu = C_ned_enu * Δr_e_ned
    sph_coord = _cartesian_to_spherical(Δr_e_enu)

    return (π/2-sph_coord.θ) > θ
end

function is_ground_facility_visible(
    sat_r_e::AbstractVector,
    gf_r_e::AbstractVector,
    gf_rot_ecef_enu::SMatrix{3,3},
    θ::Number
)
    # Check if the satellite is within the minimum elevation supported by the facility.
    Δr_e = sat_r_e - gf_r_e 
    Δr_e_enu = gf_rot_ecef_enu * Δr_e
    sph_coord = _cartesian_to_spherical(Δr_e_enu)

    return (π/2-sph_coord.θ) > θ
end

############################################################################################
#                                    Private Functions                                     #
############################################################################################

C_ned_enu = SMatrix{3,3}([0 1 0; 1 0 0; 0 0 -1]) # Conversion matrix NED to ENU

function _cartesian_to_spherical(xyz::AbstractVector)
    (;x,y,z) = xyz
    r = norm(xyz)
    θ = acos(z/r)
    φ = atan(y,x)

    return (;r=r, θ=θ, φ=φ)
end

# Only for test purposes
function _get_elevation(
    sat_r_e::AbstractVector,
    gf_lat::Number,
    gf_lon::Number,
    gf_h::Number
)
    gf_r_e = geodetic_to_ecef(gf_lat, gf_lon, gf_h)

    # Check if the satellite is within the minimum elevation supported by the facility.
    Δr_e = sat_r_e - gf_r_e 
    Δr_e_ned = ecef_to_ned(Δr_e, gf_lat, gf_lon, gf_h)
    Δr_e_enu = C_ned_enu * Δr_e_ned
    sph_coord = _cartesian_to_spherical(Δr_e_enu)

    return (π/2-sph_coord.θ)
end

function _get_elevation(
    sat_r_e::AbstractVector,
    gf_r_e::AbstractVector,
    gf_rot_ecef_enu::SMatrix{3,3}
)
    # Check if the satellite is within the minimum elevation supported by the facility.
    Δr_e = sat_r_e - gf_r_e 
    Δr_e_enu = gf_rot_ecef_enu * Δr_e
    sph_coord = _cartesian_to_spherical(Δr_e_enu)

    return (π/2-sph_coord.θ)
end

# Old function with ellipsoid error
function is_ground_facility_visible_old(
    sat_r_e::AbstractVector,
    gf_r_e::AbstractVector,
    θ::Number
)
    # Check if the satellite is within the visibility circle of the facility.
    Δr_e = sat_r_e - gf_r_e
    cos_β = dot(Δr_e / norm(Δr_e), gf_r_e / norm(gf_r_e))

    return cos_β > cos(π / 2 - θ)
end

function _get_elevation_old(
    sat_r_e::AbstractVector,
    gf_r_e::AbstractVector
)
    # Check if the satellite is within the visibility circle of the facility.
    Δr_e = sat_r_e - gf_r_e
    cos_β = dot(Δr_e / norm(Δr_e), gf_r_e / norm(gf_r_e))

    return π/2 - acos(abs(cos_β) > 1 ? 1 : cos_β)
end

### Functioning 
# """
#     is_ground_facility_visible(sat_r_e::AbstractVector, gf_r_e::AbstractVector, θ::Number) -> Bool

# Check if the satellite with position vector `sat_r_e` (ECEF) is inside the visibility circle
# of a ground facility with position vector `gf_r_e` (ECEF) and a minimum elevation angle of
# `θ` [rad].

# Notice that `sat_r_e` and `gf_r_e` must be represented in the same ECEF frame, and must have
# the same unit.

# # Returns

# - `Bool`: `true` if the satellite is inside the visibility circle, or `false` otherwise.
# """
# function get_el( # Same result as old
#     sat_r_e::AbstractVector,
#     gf_r_e::AbstractVector,
#     vgf_wgs84::Tuple,
# )
#     # Check if the satellite is within the visibility circle of the facility.
#     Δr_e = sat_r_e - gf_r_e # Corresponds to Δecef in TelecomUtils.get_visibility(gs,sat)
#     R = _ecef_to_enu_rotmat(vgf_wgs84[1],vgf_wgs84[2])
#     xyz = R * Δr_e # Bring Δecef into the ENU frame

#     (;x,y,z) = xyz
#     r = norm(xyz)
#     θ = acos(z/r)
#     # φ = atan(y,x)

#     return pi/2-θ
# end

# function is_ground_facility_visible( # Same result as old
#     sat_r_e::AbstractVector,
#     gf_r_e::AbstractVector,
#     vgf_wgs84::Tuple,
#     θ::Number
# )
#     el = get_el(sat_r_e, gf_r_e, vgf_wgs84)

#     return el > θ
# end

# function _ecef_to_enu_rotmat(lat,lon)
#     # Precompute the sines and cosines
#     sλ, cλ = sincos(lon)
#     sφ, cφ = sincos(lat)
    
#     # Generate the rotation matrix as a StaticArray
#     # Rotation matrix ECEF -> ENU [https://gssc.esa.int/navipedia/index.php/Transformations_between_ECEF_and_ENU_coordinates]
#     return SA_F64[
#         -sλ      cλ      0
#         -cλ*sφ  -sλ*sφ   cφ
#             cλ*cφ   sλ*cφ   sφ
#         ]
# end
####

