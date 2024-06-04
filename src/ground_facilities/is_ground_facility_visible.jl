## Description #############################################################################
#
# Functions to verify is a ground facility is visible given a satellite position.
#
############################################################################################

export is_ground_facility_visible

function _ecef_to_enu_rotmat(lat,lon)
    # Precompute the sines and cosines
    sλ, cλ = sincos(lon)
    sφ, cφ = sincos(lat)
    
    # Generate the rotation matrix as a StaticArray
    # Rotation matrix ECEF -> ENU [https://gssc.esa.int/navipedia/index.php/Transformations_between_ECEF_and_ENU_coordinates]
    return SA_F64[
        -sλ      cλ      0
        -cλ*sφ  -sλ*sφ   cφ
            cλ*cφ   sλ*cφ   sφ
        ]
end

"""
    is_ground_facility_visible(sat_r_e::AbstractVector, gf_r_e::AbstractVector, θ::Number) -> Bool

Check if the satellite with position vector `sat_r_e` (ECEF) is inside the visibility circle
of a ground facility with position vector `gf_r_e` (ECEF) and a minimum elevation angle of
`θ` [rad].

Notice that `sat_r_e` and `gf_r_e` must be represented in the same ECEF frame, and must have
the same unit.

# Returns

- `Bool`: `true` if the satellite is inside the visibility circle, or `false` otherwise.
"""
function is_ground_facility_visible_old(
    sat_r_e::AbstractVector,
    gf_r_e::AbstractVector,
    θ::Number
)
    # //NOTE: It does not seem to account for the ellepsoid even if gf_r_e is computed according to a ref ellipsoid. This seems to lead to a mismatch in the elevation angle computation.
    # Check if the satellite is within the visibility circle of the facility.
    Δr_e = sat_r_e - gf_r_e
    cos_β = dot(Δr_e / norm(Δr_e), gf_r_e / norm(gf_r_e))

    return cos_β > cos(π / 2 - θ)
end
# function is_ground_facility_visible(
#     sat_r_e::AbstractVector,
#     gf_r_e::AbstractVector,
#     vgf_wgs84::Tuple,
#     θ::Number
# )
#     # Check if the satellite is within the visibility circle of the facility.
#     Δr_e = sat_r_e - gf_r_e # Corresponds to Δecef in TelecomUtils.get_visibility(gs,sat)
#     # t = norm(Δr_e)
#     # normalized_Δr_e = Δr_e ./ t
#     # xyz = rv.R' * normalized_Δr_e # Bring Δecef into the ENU frame
#     R = _ecef_to_enu_rotmat(vgf_wgs84[1],vgf_wgs84[2])
#     xyz = R * Δr_e # Bring Δecef into the ENU frame

#     β = acos(dot(xyz / norm(xyz), gf_r_e / norm(gf_r_e)))

#     return (π/2 - β) > θ
# end

# function get_beta_(
#     sat_r_e::AbstractVector,
#     gf_r_e::AbstractVector,
#     vgf_wgs84::Tuple
# )
#     # Check if the satellite is within the visibility circle of the facility.
#     Δr_e = sat_r_e - gf_r_e # Corresponds to Δecef in TelecomUtils.get_visibility(gs,sat)
#     # t = norm(Δr_e)
#     # normalized_Δr_e = Δr_e ./ t
#     # xyz = rv.R' * normalized_Δr_e # Bring Δecef into the ENU frame
#     R = _ecef_to_enu_rotmat(vgf_wgs84[1],vgf_wgs84[2])
#     xyz = R * Δr_e # Bring Δecef into the ENU frame

#     β = acos(dot(xyz / norm(xyz), gf_r_e / norm(gf_r_e)))

#     return rad2deg(β)
# end

# function is_ground_facility_visible( # Same result as old
#     sat_r_e::AbstractVector,
#     gf_r_e::AbstractVector,
#     vgf_wgs84::Tuple,
#     θ::Number
# )
#     # Check if the satellite is within the visibility circle of the facility.
#     Δr_e = sat_r_e - gf_r_e # Corresponds to Δecef in TelecomUtils.get_visibility(gs,sat)
#     # t = norm(Δr_e)
#     # normalized_Δr_e = Δr_e ./ t
#     # xyz = rv.R' * normalized_Δr_e # Bring Δecef into the ENU frame
#     R = _ecef_to_enu_rotmat(vgf_wgs84[1],vgf_wgs84[2])
#     xyz_Δr_e = R * Δr_e # Bring Δecef into the ENU frame
#     xyz_gf_r_e = R * gf_r_e # Bring Δecef into the ENU frame

#     β = acos(dot(xyz_Δr_e / norm(xyz_Δr_e), xyz_gf_r_e / norm(xyz_gf_r_e)))

#     return (π/2 - β) > θ
# end
# function get_el( # Same result as old
#     sat_r_e::AbstractVector,
#     gf_r_e::AbstractVector,
#     vgf_wgs84::Tuple,
# )
#     # Check if the satellite is within the visibility circle of the facility.
#     Δr_e = sat_r_e - gf_r_e # Corresponds to Δecef in TelecomUtils.get_visibility(gs,sat)
#     # t = norm(Δr_e)
#     # normalized_Δr_e = Δr_e ./ t
#     # xyz = rv.R' * normalized_Δr_e # Bring Δecef into the ENU frame
#     R = _ecef_to_enu_rotmat(vgf_wgs84[1],vgf_wgs84[2])
#     xyz_Δr_e = R * Δr_e # Bring Δecef into the ENU frame
#     xyz_gf_r_e = R * gf_r_e # Bring Δecef into the ENU frame

#     β = acos(dot(xyz_Δr_e / norm(xyz_Δr_e), xyz_gf_r_e / norm(xyz_gf_r_e)))

#     return rad2deg(pi/2-β)
# end
function get_el( # Same result as old
    sat_r_e::AbstractVector,
    gf_r_e::AbstractVector,
    vgf_wgs84::Tuple,
)
    # Check if the satellite is within the visibility circle of the facility.
    Δr_e = sat_r_e - gf_r_e # Corresponds to Δecef in TelecomUtils.get_visibility(gs,sat)
    R = _ecef_to_enu_rotmat(vgf_wgs84[1],vgf_wgs84[2])
    xyz = R * Δr_e # Bring Δecef into the ENU frame

    (;x,y,z) = xyz
    r = norm(xyz)
    θ = acos(z/r)
    # φ = atan(y,x)

    return pi/2-θ
end

function is_ground_facility_visible( # Same result as old
    sat_r_e::AbstractVector,
    gf_r_e::AbstractVector,
    vgf_wgs84::Tuple,
    θ::Number
)
    el = get_el(sat_r_e, gf_r_e, vgf_wgs84)

    return el > θ
end
# function is_ground_facility_visible( # Same result as old
#     sat_r_e::AbstractVector,
#     gf_r_e::AbstractVector,
#     vgf_wgs84::Tuple,
#     θ::Number
# )
#     # Check if the satellite is within the visibility circle of the facility.
#     R = _ecef_to_enu_rotmat(vgf_wgs84[1],vgf_wgs84[2])
#     xyz_sat_r_e = R*sat_r_e
#     xyz_gf_r_e = R*gf_r_e
    
#     xyz_Δr_e = xyz_sat_r_e - xyz_gf_r_e

#     β = acos(dot(xyz_Δr_e / norm(xyz_Δr_e), xyz_gf_r_e / norm(xyz_gf_r_e)))

#     return (π/2 - β) > θ
# end

"""
    is_ground_facility_visible(sat_r_e::AbstractVector, gf_lat::Number, gf_lon::Number, gf_h::Number, θ::Number) -> Bool

Check if the satellite with position vector `sat_r_e` (ECEF) is inside the visibility circle
of a ground facility with latitude `gf_lat` [rad], longitude `gf_lon` [rad], altitude `gf_h`
(WGS-84), and a minimum elevation angle of `θ` [rad].

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
    return is_ground_facility_visible(sat_r_e, gf_r_e, θ)
end
