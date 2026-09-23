## Description #############################################################################
#
# Functions to verify if a ground facility is visible given a satellite position.
#
############################################################################################

export is_ground_facility_visible

"""
    is_ground_facility_visible(sat_r_e::AbstractVector, gf_lat::Number, gf_lon::Number, gf_h::Number, θ::Number) -> Bool
    is_ground_facility_visible(sat_r_e::AbstractVector, gf_r_e::AbstractVector, gf_up_e::AbstractVector, θ::Number) -> Bool
    is_ground_facility_visible(sat_r_e::AbstractVector, gf_r_e::AbstractVector, θ::Number) -> Bool
    is_ground_facility_visible(sat_r_ned::AbstractVector, θ::Number) -> Bool

Check if the satellite with position vector `sat_r_e` (ECEF) is inside the visibility circle
of a ground facility with latitude `gf_lat` [rad], longitude `gf_lon` [rad], altitude `gf_h`
[m] (WGS-84) or ECEF position `gf_r_e` [m]. The algorithm considers that the ground station
has visibility to the satellite if its elevation angle is larger than `θ` [rad].

The user can also pass the satellite position represented in the NED (North-East-Down)
reference frame `sat_r_ned` [m] at the ground station location, which increases the
performance since the algorithm performs no reference frame conversion.

If the visibility of the same ground facility must be verified for many satellite positions,
the fastest method receives the ground facility ECEF position `gf_r_e` [m] together with the
unit vector `gf_up_e` [-] that points to the local vertical (zenith) of the ground facility
represented in the ECEF reference frame. Both can be computed only once:

    gf_r_e  = geodetic_to_ecef(gf_lat, gf_lon, gf_h)
    gf_up_e = [cos(gf_lat) * cos(gf_lon), cos(gf_lat) * sin(gf_lon), sin(gf_lat)]

In this case, the algorithm performs neither reference frame conversions nor trigonometric
operations related to the ground facility position.

# Returns

- `Bool`: `true` if the satellite is inside the visibility circle, or `false` otherwise.
"""
function is_ground_facility_visible(
    sat_r_e::AbstractVector, gf_lat::Number, gf_lon::Number, gf_h::Number, θ::Number
)
    r_ned = ecef_to_ned(sat_r_e, gf_lat, gf_lon, gf_h; translate = true)
    return is_ground_facility_visible(r_ned, θ)
end

function is_ground_facility_visible(
    sat_r_e::AbstractVector, gf_r_e::AbstractVector, gf_up_e::AbstractVector, θ::Number
)
    return _is_ground_facility_visible(sat_r_e, gf_r_e, gf_up_e, sin(θ))
end

function is_ground_facility_visible(
    sat_r_e::AbstractVector, gf_r_e::AbstractVector, θ::Number
)
    gf_wgs84 = ecef_to_geodetic(gf_r_e)
    return is_ground_facility_visible(sat_r_e, gf_wgs84..., θ)
end

function is_ground_facility_visible(r_ned::AbstractVector, minimum_elevation::Number)
    # Check if the satellite is within the minimum elevation supported by the facility.
    # Using the NED vector of the satellite, w.r.t. the current ground facility, it is
    # sufficient to check the angle θ between `r_ned` and the local vertical (-Z axis),
    # then: el = π / 2 - θ.
    z = r_ned[3]
    r = norm(r_ned)
    # Since acos is decreasing, el > minimum_elevation is equivalent to
    # -z / r > cos(π / 2 - minimum_elevation) = sin(minimum_elevation).
    # This avoids the loss of precision (and the relatively expensive acos) near the
    # visibility boundary.
    return -z > r * sin(minimum_elevation)
end

############################################################################################
#                                    Private Functions                                     #
############################################################################################

"""
    _ground_facility_position_and_zenith(gf_lat::Number, gf_lon::Number, gf_h::Number) -> SVector{3, T}, SVector{3, T}

Compute the position [m] of the ground facility with latitude `gf_lat` [rad], longitude
`gf_lon` [rad], and altitude `gf_h` [m] (WGS-84), and the unit vector [-] that points to its
local vertical (zenith), which is normal to the WGS-84 ellipsoid. Both vectors are
represented in the ECEF reference frame, and they are the inputs of the function
`_is_ground_facility_visible`.

# Returns

- `SVector{3, T}`: Ground facility position represented in the ECEF reference frame [m].
- `SVector{3, T}`: Unit vector that points to the local vertical of the ground facility
    represented in the ECEF reference frame [-].
"""
function _ground_facility_position_and_zenith(gf_lat::Number, gf_lon::Number, gf_h::Number)
    gf_r_e = geodetic_to_ecef(gf_lat, gf_lon, gf_h)

    # The local vertical is the axis `-Z` of the NED reference frame, which is defined using
    # the geodetic latitude.
    sin_lat, cos_lat = sincos(gf_lat)
    sin_lon, cos_lon = sincos(gf_lon)

    T = eltype(gf_r_e)
    gf_up_e = SVector{3, T}(cos_lat * cos_lon, cos_lat * sin_lon, sin_lat)

    return gf_r_e, gf_up_e
end

"""
    _is_ground_facility_visible(sat_r_e::AbstractVector, gf_r_e::AbstractVector, gf_up_e::AbstractVector, sin_θ::Number) -> Bool

Check if the satellite with position vector `sat_r_e` [m] is visible from the ground
facility with position `gf_r_e` [m] and local vertical unit vector `gf_up_e` [-], where all
vectors are represented in the ECEF reference frame. `sin_θ` is the sine of the minimum
elevation angle [-].

This function does not perform any trigonometric operation. Hence, it is used in the
algorithms that verify the visibility of the same ground facility many times.
"""
function _is_ground_facility_visible(
    sat_r_e::AbstractVector, gf_r_e::AbstractVector, gf_up_e::AbstractVector, sin_θ::Number
)
    # Vector from the ground facility to the satellite. We do not create a temporary vector
    # here to avoid allocations if the inputs are not static vectors.
    Δx = sat_r_e[1] - gf_r_e[1]
    Δy = sat_r_e[2] - gf_r_e[2]
    Δz = sat_r_e[3] - gf_r_e[3]

    # The component of this vector along the local vertical divided by its norm is the sine
    # of the elevation angle.
    up = Δx * gf_up_e[1] + Δy * gf_up_e[2] + Δz * gf_up_e[3]
    Δr = √(Δx^2 + Δy^2 + Δz^2)

    return up > Δr * sin_θ
end
