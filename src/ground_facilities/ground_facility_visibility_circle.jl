## Description #############################################################################
#
# Functions to compute the ground facility visibility circle.
#
############################################################################################

export ground_facility_visibility_circle

"""
    ground_facility_visibility_circle(gf_wgs84::Tuple, satellite_position_norm::Number; kwargs...) -> Vector{NTuple{2, Float64}}

Compute the ground facility visibility circle from the position `gf_wgs84` (WGS84) to a
satellite in which its distance from the Earth's center is `satellite_position_norm` [m]. It
returns a vector of `NTuple{2, Float64}` where the first element is the latitude [rad] and
the second is the longitude [rad] of each point in the visibility circle.

The ground facility is specified using a tuple with its WGS84 position:

    (latitude [rad], longitude [rad], altitude [m])

# Keywords

- `azimuth_step::Number`: The step in the azimuth used to compute the visibility circle.
    (**Default**: 0.1 |> deg2rad)
- `minimum_elevation::Number`: Minimum elevation angle for communication between the
    satellite and the ground facility [rad].
    (**Default**: 10 |> deg2rad)

# Extended Help

## Examples

```julia-repl
julia> using SatelliteAnalysis, UnicodePlots

julia> gfv = ground_facility_visibility_circle((0, 0, 0), EARTH_EQUATORIAL_RADIUS + 700e3);

julia> lineplot(last.(gfv) .|> rad2deg, first.(gfv) .|> rad2deg; xlim = (-180, 180), ylim = (-90, 90))
       ┌────────────────────────────────────────┐ 
    90 │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
       │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
       │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
       │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
       │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⡴⠚⠉⡏⠓⢦⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
       │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⡞⠀⠀⠀⡇⠀⠀⢳⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
       │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢸⠀⠀⠀⠀⡇⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
       │⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⢼⠤⠤⠤⠤⡧⠤⠤⠤⡧⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤│ 
       │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢸⠀⠀⠀⠀⡇⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
       │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⢧⠀⠀⠀⡇⠀⠀⡼⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
       │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠳⢤⣀⣇⡤⠞⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
       │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
       │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
       │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
   -90 │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
       └────────────────────────────────────────┘ 
       ⠀-180⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀180⠀ 

julia> gfv = ground_facility_visibility_circle((-40 |> deg2rad, -60 |> deg2rad, 700), EARTH_EQUATORIAL_RADIUS + 700e3);

julia> lineplot(last.(gfv) .|> rad2deg, first.(gfv) .|> rad2deg; xlim = (-180, 180), ylim = (-90, 90))
       ┌────────────────────────────────────────┐ 
    90 │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
       │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
       │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
       │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
       │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
       │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
       │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
       │⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⢤⡤⠤⠶⠤⣤⠤⠤⠤⠤⡧⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤│ 
       │⠀⠀⠀⠀⠀⠀⠀⠀⠀⡴⠋⠀⠀⠀⠀⠈⠳⣄⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
       │⠀⠀⠀⠀⠀⠀⠀⠀⡞⠁⠀⠀⠀⠀⠀⠀⠀⠘⣆⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
       │⠀⠀⠀⠀⠀⠀⠀⡼⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠘⡆⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
       │⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢳⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
       │⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣸⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
       │⠀⠀⠀⠀⠀⠀⠀⠳⢤⣀⠀⠀⠀⠀⠀⠀⢀⣀⡴⠃⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
   -90 │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠉⠉⠉⠉⠉⠉⠉⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
       └────────────────────────────────────────┘ 
       ⠀-180⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀180⠀ 
```
"""
function ground_facility_visibility_circle(
    gf_wgs84::Tuple{T1, T2, T3},
    satellite_position_norm::Number;
    azimuth_step::Number = 0.1 |> deg2rad,
    minimum_elevation::Number = 10 |> deg2rad,
) where {T1 <: Number, T2 <: Number, T3 <: Number}

    gf_lat, gf_lon, gf_h = gf_wgs84

    # Vector with the azimuth angles in the ground facility that we will consider in the
    # analysis.
    vazimuth = float(-π):float(azimuth_step):float(π + azimuth_step)

    # Obtain the ground facility position in the ECEF reference frame.
    gf_ecef = geodetic_to_ecef(gf_lat, gf_lon, gf_h)

    # Distance between the Earth's center and the ground station.
    r_gf = norm(gf_ecef)

    # Compute the distance `r` from the ground facility to the satellite considering the
    # minimum elevation angle. We obtained this equation from the law of cosines.
    #
    # NOTE: Only one of the two answers make sense becanse `r` must be positive.
    sin_θ, cos_θ = sincos(minimum_elevation)
    r = -r_gf * sin_θ + √(satellite_position_norm^2 - (r_gf * cos_θ)^2)

    # Allocate the vector with the visibility circle.
    gf_visibility = NTuple{2, Float64}[]
    sizehint!(gf_visibility, length(vazimuth))

    for az in vazimuth
        # Given the local reference frame (NED) in the ground facility, compute the
        # direction to the intersection of the visibiliy region considering the minimum
        # elevation angle and the current azimuth.
        D_l_ned = angle_to_dcm(az, minimum_elevation, :ZY)
        D_ned_l = D_l_ned'
        r_ned   = D_ned_l * @SVector([r_gf, 0, 0])

        # Convert the position vector to ECEF and obtain the latitude and longitude of the
        # intersection point.
        r_ecef      = ned_to_ecef(r_ned, gf_lat, gf_lon, gf_h; translate = true)
        lat, lon, ~ = ecef_to_geodetic(r_ecef)

        # Add to the result vector.
        push!(gf_visibility, (lat, lon))
    end

    return gf_visibility
end
