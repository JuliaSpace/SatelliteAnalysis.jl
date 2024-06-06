# Additional Functions# Only for test purposes
function _get_elevation(
    sat_r_e::AbstractVector,
    gf_lat::Number,
    gf_lon::Number,
    gf_h::Number,
)
    r_ned = ecef_to_ned(sat_r_e, gf_lat, gf_lon, gf_h; translate = true)
    return _get_elevation(r_ned)
end

function _get_elevation( # //FIX: to be removed
    r_ned::AbstractVector,
)
    # Check if the satellite is within the minimum elevation supported by the facility.
    (;x,y,z) = r_ned
    r = norm(r_ned)
    θ = acos(-z/r)
    
    return (π/2 - θ)
end

# Initialization
runs = 100000
prec = pi/200
lats = rand(-pi/2:prec:pi/2, runs)
lons = rand(-pi:prec:pi, runs)
alts = rand(250e3:8000e3, runs)

vgs_lla = map(lats,lons) do lat,lon
    (lat=lat, lon=lon, alt=0.0)
end

vsat_lla = map(lats,lons,alts) do lat,lon,alt
    (lat=lat, lon=lon, alt=alt)
end



vsat_ecef = map(vsat_lla) do sat_lla
    geodetic_to_ecef(sat_lla.lat, sat_lla.lon, sat_lla.alt) # WGS84
end

@testset "Test New Functions" begin
    # Test that the elevation is actually always 90°
    for i in 1:runs
        @test _get_elevation(vsat_ecef[i], vgs_lla[i][1], vgs_lla[i][2], vgs_lla[i][3]) ≈ pi/2
        @test is_ground_facility_visible(vsat_ecef[i], vgs_lla[i][1], vgs_lla[i][2], vgs_lla[i][3], (pi/2) - 1e-6) == true
    end
end

# @testset "Test Old Functions" begin
#     # Test old Implementation to actually proof the issue related to the LTp and ellipsoid not considered for the elevation computation between satellite and ground station.
#     # Due to the issue the elevation will never be 90° (even if it should be) except for the ground station at the poles and the ground station at the equator.
    
#     # Functions
#     function is_ground_facility_visible_old(
#         sat_r_e::AbstractVector,
#         gf_r_e::AbstractVector,
#         θ::Number
#     )
#         # Check if the satellite is within the visibility circle of the facility.
#         Δr_e = sat_r_e - gf_r_e
#         cos_β = dot(Δr_e / norm(Δr_e), gf_r_e / norm(gf_r_e))

#         return cos_β > cos(π / 2 - θ)
#     end

#     function _get_elevation_old(
#         sat_r_e::AbstractVector,
#         gf_r_e::AbstractVector
#     )
#         # Check if the satellite is within the visibility circle of the facility.
#         Δr_e = sat_r_e - gf_r_e
#         cos_β = dot(Δr_e / norm(Δr_e), gf_r_e / norm(gf_r_e))

#         return π/2 - acos(abs(cos_β) > 1 ? 1 : cos_β)
#     end

#     vgs_ecef = map(vgs_lla) do gs_lla
#         geodetic_to_ecef(gs_lla.lat, gs_lla.lon, gs_lla.alt) # WGS84
#     end
#     for i in 1:runs
#         if abs(vgs_lla[i].lat) < 1e-6 || abs(abs(vgs_lla[i].lat) - pi/2) < 1e-6
#             @test abs(abs(_get_elevation_old(vsat_ecef[i], vgs_ecef[i]) - pi/2)) < 1e-6
#             @test is_ground_facility_visible_old(vsat_ecef[i], vgs_ecef[i], (pi/2) - 1e-6) == true
#         else
#             @test _get_elevation_old(vsat_ecef[i], vgs_ecef[i]) != pi/2
#             @test is_ground_facility_visible_old(vsat_ecef[i], vgs_ecef[i], (pi/2) - 1e-6) == false
#         end
#     end
# end