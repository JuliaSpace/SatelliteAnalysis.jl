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

vgs_ecef = map(vgs_lla) do gs_lla
    geodetic_to_ecef(gs_lla.lat, gs_lla.lon, gs_lla.alt) # WGS84
end

vsat_ecef = map(vsat_lla) do sat_lla
    geodetic_to_ecef(sat_lla.lat, sat_lla.lon, sat_lla.alt) # WGS84
end

@testset "Test New Functions" begin
    # Test that the elevation is actually always 90°
    for i in 1:runs
        R = SatelliteAnalysis._rotmat_ecef_to_enu(vgs_lla[i].lat, vgs_lla[i].lon);
        @test SatelliteAnalysis._get_elevation(vsat_ecef[i], vgs_ecef[i], R) ≈ pi/2
        @test is_ground_facility_visible(vsat_ecef[i], vgs_ecef[i], R, (pi/2) - 1e-6) == true
    end
end

@testset "Test Old Functions" begin
    # Test old Implementation to actually proof the issue related to the LTp and ellipsoid not considered for the elevation computation between satellite and ground station.
    # Due to the issue the elevation will never be 90° (even if it should be) except for the ground station at the poles and the ground station at the equator.
    for i in 1:runs
        if abs(vgs_lla[i].lat) < 1e-6 || abs(abs(vgs_lla[i].lat) - pi/2) < 1e-6
            @test abs(abs(SatelliteAnalysis._get_elevation_old(vsat_ecef[i], vgs_ecef[i]) - pi/2)) < 1e-6
            @test SatelliteAnalysis.is_ground_facility_visible_old(vsat_ecef[i], vgs_ecef[i], (pi/2) - 1e-6) == true
        else
            @test SatelliteAnalysis._get_elevation_old(vsat_ecef[i], vgs_ecef[i]) != pi/2
            @test SatelliteAnalysis.is_ground_facility_visible_old(vsat_ecef[i], vgs_ecef[i], (pi/2) - 1e-6) == false
        end
    end
end