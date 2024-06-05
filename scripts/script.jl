using SatelliteAnalysis
using TelecomUtils

epoch = DateTime("2024-01-01")
gs = (deg2rad(30),0.0,0.0)
mea = deg2rad(10)
step = 30
h = 500e3
i = deg2rad(80)

# Reproduce error with TEME-PEF
# acc1,gap1 = let
	
# 	orbElem = KeplerianElements(
# 				epoch |> datetime2julian, # t
# 				7878.14e3, # a
# 				0.0, # e
# 				deg2rad(60), # i
# 				0.0, # Ω
# 				0.0, # ω
# 				0.0 # f
# 			)
	
# 	orbp = Propagators.init(Val(:J2), orbElem)

# 	acc = ground_facility_accesses(orbp, gs; minimum_elevation=deg2rad(60), duration=24*3600, step=step);

# 	gap = ground_facility_gaps(orbp, gs; minimum_elevation=deg2rad(60), duration=24*3600, step=step);

# 	acc,gap
# end;

# Reproduce error with J2000-ITRF
acc2,gap2 = let
	eop_iau1980 = fetch_iers_eop(Val(:IAU1980));
	
	function _eci_to_ecef(r_i::AbstractVector, jd::Number)
		R = r_eci_to_ecef(J2000(), ITRF(), jd, eop_iau1980)
		return R*r_i
	end
	
	# orbElem = KeplerianElements(
	# 			epoch |> datetime2julian, # t
	# 			6878.14e3, # a
	# 			0.0, # e
	# 			deg2rad(20), # i
	# 			0.0, # Ω
	# 			0.0, # ω
	# 			0.0 # f
	# 		)

	orbElem = let 
        em = EllipsoidModel()
        elem = KeplerianElements(
				epoch |> datetime2julian, # t
				em.ellipsoid.a+h, # a
				0.0, # e
				i, # i
				0.0, # Ω
				0.0, # ω
				0.0 # f
			)
        elem
    end

	orbp = Propagators.init(Val(:J2), orbElem)

	acc = ground_facility_accesses(orbp, gs; minimum_elevation=mea, duration=24*3600, step=step, f_eci_to_ecef=_eci_to_ecef);

	gap = ground_facility_gaps(orbp, gs; minimum_elevation=mea, duration=24*3600, step=step, f_eci_to_ecef=_eci_to_ecef);

	acc,gap
end;

# @info """ 
# * TEME-PEF
# $(acc1)
# $(gap1)
# """
@info """
* J2000-ITRF
$(acc2)
$(gap2)
"""

# el = 0.1749 rad (10.02°) r = 1709.5km
#  az = 0.2867 rad (16.43°)
# julia> is_ground_facility_visible(r_sat,r_gs,deg2rad(10))false
# julia> is_ground_facility_visible(r_sat,r_gs,deg2rad(9.9))
# true
# julia> is_ground_facility_visible(r_sat,r_gs,deg2rad(9.95))true
# julia> is_ground_facility_visible(r_sat,r_gs,deg2rad(9.98))
# false
# julia> r_gs = locations[1].rv.ecefECEF Coordinate:
#  x = 5.528256639292835e6 y = 0.0
#  z = 3.170373735383637e6
# julia> r_sat = constellation[1].rv.pv.p_ecef
# ECEF Coordinate: x = 5.547817213164223e6
#  y = 1.6147583722314239e6 z = 3.7313600951122735e6

## TEST
r_gs = [5.528256639292835e6, 0.0, 3.170373735383637e6];
r_sat = [5.547817213164223e6, 1.6147583722314239e6, 3.7313600951122735e6];
em = EllipsoidModel();
gs = UserView(ECEF(r_gs...),em);
R_ecef_enu = _rotmat_ecef_to_enu(gs.lla.lat, gs.lla.lon);

@info "TelecomUtils ERA"
get_era(gs,ECEF(r_sat...))

@info "New Implementation"
rad2deg(SatelliteAnalysis._get_elevation(r_sat,r_gs,R_ecef_enu))
SatelliteAnalysis.is_ground_facility_visible(r_sat,r_gs,R_ecef_enu,deg2rad(9.97))
SatelliteAnalysis.is_ground_facility_visible(r_sat,r_gs,R_ecef_enu,deg2rad(9.98))

@info "Old Implementation"
rad2deg(SatelliteAnalysis._get_elevation_old(r_sat,r_gs))
SatelliteAnalysis.is_ground_facility_visible_old(r_sat,r_gs,deg2rad(9.97))
SatelliteAnalysis.is_ground_facility_visible_old(r_sat,r_gs,deg2rad(9.98))

## TEST2
@info "Equator"
sat1=SatView(LLA(0.0, 0.0, 1000e3));
gs1=UserView(LLA(0.0, 0.0, 0.0));
get_era(gs1, sat1)
rad2deg(SatelliteAnalysis._get_elevation_old(sat1.ecef, gs1.ecef))
R_ecef_enu1 = _rotmat_ecef_to_enu(gs1.lla.lat, gs1.lla.lon);
rad2deg(SatelliteAnalysis._get_elevation(sat1.ecef, gs1.ecef, R_ecef_enu1))

@info "North LAT"
sat2=SatView(LLA(deg2rad(42), deg2rad(11), 1000e3));
gs2=UserView(LLA(deg2rad(42), deg2rad(11), 0.0));
get_era(gs2, sat2)
rad2deg(SatelliteAnalysis._get_elevation_old(sat2.ecef, gs2.ecef))
R_ecef_enu2 = _rotmat_ecef_to_enu(gs2.lla.lat, gs2.lla.lon);
rad2deg(SatelliteAnalysis._get_elevation(sat2.ecef, gs2.ecef, R_ecef_enu2))

