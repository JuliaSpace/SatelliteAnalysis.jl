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