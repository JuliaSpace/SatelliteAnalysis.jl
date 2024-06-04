using SatelliteSimulationToolkit
using TelecomUtils
using SatelliteOrbitPropagation

@testset "Gap/Access Single Satellite - Benchamrk against SatelliteAnalysis.jl" begin
    ## Test against SatelliteAnalysis.jl
    ## //NOTE: A discrepancy between SatelliteAnalysis.jl and this test is expected due to the different underlying geometrical framework.
    # This leads to:
    # - A slighlty different results in terms of measured visibility angles GROUND-SAT. 
    # - As a consequence the start/end DateTime for accesses/gaps slighlty differ. 
    # - Furthermore, with our implementation we assume to have an error due to the simulation granularity: up to (2*mySim.timeStep+tol) for the access/gap duration, and up to (mySim.timeStep+tol) for the access/gap start/end (this do not influence the quality of the results for our use cases).
    # - It has been noticed that the mismatch between SatelliteSimulationToolkit.jl and SatelliteAnalysis.jl increases with the size of the timestep so high tolerance must be kept in the tests.
    
    
    eop_iau1980 = fetch_iers_eop(Val(:IAU1980));

    function _eci_to_ecef(r_i::AbstractVector, jd::Number)
        R = r_eci_to_ecef(J2000(), ITRF(), jd, eop_iau1980)
        # R = r_eci_to_ecef(TEME(), PEF(), jd)
        return R*r_i
    end

    Logging.disable_logging(Logging.Info)

    em = EllipsoidModel()

    sweepParam = (;
        h = collect(500e3:500e3:2000e3),
        i = deg2rad.(collect(20:20:80)),
        ts = collect(30:30:120),
        mea = deg2rad.(collect(10:10:60)),
        lla = map(x -> LLA(deg2rad(x),0.0,0.0), collect(0:30:90)),
        
        ## Err (J2000-ITRF)
        # h = 1000e3,
        # i = deg2rad.(40),
        # ts = 10,
        # mea = deg2rad.(10),
        # lla = [LLA(deg2rad(60),0.0,0.0)],
        ## Err (J2000-ITRF)
        # h = [500e3],
        # i = [deg2rad(20)],
        # ts = [120],
        # mea = [deg2rad(60)],
        # lla = [LLA(deg2rad(0),0.0,0.0)],
        ## Err (TEME-PEF)
        # h = [1500e3],
        # i = [deg2rad(60)],
        # ts = [120],
        # mea = [deg2rad(60)],
        # lla = [LLA(deg2rad(0),0.0,0.0)],
        ## Err (TEME-PEF)
        # h = [1500e3],
        # # i = [deg2rad(80)],
        # ts = [120],
        # mea = [deg2rad(0)],
        # lla = [LLA(deg2rad(0),0.0,0.0)],
    )
    for (h, i, ts, mea, lla) in Iterators.product(sweepParam.h, sweepParam.i, sweepParam.ts, sweepParam.mea, sweepParam.lla)
        err = Dates.Second(4*ts)
        
        mySim = SimParamInit(
            startDate = DateTime("2024-01-01"),
            duration = 24,
            timeStep = ts,
            propagator = :J2,
            rootFolder = ""
        )

        constellation = let
            satConstInit = ConstellationInit([
                SubConstellationInit(
                    h                   = h,    
                    T                   = 1,
                    P                   = 1,         
                    i                   = i,        
                    e                   = 0.0,
                    raan0               = 0.0,          
                    per                 = 0.0, 
                    geometry 		    = Auto(),
                    ID 				    = 1
                )
            ],
                    em = em
            )
            
            # Initialize the Constellation
            init = get_constellation(satConstInit; epoch=mySim.startDate)
            
            # Create the Satellites Vector
            sats = sat_from_init(init; propagator=mySim.propagator, em)
            
            sats
        end

        locations = [
            StaticLandLocation(
                UserView(lla, em),
                BasicUserClass(),
                mea,
                "none"		
            )
        ]

        args = (mySim, constellation, locations)
        kwargs = (; 
            ceilChunkSize = 1e6,
            geoArcSepAng = nothing,
            geoArcRes = nothing,
            pNGSOGridRes = 0.01,
            latVecRes = deg2rad(1),
            q = [0.05, 0.02, 0.0],
            np_quant = 500,
            np_ecdf = 100,
            saveFlag = false,
            owFlag = true,
            plotFlag = false,
            minAccessDuration = mySim.timeStep,
            statsFlag = (;
                el = true,
                doc = true,
                az = true,
                r = true,
                gapacc = true)
        )

        out = coverage_analysis(args...; kwargs...)

        accessesDF = map(locations) do thisLoc
            map(constellation) do thisSAT
                orbp = get_orbp(thisSAT)
                ground_facility_accesses(orbp, (thisLoc.rv.lla.lat, thisLoc.rv.lla.lon, thisLoc.rv.lla.alt); minimum_elevation=mea, duration=mySim.duration*3600, step=mySim.timeStep, f_eci_to_ecef=_eci_to_ecef)
            end
        end
        gapsDF = map(locations) do thisLoc
            map(constellation) do thisSAT
                orbp = get_orbp(thisSAT)
                ground_facility_gaps(orbp, (thisLoc.rv.lla.lat, thisLoc.rv.lla.lon, thisLoc.rv.lla.alt); minimum_elevation=mea, duration=mySim.duration*3600, step=mySim.timeStep, f_eci_to_ecef=_eci_to_ecef)
            end
        end

        if size(out.gapacc.dfAcc[1],1) > 0
            for k in 1:size(out.gapacc.dfAcc[1],1)
                # Measure the computation precision of the access start/end (it should be below mySim.timeStep+tol as explained, we assume 2*mySim.timeStep as max error here)
                # !(abs(out.gapacc.dfAcc[1].acc_start[k] - accessesDF[1][1].access_beginning[k]) < err) && error("ts:$(ts), h:$(h), i:$(rad2deg(i)), mea:$(rad2deg(mea)), lla:$(lla)")
                @test abs(out.gapacc.dfAcc[1].acc_start[k] - accessesDF[1][1].access_beginning[k]) < err
                # !(abs(out.gapacc.dfAcc[1].acc_end[k] - accessesDF[1][1].access_end[k]) < err) && error("ts:$(ts), h:$(h), i:$(rad2deg(i)), mea:$(rad2deg(mea)), lla:$(lla)")
                @test abs(out.gapacc.dfAcc[1].acc_end[k] - accessesDF[1][1].access_end[k]) < err
                # Measure the computation precision of the access duration (it should be below 2*mySim.timeStep+tol as explained, we assume 3*mySim.timeStep as max error here)
                # !(abs(out.gapacc.dfAcc[1].duration[k] - accessesDF[1][1].duration[k]) < err.value) && error("ts:$(ts), h:$(h), i:$(rad2deg(i)), mea:$(rad2deg(mea)), lla:$(lla)")
                @test abs(out.gapacc.dfAcc[1].duration[k] - accessesDF[1][1].duration[k]) < err.value
            end
        end
        
        if size(out.gapacc.dfGap[1],1) > 0
            for k in 1:size(out.gapacc.dfGap[1],1)
                # Measure the computation precision of the gap start/end (it should be below mySim.timeStep+tol as explained, we assume 2*mySim.timeStep as max error here)
                # !(abs(out.gapacc.dfGap[1].gap_start[k] - gapsDF[1][1].gap_beginning[k]) < err) && error("ts:$(ts), h:$(h), i:$(rad2deg(i)), mea:$(rad2deg(mea)), lla:$(lla)")
                @test abs(out.gapacc.dfGap[1].gap_start[k] - gapsDF[1][1].gap_beginning[k]) < err
                # !(abs(out.gapacc.dfGap[1].gap_end[k] - gapsDF[1][1].gap_end[k]) < err) && error("ts:$(ts), h:$(h), i:$(rad2deg(i)), mea:$(rad2deg(mea)), lla:$(lla)")
                @test abs(out.gapacc.dfGap[1].gap_end[k] - gapsDF[1][1].gap_end[k]) < err
                # Measure the computation precision of the gap duration (it should be below 2*mySim.timeStep+tol as explained, we assume 3*mySim.timeStep as max error here)
                # !(abs(out.gapacc.dfGap[1].duration[k] - gapsDF[1][1].duration[k]) < err.value) && error("ts:$(ts), h:$(h), i:$(rad2deg(i)), mea:$(rad2deg(mea)), lla:$(lla)")
                @test abs(out.gapacc.dfGap[1].duration[k] - gapsDF[1][1].duration[k]) < err.value
            end
        end
    end
end