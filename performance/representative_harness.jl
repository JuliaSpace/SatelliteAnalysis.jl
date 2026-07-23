#!/usr/bin/env julia

# Dependency-free performance harness; run manually, not as part of tests.
using SatelliteAnalysis
using SatelliteToolbox: Propagators

"""Create the representative J2 orbit propagator used by the harness."""
function sample_propagator()
    jd = date_to_jd(2021, 1, 1)
    orbit = KeplerianElements(
        jd, 7130.982e3, 0.001111, deg2rad(98.405), ltdn_to_raan(10.5, jd), deg2rad(90), 0
    )
    return Propagators.init(Val(:J2), orbit)
end

"""Warm up and measure one invocation of `f`, reporting its time and allocation."""
function measure(name, f)
    f()
    GC.gc()
    measurement = @timed f()
    println(
        rpad(name, 38),
        " time=$(round(measurement.time; digits = 6)) s allocated=$(measurement.bytes) bytes",
    )
    return measurement.value
end

orbp = sample_propagator()
measure("ground_track", () -> ground_track(orbp; step = 300, duration = 3600))
measure(
    "ground_facility_accesses",
    () -> ground_facility_accesses(
        orbp, (0, 0, 0); duration = 3600, step = 300, num_chunks = 1
    ),
)
measure("eclipse_time_summary", () -> eclipse_time_summary(orbp; num_days = 5))
measure(
    "sun-synchronous solver", () -> sun_sync_orbit_from_angular_velocity(14 * 2π / 86400)
)
