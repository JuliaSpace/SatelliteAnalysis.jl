using Documenter
using SatelliteAnalysis

makedocs(
    modules = [SatelliteAnalysis],
    format = Documenter.HTML(
        prettyurls = !("local" in ARGS),
        canonical = "https://juliaspace.github.io/SatelliteAnalysis.jl/stable/",
    ),
    sitename = "Satellite Analysis",
    authors = "Ronan Arraes Jardim Chagas",
    pages = [
        "Home" => "index.md",
        "Beta Angle" => "man/beta_angle.md",
        "Eclipse Time" => "man/eclipse_time.md",
        "Frozen Orbits" => "man/frozen_orbits.md",
        "Ground Facilities" => [
            "Ground Facility Accesses" => "man/ground_facilities/ground_facility_accesses.md",
            "Ground Facility Gaps" => "man/ground_facilities/ground_facility_gaps.md",
            "Ground Facility Visibility Circle" => "man/ground_facilities/ground_facility_visibility_circle.md"
        ],
        "Ground Track" => "man/ground_track.md",
        "Orbits" => [
            "Sun Synchronous Orbits" => "man/orbits/sun_synchronous_orbits.md",
        ],
        "Library" => "lib/library.md",
    ],
)

deploydocs(
    repo = "github.com/JuliaSpace/SatelliteAnalysis.jl.git",
    target = "build",
)
