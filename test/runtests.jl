using Test

using LinearAlgebra
using Logging
using StaticArrays
using SatelliteAnalysis

@testset "Beta Angle" verbose = true begin
    include("./beta_angle.jl")
end

@testset "Eclipse Time" verbose = true begin
    include("./eclipse_time.jl")
end

@testset "Frozen Orbits" verbose = true begin
    include("./frozen_orbits.jl")
end

@testset "Ground Facilities" verbose = true begin
    include("./ground_facility_accesses_and_gaps.jl")
    include("./ground_facility_visibility_circle.jl")
end

@testset "Ground Tracks" verbose = true begin
    include("./ground_track.jl")
end

@testset "Lighting Condition" verbose = true begin
    include("./lighting_condition.jl")
end

@testset "Miscellaneous" verbose = true begin
   include("./misc.jl")
end

@testset "Plotting" verbose = true begin
    include("./plotting.jl")
end

@testset "Orbits" verbose = true begin
    include("./sun_syncrhonous_orbits.jl")
end
