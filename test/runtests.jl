using Test

using LinearAlgebra
using StaticArrays
using SatelliteAnalysis

@testset "Beta Angle" verbose = true begin
    include("./beta_angle.jl")
end

@testset "Eclipse Time" verbose = true begin
    include("./eclipse_time.jl")
end

@testset "Lighting Condition" verbose = true begin
    include("./lighting_condition.jl")
end

# @testset "Miscellaneous" verbose = true begin
#     include("./misc.jl")
# end

@testset "Orbits" verbose = true begin
    include("./sun_syncrhonous_orbits.jl")
end
