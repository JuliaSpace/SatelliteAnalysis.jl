using Test

using SatelliteAnalysis

@testset "Beta angle" verbose = true begin
    include("./beta_angle.jl")
end

@testset "Miscellaneous" verbose = true begin
    include("./misc.jl")
end
