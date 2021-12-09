using Test

using StaticArrays
using SatelliteAnalysis
import SatelliteToolbox: date_to_jd, sun_position_i

@testset "Beta angle" verbose = true begin
    include("./beta_angle.jl")
end

@testset "Lighting condition" verbose = true begin
    include("./lighting_condition.jl")
end

@testset "Miscellaneous" verbose = true begin
    include("./misc.jl")
end
