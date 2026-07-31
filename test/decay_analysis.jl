## Description #############################################################################
#
# Tests related to the orbital decay analysis.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# == File: ./src/decay_analysis.jl =========================================================

# -- Function: decay_analysis --------------------------------------------------------------

@testset "Extension Loading" begin
    ext = Base.get_extension(SatelliteAnalysis, :SatelliteAnalysisDecayExt)
    @test !isnothing(ext)
end
