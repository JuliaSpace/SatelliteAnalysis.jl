# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==============================================================================
#
#   Tests related to the beta angle analysis.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# File: ./src/beta_angle.jl
# ==============================================================================

# Function: beta_angle
# ------------------------------------------------------------------------------

@testset "Function beta_angle" begin
  jd₀ = SatelliteAnalysis.date_to_jd(2021, 1, 1, 0, 0, 0)

  β = beta_angle(
      jd₀,
      7130.982e3,
      0.001111,
      deg2rad(98.405),
      SatelliteAnalysis.compute_RAAN_lt(jd₀, 22.5),
      5
  )

  β_expected = [
      0.43291445352510594
      0.43481148866086317
      0.43668726746800246
      0.4385397145910499
      0.44036676126144947
  ]

  @test β ≈ β_expected

  β = beta_angle(
      jd₀,
      7130.982e3,
      0.001111,
      deg2rad(98.405),
      SatelliteAnalysis.compute_RAAN_lt(jd₀, 22.5) + π,
      5
  )

  β_expected = [
      -0.31068187881262466
      -0.31295402147649054
      -0.31524372027870284
      -0.3175489312409039
      -0.3198676044177273
  ]

  @test β ≈ β_expected
end

@testset "Function beta_angle (errors)" begin
  jd₀ = SatelliteAnalysis.date_to_jd(2021, 1, 1, 0, 0, 0)

  @test_throws Exception beta_angle(
      jd₀,
      7130.982e3,
      0.001111,
      deg2rad(98.405),
      SatelliteAnalysis.compute_RAAN_lt(jd₀, 22.5),
      0
  )

  @test_throws Exception beta_angle(
      jd₀,
      7130.982e3,
      0.001111,
      deg2rad(98.405),
      SatelliteAnalysis.compute_RAAN_lt(jd₀, 22.5),
      -5
  )
end
