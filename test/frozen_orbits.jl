## Description #############################################################################
#
# Tests related to the frozen orbit computation.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# == File: ./src/frozen_orbits.jl ==========================================================

# -- Function: frozen_orbit ----------------------------------------------------------------

@testset "Function frozen_orbit" begin
    # == Default ===========================================================================

    e, ω = frozen_orbit(7130.982e3, 98.410 |> deg2rad)
    @test e ≈ 0.0011641853028456078 atol = 1e-20
    @test ω == π / 2

    # == Lower maximum degree ==============================================================

    e, ω = frozen_orbit(7130.982e3, 98.410 |> deg2rad; max_degree = 5)
    @test e ≈ 0.0011108978494835141 atol = 1e-20
    @test ω == π / 2

    # == Comparing with analytical solution ================================================

    # Let's compute the frozen orbit analytically considering only the third degree.
    J₂   =  0.0010826266835531513
    J₃   = -2.5326564853322355e-6
    e_d3 = - (6.3781363e6 / 7130.982e3) * J₃ / J₂ * sind(98.410) / 2

    e, ω = frozen_orbit(7130.982e3, 98.410 |> deg2rad; max_degree = 3)

    @test e ≈ e_d3 atol = 1e-20
    @test ω == π / 2

    # == Using lower degree than allowed ===================================================

    e, ω = frozen_orbit(7130.982e3, 98.410 |> deg2rad; max_degree = 1)

    @test e ≈ e_d3 atol = 1e-20
    @test ω == π / 2

    e, ω = frozen_orbit(7130.982e3, 98.410 |> deg2rad; max_degree = 2)

    @test e ≈ e_d3 atol = 1e-20
    @test ω == π / 2

    # == Argument of perigee at 270° =======================================================

    e, ω = frozen_orbit(7130.982e3, 64 |> deg2rad)
    @test e ≈ 0.004201259311864032 atol = 1e-20
    @test ω ≈ 4.71238898038469

    # == Using another gravity model =======================================================

    jgm3 = GravityModels.load(IcgemFile, fetch_icgem_file(:JGM3))
    e, ω = frozen_orbit(7130.982e3, 98.410 |> deg2rad; gravity_model = jgm3)
    @test e ≈ 0.001163484769069545 atol = 1e-20
    @test ω == π / 2
end

