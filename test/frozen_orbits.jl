## Description #############################################################################
#
# Tests related to the frozen orbit computation.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# == File: ./src/frozen_orbits.jl ==========================================================

# -- Function: frozen_orbit ----------------------------------------------------------------

# Gravity model that wraps another model with fully normalized coefficients, providing its
# zonal coefficients without normalization. It reports the normalization `coefficient_norm`,
# allowing us to test all the conversions in `frozen_orbit`.
struct FrozenOrbitTestGravityModel{T, M <: AbstractGravityModel{T}, N} <:
       AbstractGravityModel{T}
    model::M
    coefficient_norm::N
end

function FrozenOrbitTestGravityModel(model::AbstractGravityModel{T}, norm) where {T}
    return FrozenOrbitTestGravityModel{T, typeof(model), typeof(norm)}(model, norm)
end

function GravityModels.coefficients(
    gm::FrozenOrbitTestGravityModel, degree::Int, order::Int, time::Number
)
    clm, slm = GravityModels.coefficients(gm.model, degree, order, time)
    return clm * √(2degree + 1), slm * √(2degree + 1)
end

GravityModels.coefficient_norm(gm::FrozenOrbitTestGravityModel) = gm.coefficient_norm
GravityModels.radius(gm::FrozenOrbitTestGravityModel) = GravityModels.radius(gm.model)

function GravityModels.maximum_degree(gm::FrozenOrbitTestGravityModel)
    return GravityModels.maximum_degree(gm.model)
end

@testset "Function frozen_orbit" begin
    # == Default ===========================================================================

    e, ω = frozen_orbit(7130.982e3, 98.410 |> deg2rad)
    @test e ≈ 0.0011641853028456078 atol = 1e-20
    @test ω == π / 2

    # A non-positive degree selects the model's available maximum degree, rather than the
    # default truncation at degree 53.  This also guards against the undefined `grav_model`
    # variable.
    e_max, ω_max = frozen_orbit(7130.982e3, 98.410 |> deg2rad; max_degree = 0)
    @test e_max ≈ 0.010085085903178675 atol = 1e-20
    @test ω_max == π / 2

    # == Lower maximum degree ==============================================================

    e, ω = frozen_orbit(7130.982e3, 98.410 |> deg2rad; max_degree = 5)
    @test e ≈ 0.0011108978494835141 atol = 1e-20
    @test ω == π / 2

    # == Comparing with analytical solution ================================================

    # Let's compute the frozen orbit analytically considering only the third degree.
    J₂   = 0.0010826266835531513
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

    # == Coefficient Normalizations ========================================================

    # The result must not depend on how the gravity model normalizes its coefficients.
    # Notice that the Schmidt quasi-normalization does not modify the zonal terms.
    for coefficient_norm in (Val(:unnormalized), Val(:schmidt))
        gm = FrozenOrbitTestGravityModel(jgm3, coefficient_norm)
        e_n, ω_n = frozen_orbit(7130.982e3, 98.410 |> deg2rad; gravity_model = gm)

        @test e_n ≈ e rtol = 1e-13
        @test ω_n == ω
    end

    gm = FrozenOrbitTestGravityModel(jgm3, Val(:unknown))
    @test_throws ArgumentError frozen_orbit(
        7130.982e3, 98.410 |> deg2rad; gravity_model = gm
    )

    # == Equatorial Orbits =================================================================

    # The frozen orbit is not defined for equatorial orbits.
    @test_throws ArgumentError frozen_orbit(7130.982e3, 0; gravity_model = jgm3)
    @test_throws ArgumentError frozen_orbit(7130.982e3, π; gravity_model = jgm3)
    @test_throws ArgumentError frozen_orbit(7130.982e3, -0.1; gravity_model = jgm3)
    @test_throws ArgumentError frozen_orbit(7130.982e3, 3.2; gravity_model = jgm3)
end
