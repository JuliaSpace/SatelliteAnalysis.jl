## Description #############################################################################
#
# Compute the beta angle of a satellite.
#
## References ##############################################################################
#
# [1] Mortari, D., Wilkins, M. P., and Bruccoleri, C.  On Sun-Synchronous Orbits and
#     Associated Constellations
#
############################################################################################

export beta_angle

"""
    beta_angle(orb::KerplerianElements{Tepoch, T}, Δjd::Number; kwargs...) -> Float64

Compute the beta angle [rad] for the orbit `orb` after `Δjd` days from its epoch.

The algorithm was obtained from **[1]**.

!!! note
    It is expected that the input elements are represented in the TOD reference frame. If it
    is not the case, they can be converted using the function `orb_eci_to_eci` of
    **SatelliteToolboxTransformations.jl**.

# Keywords

- `perturbation::Symbol`: Select the perturbation terms that must be used when propagating
  the right ascencion of the ascending node. The possible values are:
    - `:J0`: Consider a Keplerian orbit.
    - `:J2`: Consider the perturbation terms up to J₂.
    - `:J4`: Consider the perturbation terms J₂, J₂², and J₄.
    (**Default**: `:J2`)

# References

- **[1]**: Mortari, D., Wilkins, M. P., and Bruccoleri, C.  On Sun-Synchronous Orbits and
    Associated Constellations

# Extended Help

The beta angle is the angle between the orbit plane and the Sun. The positive direction is
defined as that of the orbit angular momentum.

The beta angle is useful when computing the mean amount of solar radiation a satellite
receives in a particular orbit.

## Examples

```julia
julia> using SatelliteAnalysis, UnicodePlots

julia> jd₀ = date_to_jd(2021, 1, 1, 0, 0, 0)

julia> orb = KeplerianElements(
            jd₀,
            7130.982e3,
            0.001111,
            98.405 |> deg2rad,
            ltdn_to_raan(10.5, jd₀),
            90     |> deg2rad,
            0
        )

julia> β = beta_angle.(orb, 0:1:364)

julia> lineplot(0:1:364, rad2deg.(β), xlabel = "Day", ylabel = "Beta angle [°]")
                     ┌────────────────────────────────────────┐
                  30 │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
                     │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
                     │⠀⣠⠞⠉⠙⢦⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
                     │⡴⠁⠀⠀⠀⠀⠹⡄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
                     │⠀⠀⠀⠀⠀⠀⠀⠘⡄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡰⠃⠀⠀⠀│
                     │⠀⠀⠀⠀⠀⠀⠀⠀⠹⡄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⡼⠁⠀⠀⠀⠀│
                     │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠹⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣀⣀⣀⠀⠀⠀⠀⠀⠀⣠⠞⠀⠀⠀⠀⠀⠀│
   Beta angle [°]    │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠱⡀⠀⠀⠀⠀⠀⠀⠀⠀⣠⠖⠋⠀⠀⠀⠉⠓⠲⠤⠴⠚⠁⠀⠀⠀⠀⠀⠀⠀│
                     │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠹⣄⠀⠀⠀⠀⠀⣠⠞⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
                     │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠦⣄⣀⡠⠞⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
                     │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
                     │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
                     │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
                     │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
                  10 │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
                     └────────────────────────────────────────┘
                     ⠀0⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀400⠀
                     ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀Day⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
```
"""
function beta_angle(orb::KeplerianElements, Δjd::Number; perturbation::Symbol = :J2)
    # RAAN rotation rate [rad / day].
    δΩ = 86400 * raan_time_derivative(orb; perturbation = perturbation)

    # Compute the RAAN at the day d.
    Ω = orb.Ω + δΩ * Δjd

    # Obtain the epoch to compute the beta angle.
    jd = orb.t + Δjd

    # Compute the unit vector aligned with the orbit normal `n` represented in the TOD
    # reference frame.
    D_tod_o = angle_to_dcm(-orb.i, -Ω, 0, :XZX)
    n̄_tod   = D_tod_o * @SVector([0, 0, 1])

    # Compute the Sun position at noon (UT) represented in the TOD reference frame.
    s_mod     = sun_position_mod(jd)
    D_tod_mod = r_eci_to_eci(MOD(), jd, TOD(), jd)
    s̄_tod     = D_tod_mod * s_mod / norm(s_mod)

    # Compute the beta angle, which is the angle between the Sun vector and the orbit plane.
    β = π / 2 - acos(n̄_tod ⋅ s̄_tod)

    return β
end
