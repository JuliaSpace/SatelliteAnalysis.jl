# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==============================================================================
#
#   Compute the beta angle of a satellite.
#
# References
# ==============================================================================
#
# [1] Mortari, D., Wilkins, M. P., and Bruccoleri, C.  On Sun-Synchronous Orbits
#     and Associated Constellations
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export beta_angle

"""
    beta_angle(jd₀::Number, a::Number, e::Number, i::Number, Ω::Number, Δt::Integer; kwargs...)

Compute the beta angle of an orbit with semi-major axis `a` [m], eccentricity
`e`, inclination `i` [rad], and initial right ascension of the ascending node
`Ω₀` [rad]. The orbit epoch, which is also the day in which the analysis will
begin, is `jd₀` [Julian Day]. The analysis will be performed for each day during
`Δt` days. The algorithm was obtained from **[1]**.

!!! note
    It is expected that the input elements are represented in the TOD reference
    frame. If it is not the case, they can be converted using the function
    `change_oe_frame` of **SatelliteToolbox.jl**.

# Keywords

- `pert::Symbol`: Select the perturbation terms that must be used when
    propagating the right ascencion of the ascending node. The possible values
    are:
    - `:J0`: Consider a Keplerian orbit.
    - `:J2`: Consider the perturbation terms up to J2.
    - `:J4`: Consider the perturbation terms J2, J4, and J2².
    (**Default**: `:J2`)

# Returns

- An array with the beta angle [rad] for each day.

# References

- **[1]**: Mortari, D., Wilkins, M. P., and Bruccoleri, C.  On Sun-Synchronous
           Orbits and Associated Constellations

# Extended help

The beta angle is defined as the angle between the orbit plane and the Sun. The
former is aligned with the orbit angular momentum vector.

The beta angle is useful when computing the mean amount of solar radiation a
satellite receives in a particular orbit.

## Examples

```julia
julia> using SatelliteAnalysis, UnicodePlots

julia> jd₀ = SatelliteAnalysis.date_to_jd(2021, 1, 1, 0, 0, 0)

julia> β = beta_angle(
    jd₀,
    7130.982e3,
    0.001111,
    deg2rad(98.405),
    SatelliteAnalysis.compute_RAAN_lt(jd₀, 22.5),
    365
)

julia> lineplot(0:1:364, rad2deg.(β), xlabel = "Day", ylabel = "Beta angle [°]")
                     ┌────────────────────────────────────────┐
                  30 │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
                     │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
                     │⠀⣠⠞⠉⠙⢦⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
                     │⡴⠁⠀⠀⠀⠀⠹⡄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
                     │⠀⠀⠀⠀⠀⠀⠀⠹⡄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣰⠃⠀⠀⠀│
                     │⠀⠀⠀⠀⠀⠀⠀⠀⠹⡄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡰⠃⠀⠀⠀⠀│
                     │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠹⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⣀⡀⠀⠀⠀⠀⠀⠀⢀⠞⠁⠀⠀⠀⠀⠀│
   Beta angle [°]    │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠱⡀⠀⠀⠀⠀⠀⠀⠀⠀⢀⠴⠋⠉⠀⠉⠙⠒⠤⠤⠤⠖⠉⠀⠀⠀⠀⠀⠀⠀│
                     │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠱⡄⠀⠀⠀⠀⠀⢀⠔⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
                     │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠘⠦⣀⣀⣠⠔⠉⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
                     │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
                     │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
                     │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
                     │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
                  10 │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
                     └────────────────────────────────────────┘
                     ⠀0⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀400⠀
                     ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀Day⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
```
"""
function beta_angle(
    jd₀::Number,
    a::Number,
    e::Number,
    i::Number,
    Ω₀::Number,
    Δt::Integer;
    pert::Symbol = :J2
)
    # Check the input.
    Δt ≤ 0 && error("The number of days for the analysis `Δt` must be equal or higher than 0")

    # Vector of the days in which the beta angle will be computed.
    days = 1:1:Δt

    # Output vector.
    β = Vector{Float64}(undef, Δt)

    # RAAN rotation rate [rad/day].
    δΩ = 86400 * dRAAN(a, e, i, pert)

    # Loop
    @inbounds for d in days
        # Compute the RAAN at the day d.
        Ω = Ω₀ + δΩ * d

        # Compute the unit vector aligned with the orbit normal `n` represented
        # in the TOD reference frame.
        D_tod_o = angle_to_dcm(-i, -Ω, 0, :XZX)
        n̄_tod   = D_tod_o * (@SVector [0, 0, 1])

        # Compute the Sun position at noon (UT) represented in the TOD reference
        # frame.
        s_mod     = sun_position_i(jd₀ + d)
        D_tod_mod = r_eci_to_eci(TOD(), jd₀, MOD(), jd₀)
        s̄_tod     = D_tod_mod * s_mod / norm(s_mod)

        # Compute the beta angle, which is the angle between the Sun vector and
        # the orbit plane.
        β[d] = abs(π / 2 - acos(n̄_tod ⋅ s̄_tod))
    end

    return β
end
