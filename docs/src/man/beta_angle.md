Beta Angle
==========

```@meta
CurrentModule = SatelliteAnalysis
```

```@setup beta_angle
using SatelliteAnalysis
```

The beta angle is the angle between the orbit plane and the Sun, as shown in the following
figure. The positive direction is defined as that of the orbit angular momentum.

The beta angle is useful when computing the mean amount of solar radiation a satellite
receives in a particular orbit.

```@raw html
<div align="center">
  <img src="../../assets/beta_angle.png" alt="Beta Angle" width="50%"/>
</div>
```

We can compute the beta angle of an orbit using the function:

```julia
beta_angle(orb::KerplerianElements{Tepoch, T}, Δjd::Number; kwargs...) -> Float64
```

This function computes the beta angle [rad] for the orbit `orb` after `Δjd` days from its
epoch.

The algorithm was obtained from **[1]**.

!!! note
    It is expected that the input elements are represented in the TOD reference frame. If it
    is not the case, they can be converted using the function `orb_eci_to_eci` of
    **SatelliteToolboxTransformations.jl**.

The following keywords are available:

- `perturbation::Symbol`: Select the perturbation terms that must be used when propagating
  the right ascencion of the ascending node. The possible values are:
    - `:J0`: Consider a Keplerian orbit.
    - `:J2`: Consider the perturbation terms up to J₂.
    - `:J4`: Consider the perturbation terms J₂, J₂², and J₄.
    (**Default**: `:J2`)

## Examples

We will compute the beta angle of the Amazonia-1 mission for one year. The first thing we
need to do is define the orbit:

```@repl beta_angle
jd₀ = date_to_jd(2021, 1, 1)

orb = KeplerianElements(
    jd₀,
    7130.982e3,
    0.001111,
    98.405 |> deg2rad,
    ltdn_to_raan(10.5, jd₀),
    π / 2,
    0
)
```

Now, we can use the function `beta_angle` to obtain the beta angle [rad] for each day of the
year:

```@repl beta_angle
β = beta_angle.(orb, 1:365)
```

If we use **CairoMakie.jl** to plot, we obtain:

```@raw html
<div align="center">
  <img src="../../assets/beta_angle_amz1.png" alt="Beta Angle of the Amazonia-1 Mission" width="100%"/>
</div>
```

## References

- **[1]**: **Mortari, D., Wilkins, M. P., and Bruccoleri, C**. _On Sun-Synchronous Orbits and
  Associated Constellations_.
