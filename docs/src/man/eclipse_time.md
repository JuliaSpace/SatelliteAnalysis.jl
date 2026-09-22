# Eclipse Time

```@meta
CurrentModule = SatelliteAnalysis
```

```@setup eclipse_time
using SatelliteAnalysis
```

The eclipse time is the period the satellite does not receive sunlight due to the Earth
shadow. This information is paramount for mission design since it directly interferes in the
power and thermal subsystems.

We can compute the eclipse time of a satellite using the function
`eclipse_time_summary`:

```@docs; canonical = false
eclipse_time_summary
```

!!! note

    If we want to verify the current lighting condition in a satellite (sunlight, umbra, or
    penumbra), see the function [`lighting_condition`](@ref).

## Examples

We will compute the eclipse time of the Amazonia-1 mission for one year. The first thing we
need to do is define the orbit:

```@repl eclipse_time
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

The next step is to define the desired propagator:

```@repl eclipse_time
orbp = Propagators.init(Val(:J2), orb)
```

Now, we can use the function `eclipse_time_summary` to obtain the eclipse time information
for each day of the year:

```@repl eclipse_time
df = eclipse_time_summary(orbp; time_unit = :min)
```

Finally, we can use the `DataFrame` to analyze the result. For example, the maximum eclipse
time in an orbit is:

```@repl eclipse_time
maximum(df.penumbra .+ df.umbra)
```

_i.e._, 34.66 minutes.

## References

- **[1]** **Longo, C. R. O., Rickman, S. L (1995)**. _Method for the Calculation of
  Spacecraft Umbra and Penumbra Shadow Terminator Points_. **NASA Technical Paper** 3547.
