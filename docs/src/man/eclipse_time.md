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

We can compute the eclipse time of a satellite using the function:

```julia
eclipse_time_summary(orbp::OrbitPropagator; kwargs...) -> DataFrame
```

This function computes the eclipse time summary for the orbit propagator `orbp`. The summary
is computed as the total time the object stays in the sunlight, penumbra, and umbra regions
per orbit at each day. The algorithm was adapted from **[1]**.

The following keywords are available:

- `num_days::Number`: Number of days in which the analysis will be performed.
  (**Default** = 365)
- `step::Number`: The step in which the propagation will occur. Notice that this function
  has a crossing estimation to accurately estimate the transition between the regions.
  However, if this step is very large, we may miss some small regions. If it is negative, it
  will be selected as the time in which the mean anomaly advances 0.5°.
  (**Default** = -1)
- `unit::Symbol`: Select the unit in which the results will be generated. The possible
  values are:
    - `:s` for seconds (**Default**);
    - `:m` for minutes; or
    - `:h` for hours.

The function returns a `DataFrame` with three columns:

- `sunlight`: Total sunlight time per orbit at each day [`unit`].
- `penumbra`: Total penumbra time per orbit at each day [`unit`].
- `umbra`: Total umbra time per orbit at each day [`unit`].

The unit of each column is stored in the `DataFrame` using metadata.

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
df = eclipse_time_summary(orbp; unit = :m)
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
