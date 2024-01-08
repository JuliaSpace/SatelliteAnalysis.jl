Eclipse Time
============

```@meta
CurrentModule = SatelliteAnalysis
DocTestSetup = quote
    using SatelliteAnalysis
end
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

## Examples

We will compute the eclipse time of the Amazonia-1 mission for one year. The first thing we
need to do is define the orbit:

```jldoctest eclipse_time
julia> jd₀ = date_to_jd(2021, 1, 1)
2.4592155e6

julia> orb = KeplerianElements(
           jd₀,
           7130.982e3,
           0.001111,
           98.405 |> deg2rad,
           ltdn_to_raan(10.5, jd₀),
           π / 2,
           0
       )
KeplerianElements{Float64, Float64}:
           Epoch :    2.45922e6 (2021-01-01T00:00:00)
 Semi-major axis : 7130.98     km
    Eccentricity :    0.001111
     Inclination :   98.405    °
            RAAN :   78.4021   °
 Arg. of Perigee :   90.0      °
    True Anomaly :    0.0      °
```

The next step is to define the desired propagator:

```jldoctest eclipse_time
julia> orbp = Propagators.init(Val(:J2), orb)
OrbitPropagatorJ2{Float64, Float64}:
   Propagator name : J2 Orbit Propagator
  Propagator epoch : 2021-01-01T00:00:00
  Last propagation : 2021-01-01T00:00:00
```

Now, we can use the function `eclipse_time_summary` to obtain the eclipse time information
for each day of the year:

```jldoctest eclipse_time
julia> df = eclipse_time_summary(orbp; unit = :m)
365×4 DataFrame
 Row │ date        sunlight  penumbra  umbra
     │ Date        Float64   Float64   Float64
─────┼─────────────────────────────────────────
   1 │ 2021-01-01   66.2105  0.340195  33.4493
   2 │ 2021-01-02   66.2308  0.340627  33.4285
   3 │ 2021-01-03   66.2461  0.340958  33.4129
   4 │ 2021-01-04   66.2623  0.341263  33.3964
   5 │ 2021-01-05   66.2824  0.341704  33.3759
   6 │ 2021-01-06   66.2932  0.341899  33.3649
   7 │ 2021-01-07   66.3129  0.342331  33.3447
   8 │ 2021-01-08   66.3274  0.342628  33.33
  ⋮  │     ⋮          ⋮         ⋮         ⋮
 359 │ 2021-12-25   66.1274  0.338075  33.5345
 360 │ 2021-12-26   66.1477  0.338525  33.5137
 361 │ 2021-12-27   66.1595  0.338762  33.5017
 362 │ 2021-12-28   66.18    0.339203  33.4808
 363 │ 2021-12-29   66.1956  0.33955   33.4648
 364 │ 2021-12-30   66.2122  0.339889  33.4479
 365 │ 2021-12-31   66.2327  0.340339  33.4269
                               350 rows omitted
```

Finally, we can use the `DataFrame` to analyze the result. For example, the maximum eclipse
time in an orbit is:

```jldoctest eclipse_time
julia> maximum(df.penumbra .+ df.umbra)
34.66395872764383
```

_i.e._, 34.66 minutes.

## References

- **[1]** Longo, C. R. O., Rickman, S. L (1995). Method for the Calculation of Spacecraft
  Umbra and Penumbra Shadow Terminator Points. NASA Technical Paper 3547.
