Ground Facilities
=================

## Ground Facility Accesses

We can use the function:

```julia
ground_facility_accesses(orbp, [(WGS84)]; kwargs...) -> DataFrame
```

to compute the accesses of a satellite with orbit propagator `orbp` (see `Propagators.init`)
to the ground facilities defined in the vector `[(WGS84)]`. The analysis interval begins in
the propagator epoch plus `initial_time` and lasts for `duration` [s], where both are
keywords.

The ground facilities are specified using a vector of tuples with three numbers:

    Tuple{T1, T2, T3} where {T1 <: Number, T2 <: Number, T3 <: Number}

containing the WGS84 position of each ground facility `[(WGS84)]`:

    (latitude [rad], longitude [rad], altitude [m])

Those geodetic information are transformed to an ECEF vector using the function
`geodetic_to_ecef`.

The following keywords are available:

- `duration::Number`: Duration of the analysis [s].
    (**Default** = 86400)
- `f_eci_to_ecef::Function`: Function to convert the orbit propagator position represented
    in the Earth-centered inertial (ECI) reference frame to the Earth-centered, Earth-fixed
    (ECEF) reference frame. The signature must be

    ```julia
    f_eci_to_ecef(r_i::AbstractVector, jd::Number) -> AbstractVector
    ```

    and it must return the position vector `r_i` represented in the ECEF at the instant `jd`
    [Julian Day]. By default, we use TEME as the ECI and PEF as the ECEF.
    (**Default**: `_ground_facility_default_eci_to_ecef`)
- `initial_time::Number`: Initial time of the analysis after the propagator epoch [s].
    (**Default** = 0)
- `minimum_elevation::Number`: Minimum elevation angle for communication between the
    satellite and the ground facilities [rad].
    (**Default** = 10°)
- `reduction::Function`: A function that receives a boolean vector with the visibility
    between the satellite and each ground facility. It must return a boolean value
    indicating if the access must be computed or not. This is useful to merge access time
    between two or more facilities.
    (**Default** = `v -> |(v...)` *i.e.* compute the access if at least one ground
    facilities is visible)
- `step::Number`: The step [s] used to propagate the orbit. Notice that we perform a cross
    tuning to accurately obtain the access time. However, if an access is lower than the
    step, it can be neglected.
    (**Default** = 60)
- `unit::Symbol`: Select the unit in which the duration will be computed. The possible
    values are:
    - `:s` for seconds (**Default**);
    - `:m` for minutes; or
    - `:h` for hours.

This function returns a `DataFrame` with three columns:

- `access_beginning`: Time of the access beginning [UTC] encoded using `DateTime`.
- `access_end`: Time of the access end [UTC] encoded using `DateTime`.
- `duration`: Duration of the access [s].
  The unit of the column `duration` is stored in the `DataFrame` using metadata.

### Examples

Let's compute the access time of the Amazonia-1 satellite to the INPE's ground station at
Cuiabá, MT, Brazil.

First, we need to define an orbit propagator. In this case, we will use the TLE obtained
from CelesTrak when this documentation was written:

```jldoctest ground_facility_access
julia> tle_amz1 = tle"""
       AMAZONIA 1              
       1 47699U 21015A   24008.13079366  .00000299  00000+0  10693-3 0  9994
       2 47699  98.4120  87.2350 0001570  92.1147 268.0222 14.40836963150331
       """
TLE:
                      Name : AMAZONIA 1
          Satellite number : 47699
  International designator : 21015A
        Epoch (Year / Day) : 24 /   8.13079366 (2024-01-08T03:08:20.572)
        Element set number : 999
              Eccentricity :   0.00015700
               Inclination :  98.41200000 deg
                      RAAN :  87.23500000 deg
       Argument of perigee :  92.11470000 deg
              Mean anomaly : 268.02220000 deg
           Mean motion (n) :  14.40836963 revs / day
         Revolution number : 15033
                        B* :   0.00010693 1 / er
                     ṅ / 2 :     2.99e-06 rev / day²
                     n̈ / 6 :            0 rev / day³

julia> orbp = Propagators.init(Val(:SGP4), tle_amz1)
OrbitPropagatorSgp4{Float64, Float64}:
   Propagator name : SGP4 Orbit Propagator
  Propagator epoch : 2024-01-08T03:08:20.572
  Last propagation : 2024-01-08T03:08:20.572
```

Now, we can compute the access time during five days considering the INPE's station at
Cuiabá:

```jldoctest ground_facility_access
julia> ground_facility_accesses(
           orbp,
           [(-(15 + 33/60) |> deg2rad, -(56 + 04/60) |> deg2rad, 0)];
           duration = 5 * 86400,
           minimum_elevation = 5 |> deg2rad,
           unit = :m
       )
21×3 DataFrame
 Row │ access_beginning         access_end               duration  
     │ DateTime                 DateTime                 Float64   
─────┼─────────────────────────────────────────────────────────────
   1 │ 2024-01-08T03:08:20.572  2024-01-08T03:09:17.511   0.948983
   2 │ 2024-01-08T13:56:46.830  2024-01-08T14:08:49.523  12.0449
   3 │ 2024-01-08T15:39:26.373  2024-01-08T15:43:46.376   4.33338
   4 │ 2024-01-09T02:18:04.189  2024-01-09T02:30:10.279  12.1015
   5 │ 2024-01-09T04:00:17.780  2024-01-09T04:04:24.011   4.10385
   6 │ 2024-01-09T13:18:23.251  2024-01-09T13:28:33.092  10.164
   7 │ 2024-01-09T14:56:43.963  2024-01-09T15:07:18.440  10.5746
   8 │ 2024-01-10T01:39:24.600  2024-01-10T01:49:48.246  10.3941
   9 │ 2024-01-10T03:17:59.993  2024-01-10T03:28:20.957  10.3494
  10 │ 2024-01-10T12:42:43.882  2024-01-10T12:45:52.993   3.15185
  11 │ 2024-01-10T14:16:26.860  2024-01-10T14:28:34.516  12.1276
  12 │ 2024-01-11T01:03:20.129  2024-01-11T01:06:41.796   3.36112
  13 │ 2024-01-11T02:37:47.498  2024-01-11T02:49:52.306  12.0801
  14 │ 2024-01-11T13:37:23.787  2024-01-11T13:48:49.336  11.4258
  15 │ 2024-01-11T15:17:31.189  2024-01-11T15:26:04.343   8.55257
  16 │ 2024-01-12T01:58:34.688  2024-01-12T02:10:09.158  11.5745
  17 │ 2024-01-12T03:38:38.236  2024-01-12T03:46:54.916   8.278
  18 │ 2024-01-12T12:59:52.458  2024-01-12T13:07:49.894   7.95727
  19 │ 2024-01-12T14:36:24.532  2024-01-12T14:48:04.736  11.6701
  20 │ 2024-01-13T01:20:41.862  2024-01-13T01:28:55.139   8.22128
  21 │ 2024-01-13T02:57:44.892  2024-01-13T03:08:20.572  10.5947
```
