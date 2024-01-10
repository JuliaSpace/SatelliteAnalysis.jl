Ground Facility Accesses
========================

```@meta
CurrentModule = SatelliteAnalysis
```

```@setup ground_facility_access
using SatelliteAnalysis
```

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

    ```
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

## Examples

Let's compute the access time of the Amazonia-1 satellite to the INPE's ground station at
Cuiabá, MT, Brazil.

First, we need to define an orbit propagator. In this case, we will use the TLE obtained
from CelesTrak when this documentation was written:

```@repl ground_facility_access
tle_amz1 = tle"""
    AMAZONIA 1
    1 47699U 21015A   24008.13079366  .00000299  00000+0  10693-3 0  9994
    2 47699  98.4120  87.2350 0001570  92.1147 268.0222 14.40836963150331
    """

orbp = Propagators.init(Val(:SGP4), tle_amz1)
```

Now, we can compute the access time during one day considering the INPE's station at Cuiabá:

```@repl ground_facility_access
ground_facility_accesses(
    orbp,
    [(-(15 + 33 / 60) |> deg2rad, -(56 + 04 / 60) |> deg2rad, 0)];
    duration = 1 * 86400,
    minimum_elevation = 5 |> deg2rad,
    unit = :m
)
```

If we want to change the reference frames used in the analysis, we must provide a function
`f_eci_to_ecef(r_i, jd)` that converts the vector `r_i` to the desired ECEF frame at the
instant `jd` [Julian Day]. For example, let's use the more precise ITRF instead of the PEF
(default):

```@repl ground_facility_access
const eop = fetch_iers_eop(Val(:IAU1980));

function f_eci_to_ecef(r_teme, jd)
    D_itrf_teme = r_eci_to_ecef(TEME(), ITRF(), jd, eop)
    r_itrf = D_itrf_teme * r_teme
    return r_itrf
end

ground_facility_accesses(
    orbp,
    [(-(15 + 33 / 60) |> deg2rad, -(56 + 04 / 60) |> deg2rad, 0)];
    duration = 1 * 86400,
    f_eci_to_ecef = f_eci_to_ecef,
    minimum_elevation = 5 |> deg2rad,
    unit = :m
)
```

We can also perform analyses using multiple ground facilities. For example, let's find the
accumulated access if we consider the INPE's stations at Cuiabá, MT, Brazil, and Alcântara,
MA, Brazil:

```@repl ground_facility_access
ground_facility_accesses(
    orbp,
    [
        (-(15 + 33 / 60) |> deg2rad, -(56 + 04 / 60) |> deg2rad, 0)
        (-( 2 + 20 / 60) |> deg2rad, -(44 + 24 / 60) |> deg2rad, 0)
    ];
    duration = 1 * 86400,
    minimum_elevation = 5 |> deg2rad,
    unit = :m
)
```

By default, the algorithm computes the accumulated access, _i.e._, it considers the access
active if either station has visibility to the satellite. We can change this logic by
overloading the function in the keyword parameter `reduction`. For example, let's compute
only the accesses when the satellite has visibility to both ground stations at the same
time:

```@repl ground_facility_access
ground_facility_accesses(
    orbp,
    [
        (-(15 + 33 / 60) |> deg2rad, -(56 + 04 / 60) |> deg2rad, 0)
        (-( 2 + 20 / 60) |> deg2rad, -(44 + 24 / 60) |> deg2rad, 0)
    ];
    duration = 1 * 86400,
    minimum_elevation = 5 |> deg2rad,
    reduction = v -> (&)(v...),
    unit = :m
)
```
