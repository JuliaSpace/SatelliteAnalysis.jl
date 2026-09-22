# Ground Facility Accesses

```@meta
CurrentModule = SatelliteAnalysis
```

```@setup ground_facility_access
using SatelliteAnalysis
```

We can compute the accesses of a satellite to a set of ground facilities using the
function `ground_facility_accesses`:

```@docs; canonical = false
ground_facility_accesses
```

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
    time_unit = :min
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
    time_unit = :min
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
    time_unit = :min
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
    time_unit = :min
)
```
