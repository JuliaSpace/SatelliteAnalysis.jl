Ground Facility Gaps
====================

```@meta
CurrentModule = SatelliteAnalysis
```

```@setup ground_facility_gap
using SatelliteAnalysis
```

We can use the function:

```julia
ground_facility_gaps(orbp, args...; duration::Number = 86400, initial_time::Number = 0, kwargs...) -> DataFrame
```

to compute the gaps between the accesses of ground facilities. The arguments and keywords
are the same as the ones used in the function [`ground_facility_accesses`](@ref) (see
[Ground Facility Accesses](@ref)).

Notice that the gap analysis starts in the orbit propagator epoch plus `initial_time` and
lasts for `duration` [s].

This function returns a `DataFrame` with three columns:

- `access_beginning`: Time of the access beginning [UTC] encoded using `DateTime`.
- `access_end`: Time of the access end [UTC] encoded using `DateTime`.
- `duration`: Duration of the access [s].
  The unit of the column `duration` is stored in the `DataFrame` using metadata.
## Examples

Let's compute the gaps of the Amazonia-1 satellite to the INPE's ground station at Cuiabá,
MT, Brazil.

First, we need to define an orbit propagator. In this case, we will use the TLE obtained
from CelesTrak when this documentation was written:

```@repl ground_facility_gap
tle_amz1 = tle"""
    AMAZONIA 1
    1 47699U 21015A   24008.13079366  .00000299  00000+0  10693-3 0  9994
    2 47699  98.4120  87.2350 0001570  92.1147 268.0222 14.40836963150331
    """

orbp = Propagators.init(Val(:SGP4), tle_amz1)
```

Now, we can compute the gaps during one day considering the INPE's station at Cuiabá:

```@repl ground_facility_gap
ground_facility_gaps(
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

```@repl ground_facility_gap
const eop = fetch_iers_eop(Val(:IAU1980));

function f_eci_to_ecef(r_teme, jd)
    D_itrf_teme = r_eci_to_ecef(TEME(), ITRF(), jd, eop)
    r_itrf = D_itrf_teme * r_teme
    return r_itrf
end

ground_facility_gaps(
    orbp,
    [(-(15 + 33 / 60) |> deg2rad, -(56 + 04 / 60) |> deg2rad, 0)];
    duration = 1 * 86400,
    f_eci_to_ecef = f_eci_to_ecef,
    minimum_elevation = 5 |> deg2rad,
    unit = :m
)
```

We can also perform analyses using multiple ground facilities. For example, let's find the
accumulated gap if we consider the INPE's stations at Cuiabá, MT, Brazil, and Alcântara, MA,
Brazil:

```@repl ground_facility_gap
ground_facility_gaps(
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

By default, the algorithm computes the accumulated gap, _i.e._, it considers the gap active
if neither station has visibility to the satellite. We can change this logic by overloading
the function in the keyword parameter `reduction`. For example, let's compute the gap when
the satellite does not have visibility to both ground stations at the same time:

```@repl ground_facility_gap
ground_facility_gaps(
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
