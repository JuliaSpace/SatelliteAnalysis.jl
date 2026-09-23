# Frozen Orbits

```@meta
CurrentModule = SatelliteAnalysis
```

```@setup frozen_orbit
using SatelliteAnalysis
```

Due to the Earth's gravitational perturbation, the orbit of a satellite will experience
secular changes in the argument of perigee. Hence, the satellite mean altitude per latitude
will differ during the mission. This effect can be problematic, especially if we must
compare images by a camera onboard the satellite in different periods. The altitude
variation will change the resolution, leading to some problems when comparing the data.

We can avoid this problem if we compute an eccentricity ``e`` and the argument of perigee
``\omega`` that yields theoretically:

```math
\begin{equation*}
  \frac{de}{dt} = 0,\ \frac{d\omega}{dt} = 0\ .
\end{equation*}
```

This orbit is called **frozen**. Refer to **[1]** for more information.

We can compute the eccentricity and the argument of perigee of a frozen orbit using the
function `frozen_orbit`:

```@docs; canonical = false
frozen_orbit
```

## Examples

We will compute the eccentricity and argument of perigee that yields a frozen orbit using the
data from Amazonia-1 mission. First, we will use only up to degree 5, and the default gravity
model (EGM96):

```@repl frozen_orbit
frozen_orbit(7130.982e3, 98.410 |> deg2rad; max_degree = 5)

e, ω = frozen_orbit(7130.982e3, 98.410 |> deg2rad; max_degree = 5)

e

ω |> rad2deg
```

If we want to use all the 360 degrees in EGM96, which is selected by `max_degree = 0`, we
need to increase the precision of `BigFloat` to keep the accuracy:

```@repl frozen_orbit
setprecision(1024)

e, ω = frozen_orbit(7130.982e3, 98.410 |> deg2rad; max_degree = 0)

e

ω |> rad2deg
```

We can use a different gravity model as follows:

```@repl frozen_orbit
jgm3 = GravityModels.load(IcgemFile, fetch_icgem_file(:JGM3))

e, ω = frozen_orbit(7130.982e3, 98.410 |> deg2rad; max_degree = 70, gravity_model = jgm3)

e

ω |> rad2deg
```

## References

- **[1]** **Rosborough, G. W.; Ocampo, C. A (1991)**. _Influence of higher degree zonals on
  the frozen orbit geometry_. Proceedings of the **AAS/AIAA Astrodynamics Conference**,
  Durango, CO.
