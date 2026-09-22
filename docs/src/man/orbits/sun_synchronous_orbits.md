# Sun-Synchronous Orbits

```@meta
CurrentModule = SatelliteAnalysis
```

```@setup sun_sync
using SatelliteAnalysis
```

Given the Earth's gravitational potential, all the low Earth orbits (LEO) suffer from
perturbations in their elements. The right ascension of the ascending node (RAAN) is one of
the elements that see a secular perturbation. Hence, we can use this effect to design
orbits that the RAAN time derivative matches that of the Earth's orbit around the Sun. In
this case, the orbit plane keeps its geometry almost constant regarding the Sun vector.
Hence, the apparent local time of the ascending or descending node (LTAN and LTDN) will also
be almost invariable, as shown in the following figure:

```@raw html
<div align="center">
  <img src="../../../assets/sun_sync_orbit_plane.png" alt="Orbit Plane in Sun-Synchronous Orbits" width="100%"/>
</div>
```

In fact, if the Earth was a perfect sphere, no secular perturbation would be seen in RAAN.
Hence, the LTAN would have a variation of 24 hours in a year.

Considering only terms up to ``J_2``, which is the dominant effect, the RAAN time derivative
is:

```math
\frac{d\Omega}{dt} = -\frac{3}{2} J_2 \left(\frac{R_0}{p_0}\right)^2 \bar{n} \cos i_0\ ,
```

where ``R_0`` is the Earth's Equatorial radius, ``i_0`` is the orbit inclination, ``p_0 =
a_0 (1 - e_0^2)``, ``a_0`` is the orbit semi-major axis, and ``e_0`` is the orbit
eccentricity, and ``\bar{n}`` is the perturbed mean-motion, given by:

```math
\bar{n} = n_0 \left[1 + \frac{3}{4} J_2 \left(\frac{R_0}{p_0}\right)^2 \sqrt{1 - e_0^2} \left(2 - 3\sin^2 i_0\right)\right]\ ,
```

where ``n_0 = \sqrt{\mu / a^3}``, and ``\mu`` is the Earth's standard gravitational
parameter.

!!! note

    Formally, all quantities on the right-hand side of those equations must be the mean
    elements instead of the initial ones. However, we are considering only the secular
    perturbations caused by the ``J_2`` term. Thus, the semi-major axis, eccentricity, and
    inclination do not suffer from secular effects.

Finally, we can design a Sun-synchronous orbit by selecting the semi-major axis and
inclination that leads to:

```math
\frac{d\Omega}{dt} = 0.9856473598947981^\circ / \text{day}\ .
```

## Designing Sun-Synchronous Orbits from Angular Velocity

The satellite orbital angular velocity ``\omega_o`` is:

```math
\omega_{o}(a_0, i_0) = \frac{dM}{dt} + \frac{d\omega}{dt}\ ,
```

where ``M`` is the mean anomaly, and ``\omega`` is the argument of perigee. Considering only
the secular effects caused by the ``J_2`` term, one gets:

```math
\begin{aligned}
  \frac{dM}{dt} &= \bar{n}\ , \\
  \frac{d\omega}{dt} &= \frac{3}{4} J_2 \left(\frac{R_0}{p_0}\right)^2 \bar{n} \left(4 - 5\sin^2 i_0\right)\ .
\end{aligned}
```

Thus, given a desired angular velocity ``n_d``, we can find the semi-major axis and
inclination that leads to a Sun-synchronous orbit by numerically solving the system:

```math
\begin{aligned}
  \frac{d\Omega}{dt}(a_0, i_0) &= 0.9856473598947981^\circ / \text{day}\ , \\
  \omega_o(a_0, i_0) &= n_d\ .
\end{aligned}
```

We can use the function `sun_sync_orbit_from_angular_velocity`:

```@docs; canonical = false
sun_sync_orbit_from_angular_velocity
```

### Example

Let's say we want to compute the Sun-synchronous orbit for a mission that must perform
exactly 14 orbits per day:

```@repl sun_sync
n_d = 14 * (2π / 86400)
a, i, converged = sun_sync_orbit_from_angular_velocity(n_d)
a / 1000
rad2deg(i)
```

## Designing Sun-Synchronous Orbits from Semi-Major Axis

Given a desired semi-major axis ``a_d``, we can compute the inclination that turns the orbit
into a Sun-synchronous one by solving numerically:

```math
\frac{d\Omega}{dt}(a_d, i_0) = 0.9856473598947981^\circ / \text{day}\ .
```

We can use the function `sun_sync_orbit_inclination`:

```@docs; canonical = false
sun_sync_orbit_inclination
```

### Example

Let's find the inclination that turns an orbit with semi-major axis 6819 km and eccentricity
0.0015 into a Sun-synchronous one:

```@repl sun_sync
i, converged = sun_sync_orbit_inclination(6819e3, 0.0015)
rad2deg(i)
```

## Designing Sun-Synchronous Orbits from Inclination


Given a desired inclination ``i_d``, we can compute the semi-major axis that turns the orbit
into a Sun-synchronous one by solving numerically:

```math
\frac{d\Omega}{dt}(a_0, i_d) = 0.9856473598947981^\circ / \text{day}\ .
```

We can use the function `sun_sync_orbit_semi_major_axis`:

```@docs; canonical = false
sun_sync_orbit_semi_major_axis
```

### Example

Let's find the semi-major axis that turns an orbit with inclination axis 98.190° and
eccentricity 0.001987 into a Sun-synchronous one:

```@repl sun_sync
a, converged = sun_sync_orbit_semi_major_axis(98.190 |> deg2rad, 0.001987)
a / 1000
```

## Designing Sun-Synchronous, Ground-Repeating Orbits

The most common type of a Sun-synchronous orbit is a ground-repeating one. In this case, we
have a Sun-synchronous orbit, as mentioned before, whose ground track repeats after a finite
number of (solar) days. The ground-repeating condition happens if the satellite performs a
rational number of orbits per solar day:

```math
R_d = I + \frac{N}{D}\ ,
```

where ``R_d`` is the number of revolutions per day, and ``I, N, D \in \mathbb{N}``. If the
greatest common divisor of ``N`` and ``D`` is one, the ground track of such an orbit repeats
after ``ID + N`` revolutions, or ``D`` solar days.

If an orbit has ``R_d`` revolutions per solar day, we can compute its angular velocity as
follows:

```math
\omega_o = R_d \frac{2\pi}{86400}\ .
```

Hence, we can design a Sun-synchronous, ground-repeating orbit by numerically solving the
system:

```math
\begin{aligned}
  \frac{d\Omega}{dt}(a_0, i_0) &= 0.9856473598947981^\circ / \text{day}\ , \\
  \omega_o(a_0, i_0) &= R_d \frac{2\pi}{86400}\ ,
\end{aligned}
```

using the same method we described in the Section [Designing Sun-Synchronous Orbits from
Angular Velocity](@ref).

We can use the function `design_sun_sync_ground_repeating_orbit`:

```@docs; canonical = false
design_sun_sync_ground_repeating_orbit
```

### Example

Let's find all the possible orbits between 650km and 800km that repeat the ground track in,
at most, 5 days:

```@repl sun_sync
df = design_sun_sync_ground_repeating_orbit(
    1,
    5;
    minimum_altitude = 650e3,
    maximum_altitude = 800e3
)
show(df; allcols = true, allrows = true)
```

!!! info

    The designer can use the fields `adjacent_gt_distance` and `adjacent_gt_angle` to check
    whether the mission payload can operate correctly in the orbit. For example, in orbit
    \#4, the payload swath of a remote sensing satellite with a camera must be higher than
    543.811 km. Otherwise, there will be gaps in the images.
