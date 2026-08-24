# Decay Analysis

```@meta
CurrentModule = SatelliteAnalysis
```

```@setup decay_analysis
using SatelliteAnalysis
using OrdinaryDiffEqAdamsBashforthMoulton
```

The decay analysis estimates how the orbit of a satellite evolves under perturbations until
it reenters the atmosphere. This information is paramount for mission design since it
provides the expected orbital lifetime, which is required, for example, by space debris
mitigation regulations.

We can perform the decay analysis of a satellite using the function:

```julia
decay_analysis(orb::KeplerianElements; kwargs...) -> DataFrame
```

It computes the orbital decay analysis of a satellite with initial osculating mean elements
`orb` represented in the TOD reference frame, propagating the mean orbital elements with
averaged perturbations until the mean perigee altitude reaches `terminate_altitude` or the
propagation time reaches `tf`.

The model averages the following perturbations over one orbit: Earth gravity zonal harmonics
(including a closed-form J₂² correction), third-body attraction of the Sun and the Moon,
atmospheric drag using a configurable atmospheric model (NRLMSISE-00 by default), and solar
radiation pressure gated by the Earth shadow.

!!! warning

    This function **only works** after loading the package
    **OrdinaryDiffEqAdamsBashforthMoulton.jl**, which provides the default solver `VCABM`.
    Loading **OrdinaryDiffEq.jl** v6 also works because it depends on that package, but
    **OrdinaryDiffEq.jl** v7 or newer does not. In this case, the package
    **OrdinaryDiffEqAdamsBashforthMoulton.jl** must be loaded explicitly.

The following keywords are available:

- `satellite_mass::Number`: Satellite mass [kg]. This keyword is required.
- `satellite_mean_area::Number`: Mean cross-sectional area [m²] used for both the
  atmospheric drag and the solar radiation pressure. This keyword is required.
- `atmospheric_model::Any`: Callable object (a function or a callable structure) that
  returns the atmospheric density [kg/m³] at a given location and time considering a
  specific F10.7 index. It must have the signature
  `(jd_utc::Number, lat::Number, lon::Number, alt::Number, F107::Number) -> Number`
  where `jd_utc` is the Julian date in UTC and `lat`, `lon`, and `alt` are the
  geodetic latitude [rad], longitude [rad], and altitude [m] of the point where the density
  is evaluated, and `F107` is the 10.7 cm solar flux index [sfu] at that instant. The latter
  must be considered as the daily value and also the 81-day centered average. By default,
  the system uses an internal wrapper for the NRLMSISE-00 model provided by
  **AtmosphericModels.jl** with a constant geomagnetic index Ap = 9, as in STELA.
  (**Default**: `_decay_analysis__nrlmsise00`)
- `gravity_model::Union{AbstractGravityModel, Nothing}`: Gravity model used to compute the
  Earth gravitational perturbation. If it is `nothing`, the system fetches and loads the
  EGM2008 model.
  (**Default**: `nothing`)
- `num_sampling_points_per_orbit::Int`: Number of sampling points used to average the
  perturbations over one orbit.
  (**Default**: 17)
- `abstol::Number`: Absolute tolerance of the numerical integration.
  (**Default**: 1e-6)
- `C_d::Number`: Drag coefficient [-].
  (**Default**: 2.2)
- `C_r::Number`: Solar radiation pressure coefficient [-].
  (**Default**: 1.25)
- `distance_unit::Symbol`: Unit of the altitude columns in the output `DataFrame`. It can
  be `:m` for meters or `:km` for kilometers.
  (**Default**: `:km`)
- `F107::Any`: 10.7 cm solar flux index [sfu]. It can be a constant value or a callable
  object of time (in Julian days) that returns the solar flux index at that instant:
  `(jd_utc::Number) -> Number`. If it is `nothing`, the system uses the
  predicted F10.7 provided by **SpaceIndices.jl** (space index `F10predicted`), a harmonic
  model fitted to the observed data that captures the mean solar cycle behavior. In this
  case, the required space index set is initialized automatically, downloading the
  coefficient file on first use. Notice that this prediction is intended for long-term
  analyses and must not be used as a short-term forecast of the solar activity.
  (**Default**: `nothing`)
- `reltol::Number`: Relative tolerance of the numerical integration.
  (**Default**: 1e-6)
- `return_solution::Bool`: If `true`, the function also returns the raw solution of the
  numerical integration (see `SciMLBase.ODESolution`), whose state vector is the
  equinoctial orbital elements `[a, ψ, e_x, e_y, i_x, i_y]`.
  (**Default**: `false`)
- `solver`: Solver from the **OrdinaryDiffEq.jl** ecosystem used for the numerical
  integration. Notice that the user must load the package that provides the selected
  solver (for example, `Tsit5` requires **OrdinaryDiffEqTsit5.jl** or
  **OrdinaryDiffEq.jl**).
  (**Default**: `VCABM()`)
- `terminate_altitude::Number`: Mean perigee altitude [m] that terminates the analysis.
  (**Default**: 120e3)
- `tf::Number`: Maximum propagation time [s] after the orbit epoch.
  (**Default**: `30 * 365.25 * 86400`, or 30 years)
- `time_unit::Symbol`: Unit of the column `time` in the output `DataFrame`. It can be `:s`
  for seconds, `:m` for minutes, `:h` for hours, `:d` for days, or `:y` for Julian years
  (365.25 days).
  (**Default**: `:y`)

The function returns a `DataFrame` with the mean orbital element evolution during the decay
with the columns:

- `date`: Date and time of each point [UTC] encoded using `DateTime`.
- `time`: Elapsed time of each point since the beginning of the analysis [`time_unit`].
- `f107`: 10.7 cm solar flux index used by the dynamics at each point [sfu].
- `mean_elements`: Mean Keplerian elements encoded using `KeplerianElements` [SI], where
  the epoch is the point date [UTC].
- `apogee_altitude`: Mean apogee altitude [`distance_unit`].
- `perigee_altitude`: Mean perigee altitude [`distance_unit`].

The unit of each column is stored in the `DataFrame` using metadata. The `DataFrame` also
stores the table-level metadata `Satellite Mass` [kg], `Satellite Mean Area` [m²], and
`Terminate Altitude` [m], which is used, for example, by the function
[`plot_decay_analysis`](@ref).

The satellite lifetime can be obtained from the last row of the returned `DataFrame`: if
the perigee altitude reached `terminate_altitude` before `tf`, the last `date` is the decay
epoch estimation.

## Examples

We will estimate the orbital lifetime of a 100 kg satellite with a mean cross-sectional
area of 1 m² in a Sun-synchronous orbit with an altitude of 300 km. The first thing we need
to do is define the orbit:

```@repl decay_analysis
jd₀ = date_to_jd(2024, 1, 1)

orb = KeplerianElements(
    jd₀,
    EARTH_EQUATORIAL_RADIUS + 300e3,
    0.001,
    98.0 |> deg2rad,
    ltdn_to_raan(10.5, jd₀),
    90 |> deg2rad,
    0
)
```

Now, we can use the function `decay_analysis` to obtain the orbit evolution until the
reentry. Notice that we only need to provide the satellite mass and mean area: the space
index defaults to the predicted F10.7, and the system fetches the EGM2008 gravity model
automatically:

```@repl decay_analysis
df = decay_analysis(orb; satellite_mass = 100.0, satellite_mean_area = 1.0)
```

The estimated decay epoch is the date of the last point:

```@repl decay_analysis
df[end, :date]
```

We can also provide a constant F10.7, which is useful, for example, to analyze worst-case
scenarios with high solar activity:

```@repl decay_analysis
df_high = decay_analysis(
    orb;
    satellite_mass = 100.0,
    satellite_mean_area = 1.0,
    F107 = 250
)

df_high[end, :date]
```

## Plotting

If the user loads the package [Makie.jl](https://docs.makie.org/stable/), an extension is
loaded and adds the possibility to plot the decay analysis using the function
[`plot_decay_analysis`](@ref). The figure shows the evolution of the mean apogee and
perigee altitudes together with an information panel. The keyword `show_f107` also plots
the 10.7 cm solar flux index used by the dynamics using a twin y-axis:

```@example decay_analysis
using CairoMakie

CairoMakie.activate!(type = "png") # hide

fig, ax = plot_decay_analysis(df; mission_name = "My Mission", show_f107 = true)

fig
```
