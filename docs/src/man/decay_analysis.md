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

It computes the orbital decay analysis of a satellite with initial mean elements `orb`
represented in the TOD reference frame, propagating the mean orbital elements with
averaged perturbations until the mean perigee altitude reaches `terminate_altitude` or the
propagation time reaches `tf`.

By default, the input elements are treated as **mean elements** with respect to the
averaged dynamics, following the same convention of semi-analytical tools such as STELA.
Osculating elements, obtained, for example, from an instantaneous state vector, can be
used by setting the keyword `input_type` to `:osculating`, in which case they are
converted to mean elements before the propagation.

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
  returns the atmospheric density [kg/m³] at a given location and time considering a set of
  space indices. It must have the signature
  `(jd_utc::Number, lat::Number, lon::Number, alt::Number, space_indices::NamedTuple) -> Number`
  where `jd_utc` is the Julian date in UTC, `lat`, `lon`, and `alt` are the
  geodetic latitude [rad], longitude [rad], and altitude [m] of the point where the density
  is evaluated, and `space_indices` is the named tuple with the space indices at that
  instant provided by the keyword `space_indices`. If it is `nothing`, the system uses an
  internal wrapper for the NRLMSISE-00 model provided by **AtmosphericModels.jl**, which
  requires the fields `f107` (daily 10.7 cm solar flux) [sfu], `f107_avg` (81-day average
  of the 10.7 cm solar flux) [sfu], and `ap` (daily geomagnetic index) [-] in the named
  tuple. The macros [`@decay_analysis__jacchia77`](@ref) and
  [`@decay_analysis__jr1971`](@ref) provide keyword sets that select the Jacchia 1977 and
  the Jacchia-Roberts 1971 models instead (see
  [Using the Jacchia Models](@ref decay_analysis_jacchia)).
  (**Default**: `nothing`)
- `atmospheric_model_name::Union{Nothing, String}`: Name of the atmospheric model recorded
  in the metadata `Atmospheric Model` of the output `DataFrame` and shown, for example, by
  [`plot_decay_analysis`](@ref). If it is `nothing`, the name is derived from the keyword
  `atmospheric_model`: `"NRLMSISE-00"` for the default model and `"Custom (<name>)"` for
  user-provided callables.
  (**Default**: `nothing`)
- `gravity_model::Union{AbstractGravityModel, Nothing}`: Gravity model used to compute the
  Earth gravitational perturbation. If it is `nothing`, the system fetches and loads the
  EGM96 model.
  (**Default**: `nothing`)
- `input_type::Symbol`: How the input elements `orb` are interpreted. If it is `:mean`,
  they are treated as mean elements with respect to the averaged dynamics. If it is
  `:osculating`, they are treated as osculating elements and converted to mean elements
  before the propagation. Any other symbol raises an `ArgumentError`.
  (**Default**: `:mean`)
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
- `space_indices::Any`: Space indices required by the atmospheric model. It can be a
  constant `NamedTuple` used for all instants or a callable object of time (in Julian days)
  that returns the named tuple with the space indices at that instant:
  `(jd_utc::Number) -> NamedTuple`. If it is `nothing`, the system provides the named tuple
  `(f107 = ..., f107_avg = ..., ap = ...)` required by the default atmospheric model using
  the data in **SpaceIndices.jl**: the observed daily F10.7 of the previous day (space
  index `F10obs`), as prescribed by the NRLMSISE-00 documentation, the observed centered
  81-day average F10.7 (space index `F10obs_avg_center81`), and the observed daily
  geomagnetic index (space index `Ap_daily`). Outside the observed timespans, the F10.7 values fall back to the predicted observed F10.7 (space index
  `F10obs_predicted`), which is a harmonic model fitted to the observed data that
  captures the mean solar cycle behavior,
  and the geomagnetic index falls back to Ap = 9, as in STELA. In this case, the required
  space index sets are initialized automatically, downloading the data files on first use.
  Notice that the F10.7 prediction is intended for long-term analyses and must not be used
  as a short-term forecast of the solar activity.
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
- `verbose::Bool`: If `true`, a progress interface is shown in `stderr` during the
  numerical integration. In interactive terminals, a live panel shows a progress bar, the
  current mean perigee and apogee altitudes, the elapsed model time, and the elapsed wall
  time; at the end, the panel is left on the screen with the final state and a summary
  line is printed below it. Otherwise, plain progress lines are printed at every 5%,
  followed by the summary line. The progress fraction is the maximum between the time
  fraction and the perigee descent fraction, so it reaches 100% at either termination
  condition. Enabling the interface does not change the analysis result.
  (**Default**: `false`)

The function returns a `DataFrame` with the mean orbital element evolution during the decay
with the columns:

- `date`: Date and time of each point [UTC] encoded using `DateTime`.
- `time`: Elapsed time of each point since the beginning of the analysis [`time_unit`].
- `space_indices`: Named tuple with the space indices used by the dynamics at each point.
  This column has no unit metadata since its fields have heterogeneous units. The default
  source provides the fields `f107` [sfu], `f107_avg` [sfu], and `ap` [-].
- `mean_elements`: Mean Keplerian elements encoded using `KeplerianElements` [SI], where
  the epoch is the point date [UTC].
- `apogee_altitude`: Mean apogee altitude [`distance_unit`].
- `perigee_altitude`: Mean perigee altitude [`distance_unit`].

The unit of each column is stored in the `DataFrame` using metadata. The `DataFrame` also
stores the table-level metadata `Satellite Mass` [kg], `Satellite Mean Area` [m²],
`Terminate Altitude` [m], `Drag Coefficient` [-], `SRP Coefficient` [-],
`Atmospheric Model`, and `Space Indices Source`, which are used, for example, by the
function [`plot_decay_analysis`](@ref).

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
indices default to the observed and predicted values provided by **SpaceIndices.jl**, and
the system fetches the EGM96 gravity model automatically:

```@repl decay_analysis
df = decay_analysis(orb; satellite_mass = 100.0, satellite_mean_area = 1.0)
```

The estimated decay epoch is the date of the last point:

```@repl decay_analysis
df[end, :date]
```

We can also provide constant space indices, which is useful, for example, to analyze
worst-case scenarios with high solar activity:

```@repl decay_analysis
df_high = decay_analysis(
    orb;
    satellite_mass = 100.0,
    satellite_mean_area = 1.0,
    space_indices = (f107 = 250.0, f107_avg = 250.0, ap = 9.0)
)

df_high[end, :date]
```

## [Using the Jacchia Models](@id decay_analysis_jacchia)

The macros [`@decay_analysis__jacchia77`](@ref) and [`@decay_analysis__jr1971`](@ref)
provide keyword sets that select the Jacchia 1977 and the Jacchia-Roberts 1971 atmospheric
models provided by **AtmosphericModels.jl** instead of the default NRLMSISE-00. Each macro
expands to the keywords `atmospheric_model`, `atmospheric_model_name`, and
`space_indices`, hence it must be used in the keyword section of the call:

```@repl decay_analysis
df_j77 = decay_analysis(
    orb;
    satellite_mass = 100.0,
    satellite_mean_area = 1.0,
    @decay_analysis__jacchia77
)

df_j77[end, :date]
```

The Jacchia models consume the space indices `f107` (daily 10.7 cm solar flux) [sfu],
`f107_avg` (81-day average of the 10.7 cm solar flux) [sfu], and `kp` (daily geomagnetic
index Kp) [-]. The Jacchia models were derived using the F10.7 flux adjusted to 1 AU,
unlike NRLMSISE-00, which uses the observed flux at the actual Earth-Sun distance. Hence,
the default space indices source selected by the macros provides the adjusted F10.7
values and the observed Kp (space indices `F10adj`, `F10adj_avg_center81`, and
`Kp_daily`), falling back to the predicted adjusted F10.7 (space index
`F10adj_predicted`) and to Kp = 7 / 3 (equivalent to Ap = 9, as in STELA) outside the
available timespans. Keywords passed **after** the macro override the ones it provides,
so we can, for example, use the Jacchia 1977 model with constant space indices:

```@repl decay_analysis
df_j77_high = decay_analysis(
    orb;
    satellite_mass = 100.0,
    satellite_mean_area = 1.0,
    @decay_analysis__jacchia77,
    space_indices = (f107 = 250.0, f107_avg = 250.0, kp = 3.0)
)

df_j77_high[end, :date]
```

!!! note

    The Jacchia 1977 model does not have a closed-form solution, so its equations are
    numerically integrated at every density evaluation, making the analysis considerably
    slower than with the default NRLMSISE-00 model. The Jacchia-Roberts 1971 model,
    selected by [`@decay_analysis__jr1971`](@ref), is a closed-form analytic fit of the
    Jacchia model family with speed comparable to NRLMSISE-00:

```@repl decay_analysis
df_jr71 = decay_analysis(
    orb;
    satellite_mass = 100.0,
    satellite_mean_area = 1.0,
    @decay_analysis__jr1971
)

df_jr71[end, :date]
```

## Plotting

If the user loads the package [Makie.jl](https://docs.makie.org/stable/), an extension is
loaded and adds the possibility to plot the decay analysis using the function
[`plot_decay_analysis`](@ref). The figure shows the evolution of the mean apogee and
perigee altitudes, a dashed line marking the terminate altitude, and an annotation marking
the reentry. The keyword `show_dates` adds the absolute dates [UTC] to the figure: the
analysis timespan in the subtitle and the estimated reentry date in the information panel
and in the reentry annotation. The information panel
shows the satellite mass, the mean area, the time to reenter, and a card with the analysis
assumptions (atmospheric model and drag and SRP coefficients), all resolved from the
`DataFrame` metadata. The keyword `show_f107`
also plots the daily and the 81-day average 10.7 cm solar flux indices used by the dynamics
using a twin y-axis. The values are extracted from the column `space_indices` using the
keywords `f107_getter` and `f107_avg_getter`, whose defaults match the named tuple provided
by the default space indices source; passing `nothing` to a getter omits the related curve:

```@example decay_analysis
using CairoMakie

CairoMakie.activate!(type = "png", px_per_unit = 2) # hide

fig, ax = plot_decay_analysis(
    df;
    mission_name = "My Mission",
    show_dates   = true,
    show_f107    = true
)

fig
```

To export the figure in high resolution for reports, use:

```julia
save("decay_analysis.png", fig; px_per_unit = 2)
```
