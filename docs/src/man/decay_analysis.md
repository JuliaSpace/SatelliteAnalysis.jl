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

We can perform the decay analysis of a satellite using the function `decay_analysis`:

```@docs; canonical = false
decay_analysis
```

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

The macros [`@decay_analysis__jacchia77`](@ref),
[`@decay_analysis__jacchia77_stela`](@ref), and [`@decay_analysis__jr1971`](@ref) provide
keyword sets that select the Jacchia 1977 (report and STELA variants) and the
Jacchia-Roberts 1971 atmospheric models provided by **AtmosphericModels.jl** instead of
the default NRLMSISE-00. Each macro expands to the keywords `atmospheric_model`,
`atmospheric_model_name`, and `space_indices`, hence it must be used in the keyword
section of the call:

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

The macro [`@decay_analysis__jacchia77_stela`](@ref) selects the STELA variant of the
Jacchia 1977 model, which replicates the simplified assembly used by the CNES tools STELA
and PATRIUS. This variant produces total densities a few percent higher on average than
the report formulation, allowing the reproduction of decay analyses performed with those
tools: the decay time of a 500 km sun-synchronous satellite computed by STELA is
reproduced within about 1 %, whereas the report formulation yields a decay time about 7 %
longer:

```@repl decay_analysis
df_j77_stela = decay_analysis(
    orb;
    satellite_mass = 100.0,
    satellite_mean_area = 1.0,
    @decay_analysis__jacchia77_stela
)

df_j77_stela[end, :date]
```

Each macro has a function version that returns the same keywords as a named tuple:
[`decay_analysis__jacchia77_kwargs`](@ref),
[`decay_analysis__jacchia77_stela_kwargs`](@ref), and
[`decay_analysis__jr1971_kwargs`](@ref). They allow selecting the atmospheric model
programmatically, for example, when comparing the models in a loop:

```julia
for model_kwargs in (decay_analysis__jacchia77_kwargs, decay_analysis__jr1971_kwargs)
    df = decay_analysis(
        orb;
        satellite_mass      = 100.0,
        satellite_mean_area = 1.0,
        model_kwargs()...
    )

    println(df[end, :date])
end
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
