## Description #############################################################################
#
# This file contains functions for analyzing satellite decay.
#
############################################################################################

export decay_analysis, @decay_analysis__jacchia77, @decay_analysis__jr1971

"""
    decay_analysis(orb::KeplerianElements; kwargs...) -> DataFrame

Compute the orbital decay analysis of a satellite with initial osculating mean elements
`orb` represented in the TOD reference frame, propagating the mean orbital elements with
averaged perturbations until the mean perigee altitude reaches `terminate_altitude` or the
propagation time reaches `tf`.

The model averages the following perturbations over one orbit: Earth gravity zonal
harmonics (including a closed-form J₂² correction), third-body attraction of the Sun and
the Moon, atmospheric drag, and solar radiation pressure gated by the Earth shadow.

!!! warning

    This function **only works** after loading the package
    **OrdinaryDiffEqAdamsBashforthMoulton.jl**, which provides the default solver `VCABM`.
    Loading **OrdinaryDiffEq.jl** v6 also works because it depends on that package, but
    **OrdinaryDiffEq.jl** v7 or newer does not. In this case, the package
    **OrdinaryDiffEqAdamsBashforthMoulton.jl** must be loaded explicitly.

!!! note

    The default integrator configuration (`solver`, `reltol`, `abstol`, and
    `num_sampling_points_per_orbit`) is tuned for fast and accurate **decay lifetime**
    estimation: the lifetime and the altitude evolution change well below the atmospheric
    model uncertainty. However, the angular elements (RAAN, argument of perigee, and mean
    anomaly) accumulate a larger numerical error over long arcs. If accurate angles are
    required, use a tighter configuration, e.g. `solver = Tsit5()`, `reltol = 1e-8`,
    `abstol = 1e-8`, and `num_sampling_points_per_orbit = 33`.

# Keywords

- `satellite_mass::Number`: Satellite mass [kg]. This keyword is required.
- `satellite_mean_area::Number`: Mean cross-sectional area [m²] used for both the
    atmospheric drag and the solar radiation pressure. This keyword is required.
- `atmospheric_model::Any`: Callable object (a function or a callable structure) that
    returns the atmospheric density [kg/m³] at a given location and time considering a set
    of space indices. It must have the signature
    `(jd_utc::Number, lat::Number, lon::Number, alt::Number, space_indices::NamedTuple) -> Number`
    where `jd_utc` is the Julian date in UTC, `lat`, `lon`, and `alt` are the geodetic
    latitude [rad], longitude [rad], and altitude [m] of the point where the density is
    evaluated, and `space_indices` is the named tuple with the space indices at that
    instant provided by the keyword `space_indices`. If it is `nothing`, the system uses an
    internal wrapper for the NRLMSISE-00 model provided by **AtmosphericModels.jl**, which
    requires the fields `f107` (daily 10.7 cm solar flux) [sfu], `f107_avg` (81-day average
    of the 10.7 cm solar flux) [sfu], and `ap` (daily geomagnetic index) [-] in the named
    tuple. The macros [`@decay_analysis__jacchia77`](@ref) and
    [`@decay_analysis__jr1971`](@ref) provide keyword sets that select the Jacchia 1977
    and the Jacchia-Roberts 1971 models instead.
    (**Default**: `nothing`)
- `atmospheric_model_name::Union{Nothing, String}`: Name of the atmospheric model recorded
    in the metadata `Atmospheric Model` of the output `DataFrame` and shown, for example,
    by [`plot_decay_analysis`](@ref). If it is `nothing`, the name is derived from the
    keyword `atmospheric_model`: `"NRLMSISE-00"` for the default model and
    `"Custom (<name>)"` for user-provided callables.
    (**Default**: `nothing`)
- `gravity_model::Union{AbstractGravityModel, Nothing}`: Gravity model used to compute the
    Earth gravitational perturbation. If it is `nothing`, the system fetches and loads the
    EGM96 model.
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
- `space_indices::Any`: Space indices required by the atmospheric model. It can be a
    constant `NamedTuple` used for all instants or a callable object of time (in Julian
    days) that returns the named tuple with the space indices at that instant:
    `(jd_utc::Number) -> NamedTuple`. If it is `nothing`, the system provides the named
    tuple `(f107 = ..., f107_avg = ..., ap = ...)` required by the default atmospheric
    model using the data in **SpaceIndices.jl**: the observed daily F10.7 (space index
    `F10obs`), the observed last-81-day average F10.7 (space index `F10obs_avg_last81`),
    and the observed daily geomagnetic index (space index `Ap_daily`). Outside the observed
    timespans, the F10.7 values fall back to the predicted F10.7 (space index
    `F10predicted`), which is a harmonic model fitted to the observed data that captures
    the mean solar cycle behavior, and the geomagnetic index falls back to Ap = 9, as in
    STELA. In this case, the required space index sets are initialized automatically,
    downloading the data files on first use.
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
    numerical integration. In interactive terminals, a live panel shows a progress bar,
    the current mean perigee and apogee altitudes, the elapsed model time, and the
    elapsed wall time; at the end, the panel is left on the screen with the final state
    and a summary line is printed below it. Otherwise, plain progress lines are printed
    at every 5%, followed by the summary line. The progress fraction is the maximum
    between the time fraction and the perigee descent fraction, so it reaches 100% at
    either termination condition. Enabling the interface does not change the analysis
    result.
    (**Default**: `false`)

# Returns

- `DataFrame`: The mean orbital element evolution during the decay with the columns:
    - `date`: Date and time of each point [UTC] encoded using `DateTime`.
    - `time`: Elapsed time of each point since the beginning of the analysis
        [`time_unit`].
    - `space_indices`: Named tuple with the space indices used by the dynamics at each
        point. This column has no unit metadata since its fields have heterogeneous units.
        The default source provides the fields `f107` [sfu], `f107_avg` [sfu], and
        `ap` [-].
    - `mean_elements`: Mean Keplerian elements encoded using `KeplerianElements` [SI],
        where the epoch is the point date [UTC].
    - `apogee_altitude`: Mean apogee altitude [`distance_unit`].
    - `perigee_altitude`: Mean perigee altitude [`distance_unit`].
    The unit of each column is stored in the `DataFrame` using metadata. The `DataFrame`
    also stores the following table-level metadata, which is used, for example, by the
    function [`plot_decay_analysis`](@ref):
    - `Atmospheric Model`: Name of the atmospheric model used by the drag computation.
    - `Drag Coefficient`: Drag coefficient [-].
    - `Satellite Mass`: Satellite mass [kg].
    - `Satellite Mean Area`: Mean cross-sectional area [m²].
    - `Space Indices Source`: Description of the space indices source used by the
        dynamics.
    - `SRP Coefficient`: Solar radiation pressure coefficient [-].
    - `Terminate Altitude`: Mean perigee altitude that terminates the analysis [m].
    If the keyword `return_solution` is `true`, the function returns a tuple with the
    `DataFrame` and the raw `ODESolution`.

# Extended help

The satellite lifetime can be obtained from the last row of the returned `DataFrame`: if
the perigee altitude reached `terminate_altitude` before `tf`, the last `date` is the decay
epoch estimation.

## Examples

```julia-repl
julia> using SatelliteAnalysis, OrdinaryDiffEqAdamsBashforthMoulton

julia> jd₀ = date_to_jd(2024, 1, 1);

julia> orb = KeplerianElements(
           jd₀,
           EARTH_EQUATORIAL_RADIUS + 300e3,
           0.001,
           98.0 |> deg2rad,
           ltdn_to_raan(10.5, jd₀),
           90 |> deg2rad,
           0
       );

julia> df = decay_analysis(
           orb;
           satellite_mass = 100.0,
           satellite_mean_area = 1.0,
           space_indices = (f107 = 140.0, f107_avg = 140.0, ap = 9.0)
       );

julia> df[end, :date]  # ..................................... Estimation of the decay epoch
2024-02-02T18:56:08.476
```

If the keyword `space_indices` is omitted, the analysis uses the observed and predicted
indices provided by **SpaceIndices.jl**, requiring only the satellite properties:

```julia-repl
julia> df = decay_analysis(orb; satellite_mass = 100.0, satellite_mean_area = 1.0);
```
"""
function decay_analysis(::Any; kwargs...)
    return error("Load OrdinaryDiffEqAdamsBashforthMoulton.jl to use `decay_analysis`.")
end

"""
    @decay_analysis__jacchia77

Keyword set for [`decay_analysis`](@ref) that selects the Jacchia 1977 atmospheric model
provided by **AtmosphericModels.jl** instead of the default NRLMSISE-00. The macro expands
to the keywords `atmospheric_model`, `atmospheric_model_name`, and `space_indices`, hence
it must be used in the keyword section of the call:

```julia
decay_analysis(
    orb;
    satellite_mass = 100.0,
    satellite_mean_area = 1.0,
    @decay_analysis__jacchia77
)
```

The Jacchia 1977 model consumes the space indices `f107` (daily 10.7 cm solar flux) [sfu],
`f107_avg` (81-day average of the 10.7 cm solar flux) [sfu], and `kp` (daily geomagnetic
index Kp) [-] from the named tuple provided by the keyword `space_indices`. The macro also
selects a default space indices source with these fields: the observed daily F10.7 (space
index `F10obs`), the observed last-81-day average F10.7 (space index `F10obs_avg_last81`),
and the observed daily Kp (space index `Kp_daily`), falling back to the predicted F10.7
(space index `F10predicted`) and to Kp = 7 / 3 (equivalent to Ap = 9, as in STELA) outside
the observed timespans.

Keywords passed **after** the macro override the ones it provides. For example, the
following call uses the Jacchia 1977 model with constant space indices:

```julia
decay_analysis(
    orb;
    satellite_mass = 100.0,
    satellite_mean_area = 1.0,
    @decay_analysis__jacchia77,
    space_indices = (f107 = 140.0, f107_avg = 140.0, kp = 3.0)
)
```

!!! note

    The Jacchia 1977 model does not have a closed-form solution, so its equations are
    numerically integrated at every density evaluation, making the analysis considerably
    slower than with the default NRLMSISE-00 model. The macro
    [`@decay_analysis__jr1971`](@ref) selects the Jacchia-Roberts 1971 model, a
    closed-form analytic fit of the Jacchia model family with speed comparable to
    NRLMSISE-00.

!!! warning

    This macro **only works** after loading the package
    **OrdinaryDiffEqAdamsBashforthMoulton.jl**, as described in [`decay_analysis`](@ref).
"""
macro decay_analysis__jacchia77()
    # The expansion splats the named tuple returned by the setup function into the keyword
    # section of the call, adding the keywords `atmospheric_model`,
    # `atmospheric_model_name`, and `space_indices`.
    return Expr(:..., :(_decay_analysis__jacchia77_setup(nothing)))
end

"""
    _decay_analysis__jacchia77_setup(::Any) -> NamedTuple

Return the named tuple with the keywords of [`decay_analysis`](@ref) that select the
Jacchia 1977 atmospheric model, used by the macro [`@decay_analysis__jacchia77`](@ref).
The argument is a dummy value (`nothing`) that allows the decay analysis extension to
override this fallback with a more specific method instead of overwriting it, which is
forbidden during precompilation.
"""
function _decay_analysis__jacchia77_setup(::Any)
    return error(
        "Load OrdinaryDiffEqAdamsBashforthMoulton.jl to use `@decay_analysis__jacchia77`."
    )
end

"""
    @decay_analysis__jr1971

Keyword set for [`decay_analysis`](@ref) that selects the Jacchia-Roberts 1971 atmospheric
model provided by **AtmosphericModels.jl** instead of the default NRLMSISE-00. The macro
expands to the keywords `atmospheric_model`, `atmospheric_model_name`, and
`space_indices`, hence it must be used in the keyword section of the call:

```julia
decay_analysis(
    orb;
    satellite_mass = 100.0,
    satellite_mean_area = 1.0,
    @decay_analysis__jr1971
)
```

The Jacchia-Roberts 1971 model is a closed-form analytic fit of the Jacchia model family,
with speed comparable to NRLMSISE-00, unlike the Jacchia 1977 model selected by
[`@decay_analysis__jacchia77`](@ref), whose equations are numerically integrated at every
density evaluation.

The model consumes the space indices `f107` (daily 10.7 cm solar flux) [sfu], `f107_avg`
(81-day average of the 10.7 cm solar flux) [sfu], and `kp` (daily geomagnetic index Kp)
[-] from the named tuple provided by the keyword `space_indices`. The macro also selects a
default space indices source with these fields: the observed daily F10.7 (space index
`F10obs`), the observed last-81-day average F10.7 (space index `F10obs_avg_last81`), and
the observed daily Kp (space index `Kp_daily`), falling back to the predicted F10.7 (space
index `F10predicted`) and to Kp = 7 / 3 (equivalent to Ap = 9, as in STELA) outside the
observed timespans.

Keywords passed **after** the macro override the ones it provides. For example, the
following call uses the Jacchia-Roberts 1971 model with constant space indices:

```julia
decay_analysis(
    orb;
    satellite_mass = 100.0,
    satellite_mean_area = 1.0,
    @decay_analysis__jr1971,
    space_indices = (f107 = 140.0, f107_avg = 140.0, kp = 3.0)
)
```

!!! warning

    This macro **only works** after loading the package
    **OrdinaryDiffEqAdamsBashforthMoulton.jl**, as described in [`decay_analysis`](@ref).
"""
macro decay_analysis__jr1971()
    # The expansion splats the named tuple returned by the setup function into the keyword
    # section of the call, adding the keywords `atmospheric_model`,
    # `atmospheric_model_name`, and `space_indices`.
    return Expr(:..., :(_decay_analysis__jr1971_setup(nothing)))
end

"""
    _decay_analysis__jr1971_setup(::Any) -> NamedTuple

Return the named tuple with the keywords of [`decay_analysis`](@ref) that select the
Jacchia-Roberts 1971 atmospheric model, used by the macro
[`@decay_analysis__jr1971`](@ref). The argument is a dummy value (`nothing`) that allows
the decay analysis extension to override this fallback with a more specific method instead
of overwriting it, which is forbidden during precompilation.
"""
function _decay_analysis__jr1971_setup(::Any)
    return error(
        "Load OrdinaryDiffEqAdamsBashforthMoulton.jl to use `@decay_analysis__jr1971`."
    )
end
