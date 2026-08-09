## Description #############################################################################
#
# This file contains functions for analyzing satellite decay.
#
############################################################################################

export decay_analysis

"""
    decay_analysis(orb::KeplerianElements; kwargs...) -> DataFrame

Compute the orbital decay analysis of a satellite with initial osculating mean elements
`orb` represented in the TOD reference frame, propagating the mean orbital elements with
averaged perturbations until the mean perigee altitude reaches `terminate_altitude` or the
propagation time reaches `tf`.

The model averages the following perturbations over one orbit: Earth gravity zonal
harmonics (including a closed-form J₂² correction), third-body attraction of the Sun and
the Moon, atmospheric drag using the NRLMSISE-00 model, and solar radiation pressure gated
by the Earth shadow.

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
- `gravity_model::Union{AbstractGravityModel, Nothing}`: Gravity model used to compute the
    Earth gravitational perturbation. If it is `nothing`, the system fetches and loads the
    EGM2008 model.
    (**Default**: `nothing`)
- `num_sampling_points_per_orbit::Int`: Number of sampling points used to average the
    perturbations over one orbit.
    (**Default**: 17)
- `abstol::Number`: Absolute tolerance of the numerical integration.
    (**Default**: 1e-6)
- `Ap::Union{Function, Number}`: Geomagnetic index. It can be a constant value or a
    function of time (in Julian days) that returns the geomagnetic index at that instant:
    `(jd_utc::Number) -> Number`. By default, it will use a pre-defined function that
    obtains the index from the `SpaceIndices.jl` package, which must be already initialized
    with `SpaceIndices.init()`.
    (**Default**: a function that obtains the index from **SpaceIndices.jl**)
- `C_d::Number`: Drag coefficient [-].
    (**Default**: 2.2)
- `C_r::Number`: Solar radiation pressure coefficient [-].
    (**Default**: 1.25)
- `distance_unit::Symbol`: Unit of the altitude columns in the output `DataFrame`. It can
    be `:m` for meters or `:km` for kilometers.
    (**Default**: `:km`)
- `F107::Union{Function, Number}`: 10.7 cm solar flux index [sfu]. It can be a constant
    value or a function of time (in Julian days) that returns the solar flux index at that
    instant: `(jd_utc::Number) -> Number`. By default, it will use a pre-defined function
    that obtains the index from the `SpaceIndices.jl` package, which must be already
    initialized with `SpaceIndices.init()`.
    (**Default**: a function that obtains the index from **SpaceIndices.jl**)
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

# Returns

- `DataFrame`: The mean orbital element evolution during the decay with the columns:
    - `date`: Date and time of each point [UTC] encoded using `DateTime`.
    - `time`: Elapsed time of each point since the beginning of the analysis
        [`time_unit`].
    - `f107`: 10.7 cm solar flux index used by the dynamics at each point [sfu].
    - `ap`: Geomagnetic index used by the dynamics at each point [-].
    - `mean_elements`: Mean Keplerian elements encoded using `KeplerianElements` [SI],
        where the epoch is the point date [UTC].
    - `apogee_altitude`: Mean apogee altitude [`distance_unit`].
    - `perigee_altitude`: Mean perigee altitude [`distance_unit`].
    The unit of each column is stored in the `DataFrame` using metadata. The `DataFrame`
    also stores the following table-level metadata, which is used, for example, by the
    function [`plot_decay_analysis`](@ref):
    - `Satellite Mass`: Satellite mass [kg].
    - `Satellite Mean Area`: Mean cross-sectional area [m²].
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
           F107 = 140,
           Ap = 15
       );

julia> df[end, :date]  # ..................................... Estimation of the decay epoch
2024-02-01T07:12:53.537
```
"""
function decay_analysis(::Any; kwargs...)
    return error("Load OrdinaryDiffEqAdamsBashforthMoulton.jl to use `decay_analysis`.")
end
