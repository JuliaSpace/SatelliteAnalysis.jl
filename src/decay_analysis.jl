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

    This function **only works** after loading the package **OrdinaryDiffEq.jl**.

!!! note

    If `F107` or `Ap` are `nothing`, the space indices are obtained using the function
    `space_index`, which requires calling `SpaceIndices.init()` beforehand.

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
    (**Default**: 33)
- `Ap::Union{Nothing, Number}`: Geomagnetic index. If it is `nothing`, the value is
    obtained from the initialized space indices at each instant.
    (**Default**: `nothing`)
- `C_d::Number`: Drag coefficient [-].
    (**Default**: 2.2)
- `C_r::Number`: Solar radiation pressure coefficient [-].
    (**Default**: 1.25)
- `F107::Union{Nothing, Number}`: 10.7 cm solar flux index. If it is `nothing`, the value
    is obtained from the initialized space indices at each instant.
    (**Default**: `nothing`)
- `return_solution::Bool`: If `true`, the function also returns the raw solution of the
    numerical integration (see `OrdinaryDiffEq.ODESolution`), whose state vector is the
    equinoctial orbital elements `[a, ψ, e_x, e_y, i_x, i_y]`.
    (**Default**: `false`)
- `terminate_altitude::Number`: Mean perigee altitude [m] that terminates the analysis.
    (**Default**: 120e3)
- `tf::Number`: Maximum propagation time [s] after the orbit epoch.
    (**Default**: `30 * 365.25 * 86400`, or 30 years)

# Returns

- `DataFrame`: The mean orbital element evolution during the decay with the columns:
    - `date`: Date and time of each point [UTC] encoded using `DateTime`.
    - `time`: Elapsed time of each point since the beginning of the analysis [s].
    - `f107`: 10.7 cm solar flux index used by the dynamics at each point [sfu].
    - `ap`: Geomagnetic index used by the dynamics at each point [-].
    - `mean_elements`: Mean Keplerian elements encoded using `KeplerianElements` [SI],
        where the epoch is the point date [UTC].
    - `apogee_altitude`: Mean apogee altitude [m].
    - `perigee_altitude`: Mean perigee altitude [m].
    The unit of each column is stored in the `DataFrame` using metadata. If the keyword
    `return_solution` is `true`, the function returns a tuple with the `DataFrame` and the
    raw `ODESolution`.

# Extended help

The satellite lifetime can be obtained from the last row of the returned `DataFrame`: if
the perigee altitude reached `terminate_altitude` before `tf`, the last `date` is the decay
epoch estimation.

## Examples

```julia-repl
julia> using SatelliteAnalysis, OrdinaryDiffEq

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

julia> df[end, :date]  # ................................... Estimation of the decay epoch
2024-02-01T07:53:26.628
```
"""
function decay_analysis(::Any)
    error("Load OrdinaryDiffEq.jl to use `decay_analysis`.")
end
