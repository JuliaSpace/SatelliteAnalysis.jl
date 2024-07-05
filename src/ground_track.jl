## Description #############################################################################
#
# Functions to compute the satellite ground track.
#
############################################################################################

export ground_track, ground_track_inclination

"""
    ground_track(orbp::OrbitPropagator; kwargs...) -> Vector{NTuple{2, Float64}}

Compute the satellite ground track using the orbit propagator `orbp`. It returns a vector of
`NTuple{2, Float64}` where the first element is the latitude [rad] and the second is the
longitude [rad] of each point in the ground track.

# Keywords

- `add_nans::Bool`: If `true`, we add `NaN` if there is a discontinuity in the ground track
    to improve plotting.
    (**Default**: true)
- `duration::Number`: Duration of the analysis.
    (**Default**: 86400)
- `initial_time::Number`: Initial time regarding the orbit propagator `orbp` epoch [s].
    (**Default**: 0)
- `f_eci_to_ecef::Function`: Function to convert the orbit propagator position represented
    in the Earth-centered inertial (ECI) reference frame to the Earth-centered, Earth-fixed
    (ECEF) reference frame. The signature must be

    `f_eci_to_ecef(r_i::AbstractVector, jd::Number) -> AbstractVector`

    and it must return the position vector `r_i` represented in the ECEF at the instant `jd`
    [Julian Day]. By default, we use TEME as the ECI and PEF as the ECEF.
    (**Default**: `_ground_track_default_eci_to_ecef`)
- `step::Union{Nothing, Number}`: Step for the computation. If `nothing`, we will roughly
    compute the step to approximate 1° in the mean anomaly.
    (**Default**: `nothing`)
- `track_types::Symbol`: A symbol describing what kind of track types we must add to the
    output vector. It can be `:ascending` for only ascending passages, `:descending` for
    only descending passages, or `:all` for both.
    (**Default**: `:all`)

# Extended Help

## Examples

```julia-repl
julia> using SatelliteAnalysis, UnicodePlots

julia> jd₀ = date_to_jd(2024, 1, 1);

julia> orb = KeplerianElements(
           jd₀,
           7130.982e3,
           0.001111,
           98.405 |> deg2rad,
           ltdn_to_raan(10.5, jd₀),
           π / 2,
           0
       );

julia> orbp = Propagators.init(Val(:J2), orb);

julia> gt = ground_track(orbp);

julia> lineplot(last.(gt), first.(gt))
      ┌────────────────────────────────────────┐
    2 │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
      │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
      │⠀⠀⠀⠀⣸⡿⠿⣻⠿⣿⠿⣟⡯⢟⡿⢿⣿⡿⣿⡟⣿⡿⢿⣿⢿⣿⢻⣿⠿⣿⡟⣿⡿⣿⡿⢿⠀⠀⠀⠀│
      │⠀⠀⠀⠀⠞⠹⡼⠙⣶⠋⢶⠋⢣⠏⢳⡞⠀⣿⡁⡽⡇⣹⡇⢨⢏⢈⢯⠀⣿⡀⡽⡁⣹⡇⢸⣯⠀⠀⠀⠀│
      │⠀⠀⠀⠀⢦⢰⢳⢠⢿⡀⡟⡄⡞⡇⡸⣇⢸⠁⢷⠃⣷⡇⢸⡏⠸⡼⠈⣾⠀⣷⠃⢣⠇⢹⣏⠏⠀⠀⠀⠀│
      │⠀⠀⠀⠀⢸⡞⠘⣾⠀⣿⠁⢧⠃⢹⡇⢸⡏⠀⣼⠀⣿⡆⢰⡇⢀⣇⠀⣷⠀⣾⠀⣸⡀⢸⣿⠀⠀⠀⠀⠀│
      │⠀⠀⠀⠀⢀⡇⠀⣇⠀⣿⠀⢸⠀⢸⡁⢨⡇⠀⡇⡇⡏⡇⣸⢳⢸⢸⢰⢻⣀⡏⡆⡇⡇⡜⣯⠀⠀⠀⠀⠀│
      │⠤⠤⠤⠤⢼⢧⢤⢿⠤⣿⠤⡿⡦⡼⡧⢼⣧⢴⠥⢧⡧⢼⡧⢼⡾⠼⣼⠤⣿⠤⣧⠧⢷⣧⢿⠤⠤⠤⠤⠤│
      │⠀⠀⠀⠀⡼⢸⢸⠸⣼⠁⣇⠇⡇⡇⢳⡞⢸⢸⠀⢸⡇⢸⡇⠈⡇⠀⡏⠀⣿⠀⢸⠀⢸⢹⠘⡆⠀⠀⠀⠀│
      │⠀⠀⠀⠀⠁⠈⡏⠀⡿⠀⢿⠀⢹⠃⢸⡇⠘⡇⠀⡜⡇⣸⡇⢸⢧⢠⢿⠀⣿⠀⡾⡄⡸⡟⠀⣷⠀⠀⠀⠀│
      │⠀⠀⠀⠀⢆⢰⢳⢀⢿⠀⡿⡄⡼⡆⣸⣇⢰⢧⢠⠇⣧⡇⢹⡞⠸⡼⠘⣾⠁⣷⠃⢧⢧⢷⢀⡟⠀⠀⠀⠀│
      │⠀⠀⠀⠀⢘⣎⢈⣞⠀⣿⡁⣳⡃⣹⡇⢸⣏⢘⣞⢀⡿⢆⡼⢧⣰⢳⣠⠿⣀⠾⣄⡞⣆⢈⣿⠀⠀⠀⠀⠀│
      │⠀⠀⠀⠀⢹⣾⣿⣾⣿⣶⣿⣧⣿⣷⣿⣷⣾⣿⣼⣿⣷⣾⣷⣾⣵⣲⣽⣶⣿⣶⣯⣶⣼⣿⣾⣿⠀⠀⠀⠀│
      │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
   -2 │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
      └────────────────────────────────────────┘
      ⠀-4⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀4⠀

julia> gt = ground_track(orbp, track_types = :ascending);

julia> lineplot(last.(gt), first.(gt))
      ┌────────────────────────────────────────┐
    2 │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
      │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
      │⠀⠀⠀⠀⢸⡛⠽⡛⠿⣟⠫⣟⠫⢍⠛⠻⢿⡻⢽⡛⡿⡛⠿⣟⠯⣟⠫⢟⠻⢿⡛⢽⡛⠽⡛⠯⠀⠀⠀⠀│
      │⠀⠀⠀⠀⠀⠹⡄⠘⡆⠈⢆⠈⢣⠀⢳⡀⠀⢳⡀⠹⡇⠙⡆⠈⢆⠈⢧⠀⢳⡀⠱⡀⠙⡄⠘⣆⠀⠀⠀⠀│
      │⠀⠀⠀⠀⢦⠀⢱⠀⢸⡀⠘⡄⠈⡇⠀⣇⠀⠀⢧⠀⣷⠀⢸⠀⠸⡄⠈⡆⠀⣇⠀⢣⠀⢹⠀⠈⠀⠀⠀⠀│
      │⠀⠀⠀⠀⢸⡀⠘⡆⠀⡇⠀⢇⠀⢹⠀⢸⠀⠀⢸⠀⡟⡆⠀⡇⠀⣇⠀⢳⠀⢸⠀⠸⡀⠈⡇⠀⠀⠀⠀⠀│
      │⠀⠀⠀⠀⠀⡇⠀⣇⠀⢳⠀⢸⠀⠸⡀⠈⡇⠀⠀⡇⡇⡇⠀⢳⠀⢸⠀⢸⡀⠈⡆⠀⡇⠀⢧⠀⠀⠀⠀⠀│
      │⠤⠤⠤⠤⠤⢧⠤⢼⠤⢼⠤⠼⡦⠤⡧⠤⣧⠤⠤⢧⡧⢼⠤⢼⠤⠼⡤⠤⡧⠤⡧⠤⢷⠤⢼⠤⠤⠤⠤⠤│
      │⠀⠀⠀⠀⠀⢸⠀⠸⡄⠀⡇⠀⡇⠀⢳⠀⢸⠀⠀⢸⡇⠸⡄⠈⡇⠀⡇⠀⢧⠀⢸⠀⢸⠀⠘⡆⠀⠀⠀⠀│
      │⠀⠀⠀⠀⠀⠈⡇⠀⡇⠀⢧⠀⢹⠀⢸⡀⠘⡆⠀⠘⡇⠀⡇⠀⢧⠀⢹⠀⢸⠀⠸⡄⠈⡇⠀⣇⠀⠀⠀⠀│
      │⠀⠀⠀⠀⢆⠀⢳⠀⢸⠀⠸⡄⠘⡆⠀⣇⠀⢧⠀⠀⣧⠀⢹⠀⠸⡀⠘⡆⠀⣇⠀⢧⠀⢳⠀⠘⠀⠀⠀⠀│
      │⠀⠀⠀⠀⠘⣆⠈⢆⠀⢧⡀⢳⡀⠹⡄⠘⣄⠘⣆⠀⡏⢆⠀⢧⠀⢳⡀⠹⡀⠸⣄⠘⣆⠈⢧⠀⠀⠀⠀⠀│
      │⠀⠀⠀⠀⠐⣬⣗⣮⣷⣦⣵⣢⣽⣲⣽⣶⣬⣗⣬⣗⣧⣬⣓⣦⣵⣢⣽⣲⣽⣶⣬⣖⣬⣗⣮⡵⠀⠀⠀⠀│
      │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
   -2 │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
      └────────────────────────────────────────┘
      ⠀-4⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀4⠀

julia> gt = ground_track(orbp, track_types = :descending);

julia> lineplot(last.(gt), first.(gt))
      ┌────────────────────────────────────────┐
    2 │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
      │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
      │⠀⠀⠀⠀⣸⡽⠟⣻⠽⣻⠿⢛⡯⢛⡯⢟⡻⠝⣻⠝⣿⠽⢛⠯⢛⡯⢛⡿⠟⡻⠝⣻⠽⣻⠿⢓⠀⠀⠀⠀│
      │⠀⠀⠀⠀⠞⠀⡼⠁⣰⠃⢰⠋⢠⠏⢀⡞⠀⡜⠁⡼⡇⣰⠃⢠⠏⢀⠎⠀⡞⠀⡼⠁⣰⠃⢰⣫⠀⠀⠀⠀│
      │⠀⠀⠀⠀⠀⢰⠃⢠⠇⠀⡏⠀⡞⠀⡸⠀⢸⠁⢰⠃⣇⡇⠀⡏⠀⡼⠀⣸⠀⢰⠃⢠⠇⠀⣏⠇⠀⠀⠀⠀│
      │⠀⠀⠀⠀⠀⡞⠀⣸⠀⢸⠁⢠⠃⠀⡇⠀⡏⠀⡼⠀⣿⠀⢰⠁⢀⡇⠀⡇⠀⡞⠀⣸⠀⢸⢸⠀⠀⠀⠀⠀│
      │⠀⠀⠀⠀⢀⡇⠀⡇⠀⡜⠀⢸⠀⢸⠁⢠⠇⠀⡇⠀⡏⠀⣸⠀⢸⠀⢰⠃⢀⡇⠀⡇⠀⡜⡏⠀⠀⠀⠀⠀│
      │⠤⠤⠤⠤⢼⠤⢤⠧⠤⡧⠤⡯⠤⡼⠤⢼⠤⢴⠥⢤⡧⠤⡧⠤⡾⠤⢼⠤⢼⠤⢤⠧⠤⣧⠧⠤⠤⠤⠤⠤│
      │⠀⠀⠀⠀⡼⠀⢸⠀⢸⠁⢀⠇⠀⡇⠀⡞⠀⢸⠀⢸⡇⢠⠃⠀⡇⠀⡏⠀⡼⠀⢸⠀⢸⢹⠀⠀⠀⠀⠀⠀│
      │⠀⠀⠀⠀⠁⠀⡎⠀⡸⠀⢸⠀⢰⠃⢀⡇⠀⡇⠀⡜⡇⣸⠀⢸⠁⢠⠇⠀⡇⠀⡎⠀⡸⡞⠀⣰⠀⠀⠀⠀│
      │⠀⠀⠀⠀⠀⢰⠃⢀⠇⠀⡏⠀⡼⠀⣸⠀⢰⠁⢠⠇⣇⡇⠀⡞⠀⡼⠀⢸⠁⢰⠃⢀⢧⠇⢀⡇⠀⠀⠀⠀│
      │⠀⠀⠀⠀⢀⠎⢀⡞⠀⡼⠁⡰⠃⣠⠃⢠⠏⢀⡞⢀⡿⠀⡼⠁⣰⠃⢠⠇⢀⠎⢀⡞⠀⢀⡜⠀⠀⠀⠀⠀│
      │⠀⠀⠀⠀⢩⣖⣯⣔⣮⣴⣾⣥⣺⣥⣲⣥⣖⣯⣔⣯⣷⣾⣵⣺⣥⣲⣥⣶⣯⣖⣋⣤⣔⣯⣔⣎⠀⠀⠀⠀│
      │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
   -2 │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
      └────────────────────────────────────────┘
      ⠀-4⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀4⠀
```
"""
function ground_track(
    orbp::OrbitPropagator;
    add_nans::Bool = true,
    duration::Number = 86400,
    initial_time::Number = 0,
    f_eci_to_ecef::Function = _ground_track_default_eci_to_ecef,
    step::Union{Nothing, Number} = nothing,
    track_types::Symbol = :all
)

    # == Prepare the Inputs ================================================================

    epoch = Propagators.epoch(orbp)

    if isnothing(step)
        # If the user did not specify the step, let's obtain the step roughly equivalent to
        # 5° in the mean anomaly.
        orb = Propagators.mean_elements(orbp)

        # If the propagator cannot return the mean elements, we will compute the orbit
        # period by converting the osculating elements to Keplerian elements.
        if isnothing(orb)
            r_i, v_i = Propagators.propagate(orbp, 0)
            orb = rv_to_kepler(r_i, v_i, epoch)
        end

        ang_speed = orbital_angular_velocity(orb)

        step = deg2rad(1) / ang_speed
    end

    # == Compute the Ground Track ==========================================================

    # Time vector to compute the ground track.
    vt = float(initial_time):float(step):float(initial_time + duration)

    # Allocate the output vector.
    gt = NTuple{2, Float64}[]
    sizehint!(gt, length(vt))

    # We need to compute the first point previously to analyze if we are in an ascending or
    # descending path.
    t = first(vt)
    r_i, ~ = Propagators.propagate!(orbp, t)
    r_ecef = f_eci_to_ecef(r_i, epoch + t / 86400)
    lat_k_1, lon_k_1, ~ = ecef_to_geodetic(r_ecef)

    for t in vt[2:end]
        # Propagate the orbit and convert to ECEF.
        r_i, ~ = Propagators.propagate!(orbp, t)
        r_ecef = f_eci_to_ecef(r_i, epoch + t / 86400)

        # Compute the latitude and longitude.
        lat_k, lon_k, ~ = ecef_to_geodetic(r_ecef)

        # Check if this is an ascending or descending path.
        current_type = (lat_k - lat_k_1) > 0 ? :ascending : :descending

        # Check if we need to add this point to the groun track vector.
        if (track_types != current_type) && (track_types != :all)
            lat_k_1 = lat_k
            lon_k_1 = lon_k
            continue
        end

        isempty(gt) && push!(gt, (lat_k_1, lon_k_1))

        # Check if we need to add NaNs to improve plotting.
        if add_nans
            gt_k_1 = last(gt)
            Δlat = lat_k - first(gt_k_1)
            Δlon = lon_k - last(gt_k_1)

            ((abs(Δlat) > π / 2) || (abs(Δlon) > π)) && push!(gt, (NaN, NaN))
        end

        push!(gt, (lat_k, lon_k))
        lat_k_1 = lat_k
        lon_k_1 = lon_k
    end

    return gt
end

"""
    ground_track_inclination(a::Number, e::Number, i::Number; kwargs...) -> T
    ground_track_inclination(orb::Orbit{Tepoch, T}); kwargs...) where {Tepoch <: Number, T <: Number} -> T

Compute the ground track inclination at the Equator [rad] in an orbit with semi-major axis
`a` [m], eccentricity `e` [ ], and inclination `i` [rad]. The orbit can also be specified by
`orb` (see `Orbit`).

!!! note
    The output type `T` in the first signature is obtained by promoting the inputs to a
    float type.

!!! warning
    The algorithm here assumes a small orbit eccentricity.

# Keywords

- `perturbation::Symbol`: Symbol to select the perturbation terms that will be used. It can
    be `:J0`, `:J2`, or `:J4`.
    (**Default**: `:J2`)
- `m0::Number`: Standard gravitational parameter for Earth [m³ / s²].
    (**Default**: `GM_EARTH`)
- `J2::Number`: J₂ perturbation term.
    (**Default**: `EGM_2008_J2`)
- `J4::Number`: J₄ perturbation term.
    (**Default**: `EGM_2008_J4`)
- `R0::Number`: Earth's equatorial radius [m].
    (**Default**: `EARTH_EQUATORIAL_RADIUS`)
- `we::Number`: Earth's angular speed [rad / s].
    (**Default**: `EARTH_ANGULAR_SPEED`)

# Extended Help

We define the ground track inclination as the angle that the ground track has with respect
to the Equator. This information is important to compute, for example, the required swath
for a remote sensing satellite to cover the entire Earth.

The ground track inclination `i_gt` is given by:

```
            ┌                       ┐
            │      ω_s ⋅ sin i      │
i_gt = atan │ ───────────────────── │
            │ ω_s ⋅ cos i - ω_e + Ω̇ │
            └                       ┘
```

where `ω_s = n + ω̇` is the satellite angular velocity, `n` is the perturbed mean motion,
`i` is the orbit inclination, `ω` is the orbit argument of perigee, `Ω` is the orbit right
ascension of the ascending node, and `ω_e` is the Earth's angular rate.

Formally, we should use the satellite instantaneous angular speed at the Equator instead of
the mean angular speed `ω_s`. However, given the perturbations caused by the Earth's
gravitational potential, the former is not simple to compute. This calculation would
required to implement an orbit propagator here. Thus, we simplify it by assuming that the
orbit eccentricity is small. This assumption is reasonable given the missions that would
benefit from the computation of the ground track inclination. In this case, the orbit
angular speed is almost constant and equal to `ω_s`.

## Examples

```julia-repl
julia> using SatelliteAnalysis

julia> ground_track_inclination(7130.982e3, 0.00111, 98.410 |> deg2rad) |> rad2deg
102.30052101661899

julia> jd₀ = date_to_jd(2021, 1, 1)
2.4592155e6

julia> orb = KeplerianElements(
           jd₀,
           7130.982e3,
           0.001111,
           98.410 |> deg2rad,
           ltdn_to_raan(10.5, jd₀),
           π / 2,
           0
       )
KeplerianElements{Float64, Float64}:
           Epoch :    2.45922e6 (2021-01-01T00:00:00)
 Semi-major axis : 7130.98     km
    Eccentricity :    0.001111
     Inclination :   98.41     °
            RAAN :   78.4021   °
 Arg. of Perigee :   90.0      °
    True Anomaly :    0.0      °

julia> ground_track_inclination(orb) |> rad2deg
102.30052101658998
```
"""
function ground_track_inclination(
    a::T1,
    e::T2,
    i::T3;
    J2::Number = EGM_2008_J2,
    J4::Number = EGM_2008_J4,
    R0::Number = EARTH_EQUATORIAL_RADIUS,
    m0::Number = GM_EARTH,
    perturbation::Symbol = :J2,
    we::Number = EARTH_ANGULAR_SPEED
) where {T1 <: Number, T2 <: Number, T3 <: Number}
    T   = float(promote_type(T1, T2, T3))
    μ   = T(m0)
    ω_e = T(we)

    # Satellite mean angular velocity [rad / s].
    ω_s = orbital_angular_velocity(
        a,
        e,
        i;
        perturbation = perturbation,
        m0 = m0,
        R0 = R0,
        J2 = J2,
        J4 = J4
    )

    # The ground track inclination depends on the RAAN time derivative [rad / s].
    ∂Ω_∂t = raan_time_derivative(
        a,
        e,
        i;
        perturbation = perturbation,
        J2 = J2,
        J4 = J4,
        m0 = m0,
        R0 = R0
    )

    # Finally, we can compute the ground track inclination at Equator.
    sin_i, cos_i = sincos(T(i))
    i_t = atan(ω_s * sin_i, ω_s * cos_i - ω_e + ∂Ω_∂t)

    return i_t
end

function ground_track_inclination(orb::Orbit; kwargs...)
    # Convert first to Keplerian elements.
    k = convert(KeplerianElements, orb)
    return ground_track_inclination(k.a, k.e, k.i; kwargs...)
end

############################################################################################
#                                    Private Functions                                     #
############################################################################################

# Default function to convert `r_i` from the ECI reference frame to ECEF reference frame at
# the instant `jd`.
function _ground_track_default_eci_to_ecef(r_eci::AbstractVector, jd::Number)
    D_ecef_eci = r_eci_to_ecef(TEME(), PEF(), jd)
    return D_ecef_eci * r_eci
end
