## Description #############################################################################
#
# Functions to compute the satellite ground track.
#
############################################################################################

export ground_track

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

    ```julia
    f_eci_to_ecef(r_i::AbstractVector, jd::Number) -> AbstractVector
    ```

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

```julia
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
    vt = initial_time:step:duration

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

############################################################################################
#                                    Private Functions                                     #
############################################################################################

# Default function to convert `r_i` from the ECI reference frame to ECEF reference frame at
# the instant `jd`.
function _ground_track_default_eci_to_ecef(r_eci::AbstractVector, jd::Number)
    D_ecef_eci = r_eci_to_ecef(TEME(), PEF(), jd)
    return D_ecef_eci * r_eci
end
