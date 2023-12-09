## Description #############################################################################
#
# Functions to compute the accesses and gaps between the satellite and ground facilities.
#
############################################################################################

export ground_facility_accesses, ground_facility_gaps

"""
    ground_facility_accesses(orbp, [(WGS84)]; kwargs...) -> DataFrame

Compute the accesses of a satellite with orbit propagator `orbp` (see `Propagators.init`) to
the ground facilities defined in the vector `vgs_r_e`. The analysis interval begins in the
propagator epoch plus `initial_time` and lasts for `duration` [s], where both are keywords.

The ground facilities are specified using a vector of tuples with three numbers:

    Tuple{T1, T2, T3} where {T1 <: Number, T2 <: Number, T3 <: Number}

containing the WGS84 position of each ground facility `[(WGS84)]`:

    (latitude [rad], longitude [rad], altitude [m])

Those geodetic information are transformed to an ECEF vector using the function
`geodetic_to_ecef`.

# Keywords

- `duration::Number`: Duration of the analysis [s].
    (**Default** = 86400)
- `f_eci_to_ecef::Function`: Function to convert the orbit propagator position represented
    in the Earth-centered inertial (ECI) reference frame to the Earth-centered, Earth-fixed
    (ECEF) reference frame. The signature must be

    ```julia
    f_eci_to_ecef(r_i::AbstractVector, jd::Number) -> AbstractVector
    ```

    and it must return the position vector `r_i` represented in the ECEF at the instant `jd`
    [Julian Day]. By default, we use TEME as the ECI and PEF as the ECEF.
    (**Default**: `_ground_facility_default_eci_to_ecef`)
- `initial_time::Number`: Initial time of the analysis after the propagator epoch [s].
    (**Default** = 0)
- `minimum_elevation::Number`: Minimum elevation angle for communication between the
    satellite and the ground facilities [rad].
    (**Default** = 10°)
- `reduction::Function`: A function that receives a boolean vector with the visibility
    between the satellite and each ground facility. It must return a boolean value
    indicating if the access must be computed or not. This is useful to merge access time
    between two or more facilities.
    (**Default** = `v -> |(v...)` *i.e.* compute the access if at least one ground
    facilities is visible)
- `step::Number`: The step [s] used to propagate the orbit. Notice that we perform a cross
    tuning to accurately obtain the access time. However, if an access is lower than the
    step, it can be neglected.
    (**Default** = 60)
- `unit::Symbol`: Select the unit in which the duration will be computed. The possible
    values are:
    - `:s` for seconds (**Default**);
    - `:m` for minutes; or
    - `:h` for hours.

# Returns

- `DataFrame`: The function returns a `DataFrame` with three columns:
    - `access_beginning`: Time of the access beginning [UTC] encoded using `DateTime`.
    - `access_end`: Time of the access end [UTC] encoded using `DateTime`.
    - `duration`: Duration of the access [s].
    The unit of the column `duration` is stored in the `DataFrame` using metadata.

# Extended Help

## Examples

```julia-repl
julia> using SatelliteAnalysis

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

julia> ground_facility_accesses(orbp, (0, 0, 0))
2×3 DataFrame
 Row │ access_beginning         access_end               duration 
     │ DateTime                 DateTime                 Float64  
─────┼────────────────────────────────────────────────────────────
   1 │ 2024-01-01T10:20:03.136  2024-01-01T10:30:02.971   599.835
   2 │ 2024-01-01T22:49:55.910  2024-01-01T22:59:23.470   567.56

julia> ground_facility_accesses(orbp, (0, 0, 0); unit = :m)
2×3 DataFrame
 Row │ access_beginning         access_end               duration 
     │ DateTime                 DateTime                 Float64  
─────┼────────────────────────────────────────────────────────────
   1 │ 2024-01-01T10:20:03.136  2024-01-01T10:30:02.971   9.99725
   2 │ 2024-01-01T22:49:55.910  2024-01-01T22:59:23.470   9.45933
```
"""
function ground_facility_accesses(
    orbp::OrbitPropagator,
    gf_wgs84::Tuple{T1, T2, T3};
    kwargs...
) where {T1 <: Number, T2 <: Number, T3 <: Number}
    return ground_facility_accesses(orbp, [gf_wgs84]; kwargs...)
end

function ground_facility_accesses(
    orbp::OrbitPropagator,
    vgf_wgs84::AbstractVector{T};
    duration::Number = 86400,
    f_eci_to_ecef::Function = _ground_facilities_default_eci_to_ecef,
    initial_time::Number = 0,
    minimum_elevation::Number = 10 |> deg2rad,
    reduction::Function = v -> |(v...),
    step::Number = 60,
    unit::Symbol = :s
) where {T<:Tuple{T1, T2, T3} where {T1<:Number, T2<:Number, T3<:Number}}

    # Time vector of the analysis.
    t = float(initial_time):float(step):float(initial_time + duration)

    # Get the epoch of the propagator.
    jd₀ = Propagators.epoch(orbp)

    # Convert the ground facilities positions to an ECEF vector to save
    # computational burden.
    vgs_r_e = [geodetic_to_ecef(gf_wgs84...) for gf_wgs84 in vgf_wgs84]

    # Vector that will contain the accesses.
    accesses = NTuple{2, DateTime}[]

    # State to help the computation.
    state = :initial

    # Pre-allocate the visibility vector to avoid a huge number of allocation.
    visibility = zeros(Bool, length(vgs_r_e))

    # Lambda function to check the reduced visibility.
    function f(t)::Bool
        r_i, ~ = Propagators.propagate!(orbp, t)
        r_e = f_eci_to_ecef(r_i, jd₀ + t / 86400)

        @inbounds for i in eachindex(visibility)
            visibility[i] = is_ground_facility_visible(r_e, vgs_r_e[i], minimum_elevation)
        end

        return reduction(visibility)
    end

    access_beg  = DateTime(now())
    access_end  = DateTime(now())
    vaccess_beg = DateTime[]
    vaccess_end = DateTime[]

    for k in t
        # Check the initial state of the reduced visibility.
        visible = f(k)

        # Handle the initial case.
        if state == :initial
            if visible
                access_beg = jd_to_date(DateTime, jd₀) + Dates.Second(round(Int, inital_time))
                state = :visible
            else
                state = :not_visible
            end

        # Handle transitions.
        elseif (state == :not_visible) && visible
            # Refine to find the edge.
            k₀ = k - step
            k₁ = k
            kc = find_crossing(f, k₀, k₁, false, true)

            state = :visible
            access_beg = jd_to_date(DateTime, jd₀ + kc / 86400)

        elseif (state == :visible) && !visible
            # Refine to find the edge.
            k₀ = k - step
            k₁ = k
            kc = find_crossing(f, k₀, k₁, true, false)

            state = :not_visible
            access_end = jd_to_date(DateTime, jd₀ + kc / 86400)

            push!(vaccess_beg, access_beg)
            push!(vaccess_end, access_end)
        end
    end

    # If the analysis finished during an access, then just add the end of the interval as
    # the end of the access.
    if state == :visible
        access_end = jd_to_date(DateTime, jd₀ + (initial_time + duration) / 86400)
        push!(vaccess_beg, access_beg)
        push!(vaccess_end, access_end)
    end

    # Compute the access duration and convert to the desired unit.
    duration = Dates.value.(vaccess_end .- vaccess_beg) ./ 1000

    if unit == :h
        duration ./= 3600
    elseif unit == :m
        duration ./= 60
    else
        # If the symbol is not known, we must use seconds.
        unit = :s
    end

    # Create the DataFrame and write the metadata.
    df = DataFrame(
        :access_beginning => vaccess_beg,
        :access_end       => vaccess_end,
        :duration         => duration
    )

    metadata!(df, "Description", "Accesses to the ground facilities.")
    colmetadata!(df, :duration, "Unit", unit)

    return df
end

"""
    ground_facility_gaps(orbp, args...; duration::Number = 86400, initial_time::Number = 0, kwargs...) -> DataFrame

Compute the gaps between the accesses of ground facilities. The arguments and keywords are
the same as the ones used in the function [`ground_facility_accesses`](@ref).

Notice that the gap analysis starts in the orbit propagator epoch plus `initial_time` and
lasts for `duration` [s].

# Returns

- `DataFrame`: The function returns a `DataFrame` with three columns:
    - `gap_beginning`: Time of the access beginning [UTC] encoded using `DateTime`.
    - `gap_end`: Time of the access end [UTC] encoded using `DateTime`.
    - `duration`: Duration of the access [s].
    The unit of the column `duration` is stored in the `DataFrame` using metadata.

# Extended Help

## Examples

```julia-repl
julia> using SatelliteAnalysis

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

julia> ground_facility_gaps(orbp, (0, 0, 0))
3×3 DataFrame
 Row │ gap_beginning            gap_end                  duration 
     │ DateTime                 DateTime                 Float64  
─────┼────────────────────────────────────────────────────────────
   1 │ 2024-01-01T00:00:00      2024-01-01T10:20:03.136  37203.1
   2 │ 2024-01-01T10:30:02.971  2024-01-01T22:49:55.910  44392.9
   3 │ 2024-01-01T22:59:23.470  2024-01-02T00:00:00       3636.53

julia> ground_facility_gaps(orbp, (0, 0, 0); unit = :m)
3×3 DataFrame
 Row │ gap_beginning            gap_end                  duration 
     │ DateTime                 DateTime                 Float64  
─────┼────────────────────────────────────────────────────────────
   1 │ 2024-01-01T00:00:00      2024-01-01T10:20:03.136  620.052
   2 │ 2024-01-01T10:30:02.971  2024-01-01T22:49:55.910  739.882
   3 │ 2024-01-01T22:59:23.470  2024-01-02T00:00:00       60.6088
```
"""
function ground_facility_gaps(
    orbp::OrbitPropagator,
    gf_wgs84::Tuple{T1, T2, T3};
    kwargs...
) where {T1 <: Number, T2 <: Number, T3 <: Number}
    return ground_facility_gaps(orbp, [gf_wgs84]; kwargs...)
end

function ground_facility_gaps(
    orbp,
    vgf_wgs84::AbstractVector{T};
    duration::Number = 86400,
    f_eci_to_ecef::Function = _ground_facilities_default_eci_to_ecef,
    initial_time::Number = 0,
    minimum_elevation::Number = 10 |> deg2rad,
    reduction::Function = v -> |(v...),
    step::Number = 60,
    unit::Symbol = :s
) where {T<:Tuple{T1, T2, T3} where {T1<:Number, T2<:Number, T3<:Number}}

    # Get the epoch of the propagator.
    jd₀ = Propagators.epoch(orbp)
    dt₀ = jd_to_date(DateTime, jd₀) + Dates.Second(round(Int, initial_time))

    # Compute the list of ground facility accesses.
    dfa = ground_facility_accesses(
        orbp,
        vgf_wgs84;
        duration,
        f_eci_to_ecef,
        initial_time,
        minimum_elevation,
        reduction
    )

    # Compute the last propagation instant.
    jd₁ = jd₀ + (initial_time + duration) / 86400
    dt₁ = jd_to_date(DateTime, jd₁)

    # Compute the gaps between accesses.
    vgap_beg = DateTime[]
    vgap_end = DateTime[]

    # If the number of accesses is 0, return the entire interval.
    num_rows, num_cols = size(dfa)

    if num_rows == 0
        push!(vgap_beg, dt₀)
        push!(vgap_end, dt₁)
    else
        # Check if the simulation did not start under the visibility of a ground facility.
        access_beginning = dfa.access_beginning
        access_end       = dfa.access_end

        if first(access_beginning) != dt₀
            push!(vgap_beg, dt₀)
            push!(vgap_end, access_beginning |> first)
        end

        @inbounds for k in 1:(num_rows - 1)
            push!(vgap_beg, access_end[k])
            push!(vgap_end, access_beginning[k + 1])
        end

        # Check if the simulation did not end under the visibility of a ground facility.
        if last(access_end) != dt₁
            push!(vgap_beg, last(access_end))
            push!(vgap_end, dt₁)
        end
    end

    # Compute the access duration and convert to the desired unit.
    duration = Dates.value.(vgap_end .- vgap_beg) ./ 1000

    if unit == :h
        duration ./= 3600
    elseif unit == :m
        duration ./= 60
    else
        # If the symbol is not known, we must use seconds.
        unit = :s
    end

    # Create the DataFrame and write the metadata.
    dfg = DataFrame(
        :gap_beginning => vgap_beg,
        :gap_end       => vgap_end,
        :duration      => duration
    )

    metadata!(dfg, "Description", "Gaps to the ground facilities.")
    colmetadata!(dfg, :duration, "Unit", unit)

    return dfg
end

############################################################################################
#                                    Private Functions                                     #
############################################################################################

# Default function to convert `r_i` from the ECI reference frame to ECEF reference frame at
# the instant `jd`.
function _ground_facilities_default_eci_to_ecef(r_eci::AbstractVector, jd::Number)
    D_ecef_eci = r_eci_to_ecef(TEME(), PEF(), jd)
    return D_ecef_eci * r_eci
end
