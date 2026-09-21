## Description #############################################################################
#
# Functions to compute the accesses and gaps between the satellite and ground facilities.
#
############################################################################################

export ground_facility_accesses, ground_facility_gaps

"""
    ground_facility_accesses(orbp, [(WGS84)]; kwargs...) -> DataFrame

Compute the accesses of a satellite with orbit propagator `orbp` (see `Propagators.init`) to
the ground facilities defined in the vector `[(WGS84)]`. The analysis interval begins in the
propagator epoch plus `initial_time` and lasts for `duration` [s], where both are keywords.

The ground facilities are specified using a vector of tuples with three numbers:

    Tuple{T1, T2, T3} where {T1 <: Number, T2 <: Number, T3 <: Number}

containing the WGS84 position of each ground facility `[(WGS84)]`:

    (latitude [rad], longitude [rad], altitude [m])

Those geodetic information are transformed to an ECEF vector using the function
`geodetic_to_ecef`.

!!! warning

    This function computes the accesses using multiple threads. Hence, the function
    `f_eci_to_ecef` must be thread safe.

# Keywords

- `duration::Number`: Duration of the analysis [s].
    (**Default** = 86400)
- `f_eci_to_ecef::Function`: Function to convert the orbit propagator position represented
    in the Earth-centered inertial (ECI) reference frame to the Earth-centered, Earth-fixed
    (ECEF) reference frame. The signature must be

    ```
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
- `num_chunks::Number`: Number of chunks the algorithm will divide the time vector to
    compute the accesses.
    (**Default** = `Threads.nthreads()`)
- `reduction::Function`: A function that receives a boolean vector with the visibility
    between the satellite and each ground facility. It must return a boolean value
    indicating if the access must be computed or not. This is useful to merge access time
    between two or more facilities.
    (**Default** = `any` *i.e.* compute the access if at least one ground
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

## Throws

- `ArgumentError`: If `unit` is not `:s`, `:m`, or `:h`.

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
    orbp::OrbitPropagator, gf_wgs84::Tuple{T1, T2, T3}; kwargs...
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
    num_chunks::Integer = Threads.nthreads(),
    reduction::R = any,
    step::Number = 60,
    unit::Symbol = :s,
) where {
    T <: Tuple{T1, T2, T3} where {T1 <: Number, T2 <: Number, T3 <: Number}, R <: Function
}

    # Factor to convert the time from seconds to the selected unit. Notice that this
    # function also validates the input.
    time_factor = _time_unit_factor(unit)

    # Time vector of the analysis.
    vt = float(initial_time):float(step):float(initial_time + duration)

    # Create the chunks with the time to be computed by each thread.
    # There is no useful work in creating more tasks than propagation instants. Clamp the
    # requested count also so extreme values cannot create invalid chunk indices.
    num_chunks = min(max(num_chunks, 1), max(length(vt), 1))
    vt_chunks = _gf_access_time_vector_partition(vt, num_chunks) |> collect

    @debug begin
        dtf = dateformat"yyyy-mm-ddTHH:MM:SS.sss"

        # Create a vector with the instants computed by each chunk.
        dt₀ = julian2datetime(Propagators.epoch(orbp))

        str_vt_chunks = [
            Dates.format(
                dt₀ + Dates.Microsecond(round(Int, vt_chunks[k][begin] * 1e6)), dtf
            ) *
            " -- " *
            Dates.format(dt₀ + Dates.Microsecond(round(Int, vt_chunks[k][end] * 1e6)), dtf)
            for k in eachindex(vt_chunks)
        ]

        debug_msg = "Computing ground facility accesses using $(num_chunks) chunks:\n\n"

        for k in eachindex(str_vt_chunks)
            str_k = lpad(string(k), floor(Int, log10(num_chunks)) + 1)
            debug_msg *= "Chunk $str_k: $(str_vt_chunks[k])\n"
        end

        debug_msg
    end

    # Create the tasks to compute by each thread.
    tasks = map(eachindex(vt_chunks)) do c
        chunk_vt = vt_chunks[c]

        # A propagation modified the propagator structure. Hence, we need to copy the
        # structure for each thread to avoid racing conditions.
        chunk_orbp = c == 1 ? orbp : deepcopy(orbp)

        Threads.@spawn begin
            _ground_facility_access_chunk(
                chunk_orbp,
                chunk_vt,
                vgf_wgs84;
                f_eci_to_ecef     = f_eci_to_ecef,
                minimum_elevation = minimum_elevation,
                reduction         = reduction,
            )
        end
    end

    # Fetch the accesses computed in each chunk. The access beginnings and ends are
    # represented in seconds since the propagator epoch.
    chunk_accesses = fetch.(tasks)

    # == Merge the Accesses ================================================================

    Tt = eltype(vt)

    vaccess_beg_s = Tt[]
    vaccess_end_s = Tt[]

    # Two consecutive chunks share the instant in their boundary. Hence, if a chunk ended
    # during an access, and the next one started during an access, we have a single access
    # divided into two chunks that must be merged.
    previous_ended_visible = false

    for (vchunk_beg_s, vchunk_end_s, started_visible, ended_visible) in chunk_accesses
        k₀ = firstindex(vchunk_beg_s)

        if previous_ended_visible && started_visible && !isempty(vaccess_end_s)
            vaccess_end_s[end] = vchunk_end_s[k₀]
            k₀ += 1
        end

        append!(vaccess_beg_s, @view(vchunk_beg_s[k₀:end]))
        append!(vaccess_end_s, @view(vchunk_end_s[k₀:end]))

        previous_ended_visible = ended_visible
    end

    # Convert the access beginnings and ends to `DateTime`. Notice that all the conversions
    # must use the same function to avoid inconsistencies caused by rounding.
    dt₀ = julian2datetime(Propagators.epoch(orbp))

    vaccess_beg = _gf_seconds_to_datetime.(dt₀, vaccess_beg_s)
    vaccess_end = _gf_seconds_to_datetime.(dt₀, vaccess_end_s)

    # Compute the access duration and convert to the desired unit.
    vaccess_duration = Dates.value.(vaccess_end .- vaccess_beg) ./ 1000

    vaccess_duration .*= time_factor

    # Create the DataFrame and write the metadata.
    df = DataFrame(
        :access_beginning => vaccess_beg,
        :access_end       => vaccess_end,
        :duration         => vaccess_duration,
    )

    # The style `:note` makes the metadata propagate through DataFrame transformations.
    metadata!(df, "Description", "Accesses to the ground facilities."; style = :note)
    colmetadata!(df, :duration, "Unit", unit; style = :note)

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
    - `gap_beginning`: Time of the gap beginning [UTC] encoded using `DateTime`.
    - `gap_end`: Time of the gap end [UTC] encoded using `DateTime`.
    - `duration`: Duration of the gap [s].
    The unit of the column `duration` is stored in the `DataFrame` using metadata.

# Extended Help

## Throws

- `ArgumentError`: If `unit` is not `:s`, `:m`, or `:h`.

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
    orbp::OrbitPropagator, gf_wgs84::Tuple{T1, T2, T3}; kwargs...
) where {T1 <: Number, T2 <: Number, T3 <: Number}
    return ground_facility_gaps(orbp, [gf_wgs84]; kwargs...)
end

function ground_facility_gaps(
    orbp::OrbitPropagator,
    vgf_wgs84::AbstractVector{T};
    duration::Number = 86400,
    initial_time::Number = 0,
    unit::Symbol = :s,
    kwargs...,
) where {T <: Tuple{T1, T2, T3} where {T1 <: Number, T2 <: Number, T3 <: Number}}

    # Factor to convert the time from seconds to the selected unit. Notice that this
    # function also validates the input.
    time_factor = _time_unit_factor(unit)

    # Compute the beginning and the end of the analysis. Notice that we must use the same
    # function used to convert the access instants to avoid inconsistencies caused by
    # rounding.
    dt_epoch = julian2datetime(Propagators.epoch(orbp))

    dt₀ = _gf_seconds_to_datetime(dt_epoch, initial_time)
    dt₁ = _gf_seconds_to_datetime(dt_epoch, initial_time + duration)

    # Compute the list of ground facility accesses. All the other keywords are forwarded to
    # the function that computes the accesses.
    dfa = ground_facility_accesses(orbp, vgf_wgs84; duration, initial_time, kwargs...)

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

    duration .*= time_factor

    # Create the DataFrame and write the metadata.
    dfg = DataFrame(:gap_beginning => vgap_beg, :gap_end       => vgap_end, :duration      => duration)

    # The style `:note` makes the metadata propagate through DataFrame transformations.
    metadata!(dfg, "Description", "Gaps to the ground facilities."; style = :note)
    colmetadata!(dfg, :duration, "Unit", unit; style = :note)

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

# Return a generator that contains the time partition of vector `vt` into `np` parts.
#
# This code was adapted from the one in the blog post:
#
#   https://blog.glcs.io/parallel-processing
function _gf_access_time_vector_partition(vt::AbstractVector, np::Integer)
    len_vt = length(vt)

    # An empty time vector has no partitions. This also keeps the helper well-defined for
    # callers requesting an extreme number of chunks.
    len_vt == 0 && return Base.Generator(identity, 1:0)

    len, rem = divrem(len_vt, np)

    # Treat the case in which we want more partitions than the number of elements.
    if len == 0
        np = len_vt
    end

    Base.Generator(1:np) do p
        i₀ = firstindex(vt) + (p - 1) * len
        i₁ = p < np ? i₀ + len : i₀ + len - 1

        i₀ += p <= rem ? p - 1 : rem
        i₁ += p <= rem ? p : rem

        chunk = vt[i₀:i₁]
        return chunk
    end
end

# Compute the ground facility access for a specific time chunk.
function _ground_facility_access_chunk(
    orbp::OrbitPropagator,
    vt::StepRangeLen,
    vgf_wgs84::AbstractVector{T};
    f_eci_to_ecef::Function = _ground_facilities_default_eci_to_ecef,
    minimum_elevation::Number = 10 |> deg2rad,
    reduction::R = any,
) where {
    T <: Tuple{T1, T2, T3} where {T1 <: Number, T2 <: Number, T3 <: Number}, R <: Function
}

    # Get the epoch of the propagator.
    jd₀ = Propagators.epoch(orbp)

    # Get the step in the time chunk.
    Δt = step(vt)

    # State to help the computation.
    state = :initial

    # Pre-allocate the visibility vector for custom reductions only. The built-in
    # reductions evaluate visibility directly and do not need one.
    visibility =
        (reduction === any || reduction === all) ? nothing : zeros(Bool, length(vgf_wgs84))

    # Lambda function to check the reduced visibility.
    function f(t)::Bool
        r_i, ~ = Propagators.propagate!(orbp, t)
        r_e = f_eci_to_ecef(r_i, jd₀ + t / 86400)

        # `any` and `all` are the common reductions. Evaluate them directly so that each
        # propagation instant does not require constructing/filling a Bool vector. Custom
        # reductions retain the historical vector-based API.
        if reduction === any
            @inbounds for gf in vgf_wgs84
                is_ground_facility_visible(r_e, gf..., minimum_elevation) && return true
            end
            return false
        end

        if reduction === all
            @inbounds for gf in vgf_wgs84
                is_ground_facility_visible(r_e, gf..., minimum_elevation) || return false
            end
            return true
        end

        @inbounds for i in eachindex(visibility)
            visibility[i] = is_ground_facility_visible(
                r_e, vgf_wgs84[i]..., minimum_elevation
            )
        end

        return reduction(visibility)
    end

    # Beginning of the current access and vectors with the beginning and the end of all the
    # accesses in this chunk. All the values are in seconds since the propagator epoch.
    Tt = eltype(vt)

    access_beg_s  = first(vt)
    vaccess_beg_s = Tt[]
    vaccess_end_s = Tt[]

    # Flag indicating whether the reduced visibility was `true` at the chunk beginning.
    started_visible = false

    for k in vt
        # Check the initial state of the reduced visibility.
        visible = f(k)

        # Handle the initial case.
        if state == :initial
            if visible
                access_beg_s    = k
                started_visible = true
                state           = :visible
            else
                state = :not_visible
            end

            # Handle transitions.
        elseif (state == :not_visible) && visible
            # Refine to find the edge.
            access_beg_s = find_crossing(f, k - Δt, k, false, true)
            state = :visible

        elseif (state == :visible) && !visible
            # Refine to find the edge.
            access_end_s = find_crossing(f, k - Δt, k, true, false)
            state = :not_visible

            push!(vaccess_beg_s, access_beg_s)
            push!(vaccess_end_s, access_end_s)
        end
    end

    # If the analysis finished during an access, then just add the end of the interval as
    # the end of the access.
    ended_visible = state == :visible

    if ended_visible
        push!(vaccess_beg_s, access_beg_s)
        push!(vaccess_end_s, last(vt))
    end

    return vaccess_beg_s, vaccess_end_s, started_visible, ended_visible
end

# Convert the instant `t` [s], measured from the epoch `dt₀`, to `DateTime` by rounding it to
# the nearest millisecond.
function _gf_seconds_to_datetime(dt₀::DateTime, t::Number)
    return dt₀ + Dates.Millisecond(round(Int, 1000t))
end
