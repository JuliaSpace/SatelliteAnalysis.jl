# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#   Functions to compute the accesses and gaps between the satellite and ground facilities.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export ground_facility_accesses, ground_facility_gaps

"""
    ground_facility_accesses(orbp, [(WGS84)], Δt, eci, ecef, vargs...; kwargs...) -> DataFrame

Compute the accesses of a satellite with orbit propagator `orbp` (see `Propagators.init`) to
the ground facilities defined in the vector `vgs_r_e`. The analysis interval begins in the
propagator epoch plus `t_0` and lasts for `Δt` [s].

The ground facilities are specified using a vector of tuples with three numbers:

    Tuple{T1, T2, T3} where {T1 <: Number, T2 <: Number, T3 <: Number}

containing the WGS84 position of each ground facility `[(WGS84)]`:

    (latitude [rad], longitude [rad], altitude [m])

Those geodetic information are transformed to an ECEF vector using the function
`geodetic_to_ecef`.

# Arguments

- `eci`: Earth-Centered Inertial frame in which the state vector of the propagator is
    represented.
- `ecef`: Earth-Centered, Earth-fixed frame to be used for the analysis. It must be the same
    frame used to compute the ground facility position.
- `vargs...`: List of additional arguments to be passed to the function `r_eci_to_ecef` when
    converting the ECI frame to the ECEF.

# Keywords

- `minimum_elevation::Number`: Minimum elevation angle for communication between the
    satellite and the ground facilities [rad]. (**Default** = 10°)
- `reduction::Function`: A function that receives a boolean vector with the visibility
    between the satellite and each ground facility. It must return a boolean value
    indicating if the access must be computed or not. This is useful to merge access time
    between two or more facilities. (**Default** = `v -> |(v...)` *i.e.* compute the access
    if at least one ground facilities is visible)
- `step::Number`: The step [s] used to propagate the orbit. Notice that we perform a cross
    tuning to accurately obtain the access time. However, if an access is lower than the
    step, it can be neglected. (**Default** = 60)
- `t_0::Number`: Initial time of the analysis after the propagator epoch [s].
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
"""
function ground_facility_accesses(
    orbp::OrbitPropagator,
    gs_wgs84::Tuple{T1, T2, T3},
    vargs::Vararg{Any, N};
    kwargs...
) where {N, T1 <: Number, T2 <: Number, T3 <: Number}
    return ground_facility_accesses(orbp, [gs_wgs84], vargs...; kwargs...)
end

function ground_facility_accesses(
    orbp::OrbitPropagator,
    vgs_wgs84::AbstractVector{T},
    Δt::Number,
    eci::Union{T_ECIs, T_ECIs_IAU_2006},
    ecef::Union{T_ECEFs, T_ECEFs_IAU_2006},
    vargs::Vararg{Any, N};
    minimum_elevation::Number = 10 |> deg2rad,
    reduction::Function = v -> |(v...),
    step::Number = 60,
    t_0::Number = 0,
    unit::Symbol = :s
) where {N, T<:Tuple{T1, T2, T3} where {T1<:Number, T2<:Number, T3<:Number}}

    # Time vector of the analysis.
    t = float(t_0):float(step):float(Δt + t_0)

    # Get the epoch of the propagator.
    jd₀ = Propagators.epoch(orbp)

    # Convert the ground facilities positions to an ECEF vector to save
    # computational burden.
    vgs_r_e = [geodetic_to_ecef(gs_wgs84...) for gs_wgs84 in vgs_wgs84]

    # Vector that will contain the accesses.
    accesses = NTuple{2, DateTime}[]

    # State to help the computation.
    state = :initial

    # Pre-allocate the visibility vector to avoid a huge number of allocation.
    visibility = zeros(Bool, length(vgs_r_e))

    # Lambda function to check the reduced visibility.
    function f(t)::Bool
        r_i, v_i  = Propagators.propagate!(orbp, t)
        r_e       = r_eci_to_ecef(DCM, eci, ecef, jd₀ + t / 86400, vargs...) * r_i

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
                access_beg = jd_to_date(DateTime, jd₀) + Dates.Second(t_0)
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

    # If the analysis finished during an access, then just add the end of the
    # interval as the end of the access.
    if state == :visible
        access_end = jd_to_date(DateTime, jd₀ + (Δt + t_0) / 86400)
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
    ground_facility_gaps(T, orbp, args...; t_0::Number = 0, kwargs...) -> [Tuple, DataFrame]

Compute the gaps between the accesses of ground facilities. The arguments and keywords are
the same as the ones used in the function [`ground_facility_accesses`](@ref).

Notice that the gap analysis starts in the orbit propagator epoch plus `t_0` and lasts for
`Δt` [s].

# Returns

- `DataFrame`: The function returns a `DataFrame` with three columns:
    - `access_beginning`: Time of the access beginning [UTC] encoded using `DateTime`.
    - `access_end`: Time of the access end [UTC] encoded using `DateTime`.
    - `duration`: Duration of the access [s].
    The unit of the column `duration` is stored in the `DataFrame` using metadata.
"""
function ground_facility_gaps(
    orbp::OrbitPropagator,
    gs_wgs84::Tuple{T1, T2, T3},
    vargs::Vararg{Any, N};
    kwargs...
) where {N, T1 <: Number, T2 <: Number, T3 <: Number}
    return ground_facility_gaps(orbp, [gs_wgs84], vargs...; kwargs...)
end

function ground_facility_gaps(
    orbp,
    vgs_wgs84::AbstractVector{T},
    Δt::Number,
    eci::Union{T_ECIs, T_ECIs_IAU_2006},
    ecef::Union{T_ECEFs, T_ECEFs_IAU_2006},
    vargs::Vararg{Any, N};
    minimum_elevation::Number = 10 |> deg2rad,
    reduction::Function = v -> |(v...),
    step::Number = 60,
    t_0::Number = 0,
    unit::Symbol = :s
) where {N, T<:Tuple{T1, T2, T3} where {T1<:Number, T2<:Number, T3<:Number}}

    # Get the epoch of the propagator.
    jd₀ = Propagators.epoch(orbp)
    dt₀ = jd_to_date(DateTime, jd₀) + Dates.Second(t_0)

    # Compute the list of ground facility accesses.
    dfa = ground_facility_accesses(
        orbp,
        vgs_wgs84,
        Δt,
        eci,
        ecef,
        vargs...;
        minimum_elevation,
        reduction,
        step,
        t_0
    )

    # Compute the last propagation instant.
    jd₁ = jd₀ + (Δt + t_0) / 86400
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
