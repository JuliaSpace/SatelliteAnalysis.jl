# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==============================================================================
#
#   Functions to compute the accesses and gaps between the satellite and ground
#   facilities.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export ground_facility_accesses, ground_facility_gaps

"""
    ground_facility_accesses(T, orbp, [(WGS84)], Δt, eci, ecef, vargs...; kwargs...)

Compute the accesses of a satellite with orbit propagator `orbp` (see
[`init_orbit_propagator`](@ref)) to the ground facilities defined in the vector
`vgs_r_e`. The analysis interval begins in the propagator epoch plus `t_0` and
lasts for `Δt` [s].

The ground facilities are specified using a vector of tuples with three numbers
(`Tuple{T1, T2, T3} where {T1 <: Number, T2 <: Number, T3 <:
Number}`) containing the WGS84 position of each ground facility `[(WGS84)]`:

    (latitude [rad], longitude [rad], altitude [m])

Those geodetic information are transformed to an ECEF vector using the function
[`geodetic_to_ecef`](@ref).

The return depends on the parameter `T`:

- `T = Tuple`: The function return a vector of tuples. Each element represent an
    access between the satellite and the ground facility. The beginning of the
    access [UTC] is in the first element in the tuple whereas the end is in the
    second element. Both are represented using [`DateTime`](@ref).
- `T = DataFrame`: The function returns a `DataFrame` with three columns:
    - `access_beginning`: Time of the access beginning [UTC] encoded using
        [`DateTime`](@ref).
    - `access_end`: Time of the access end [UTC] encoded using
        [`DateTime`](@ref).
    - `duration`: Duration of the access [s].

# Arguments

- `eci`: Earth-Centered Inertial frame in which the state vector of the
    propagator is represented.
- `ecef`: Earth-Centered, Earth-fixed frame to be used for the analysis. It
    must be the same frame used to compute the ground facility position.
- `vargs...`: List of additional arguments to be passed to the function
    [`r_eci_to_ecef`](@ref) when converting the ECI frame to the ECEF.

# Keywords

- `θ::Number`: Minimum elevation angle for communication between the satellite
    and the ground facilities [rad]. (**Default** = 10ᵒ)
- `reduction::Function`: A function that receives a boolean vector with the
    visibility between the satellite and each ground facility. It must return a
    boolean value indicating if the access must be computed or not. This is
    useful to merge access time between two or more facilities.
    (**Default** = `v -> |(v...)` *i.e.* compute the access if at least one
    ground facilities is visible)
- `step::Number`: The step [s] used to propagate the orbit. Notice that we
    perform a cross tuning to accurately obtain the access time. However, if an
    access is lower than the step, it can be neglected. (**Default** = 60)
- `t_0::Number`: Initial time of the analysis after the propagator epoch [s].
"""
function ground_facility_accesses(
    T::DataType,
    orbp::OrbitPropagator,
    gs_wgs84::Tuple{T1, T2, T3},
    vargs...;
    kwargs...
) where {T1 <: Number, T2 <: Number, T3 <: Number}
    return ground_facility_accesses(T, orbp, [gs_wgs84], vargs...; kwargs...)
end

function ground_facility_accesses(
    ::Type{DataFrame},
    orbp::OrbitPropagator,
    vgs_wgs84::AbstractVector{T},
    Δt::Number,
    eci::Union{T_ECIs, T_ECIs_IAU_2006},
    ecef::Union{T_ECEFs, T_ECEFs_IAU_2006},
    vargs...;
    kwargs...
) where T<:Tuple{T1, T2, T3} where {T1<:Number, T2<:Number, T3<:Number}

    accesses = ground_facility_accesses(
        Tuple,
        orbp,
        vgs_wgs84,
        Δt,
        eci,
        ecef,
        vargs...;
        kwargs...
    )

    access_beginning = map(x -> first(x), accesses)
    access_end       = map(x -> last(x),  accesses)

    return DataFrame(
        :access_beginning => access_beginning,
        :access_end       => access_end,
        :duration         => Dates.value.(access_end .- access_beginning) ./ 1000
    )
end

function ground_facility_accesses(
    ::Type{Tuple},
    orbp::OrbitPropagator,
    vgs_wgs84::AbstractVector{T},
    Δt::Number,
    eci::Union{T_ECIs, T_ECIs_IAU_2006},
    ecef::Union{T_ECEFs, T_ECEFs_IAU_2006},
    vargs...;
    θ::Number = 10 |> deg2rad,
    reduction::Function = v->|(v...),
    step::Number = 60,
    t_0::Number = 0
) where T<:Tuple{T1, T2, T3} where {T1<:Number, T2<:Number, T3<:Number}

    # Time vector of the analysis.
    t = float(t_0):float(step):float(Δt + t_0)

    # Get the epoch of the propagator.
    jd₀ = get_epoch(orbp)

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
        r_i, v_i  = propagate!(orbp, t)
        r_e       = r_eci_to_ecef(DCM, eci, ecef, jd₀ + t / 86400, vargs...) * r_i

        @inbounds for i in eachindex(visibility)
            visibility[i] = is_ground_facility_visible(r_e, vgs_r_e[i], θ)
        end

        return reduction(visibility)
    end

    access_beg = DateTime(now())
    access_end = DateTime(now())

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

            push!(accesses, (access_beg, access_end))
        end
    end

    # If the analysis finished during an access, then just add the end of the
    # interval as the end of the access.
    if state == :visible
        access_end = jd_to_date(DateTime, jd₀ + (Δt + t_0) / 86400)
        push!(accesses, (access_beg, access_end))
    end

    return accesses
end

"""
    ground_facility_gaps(T, orbp, args...; t_0::Number = 0, kwargs...)

Compute the gaps between the accesses of ground facilities. The arguments and
keywords are the same as the ones used in the function
[`ground_facility_accesses`](@ref).

Notice that the gap analysis starts in the orbit propagator epoch plus `t_0` and
lasts for `Δt` [s].

The return depends on the parameter `T`:

- `T = Tuple`: The function return a vector of tuples. Each element represent an
    access gap between the satellite and the ground facility. The beginning of
    the gap [UTC] is in the first element in the tuple whereas the end is in the
    second element. Both are represented using [`DateTime`](@ref).
- `T = DataFrame`: The function returns a `DataFrame` with three columns:
    - `gap_beginning`: Time of the gap beginning [UTC] encoded using
        [`DateTime`](@ref).
    - `gap_end`: Time of the gap end [UTC] encoded using [`DateTime`](@ref).
    - `duration`: Duration of the gap [s].
"""
function ground_facility_gaps(
    T::DataType,
    orbp::OrbitPropagator,
    gs_wgs84::Tuple{T1, T2, T3},
    vargs...;
    kwargs...
) where {T1 <: Number, T2 <: Number, T3 <: Number}
    return ground_facility_gaps(T, orbp, [gs_wgs84], vargs...; kwargs...)
end

function ground_facility_gaps(
    ::Type{DataFrame},
    orbp::OrbitPropagator,
    vgs_wgs84::AbstractVector{T},
    Δt::Number,
    eci::Union{T_ECIs, T_ECIs_IAU_2006},
    ecef::Union{T_ECEFs, T_ECEFs_IAU_2006},
    vargs...;
    kwargs...
) where T<:Tuple{T1, T2, T3} where {T1<:Number, T2<:Number, T3<:Number}

    accesses = ground_facility_gaps(
        Tuple,
        orbp,
        vgs_wgs84,
        Δt,
        eci,
        ecef,
        vargs...;
        kwargs...
    )

    gap_beginning = map(x -> first(x), accesses)
    gap_end       = map(x -> last(x),  accesses)

    return DataFrame(
        :gap_beginning => gap_beginning,
        :gap_end       => gap_end,
        :duration      => Dates.value.(gap_end .- gap_beginning) ./ 1000
    )
end

function ground_facility_gaps(
    ::Type{Tuple},
    orbp,
    vgs_wgs84::AbstractVector{T},
    Δt::Number,
    eci::Union{T_ECIs, T_ECIs_IAU_2006},
    ecef::Union{T_ECEFs, T_ECEFs_IAU_2006},
    vargs...;
    θ::Number = 10 |> deg2rad,
    reduction::Function = v->|(v...),
    step::Number = 60,
    t_0::Number = 0
) where T<:Tuple{T1, T2, T3} where {T1<:Number, T2<:Number, T3<:Number}

    # Get the epoch of the propagator.
    jd₀ = get_epoch(orbp)
    dt₀ = jd_to_date(DateTime, jd₀) + Dates.Second(t_0)

    # Compute the list of ground facility accesses.
    accesses = ground_facility_accesses(
        Tuple,
        orbp,
        vgs_wgs84,
        Δt,
        eci,
        ecef,
        vargs...;
        θ,
        reduction,
        step,
        t_0
    )

    # Compute the last propagation instant.
    jd₁ = jd₀ + (Δt + t_0) / 86400
    dt₁ = jd_to_date(DateTime, jd₁)

    # Compute the gaps between accesses.
    gaps = NTuple{2, DateTime}[]

    # If the number of accessess is 0, then just return.
    num_accesses = length(accesses)

    num_accesses == 0 && return gaps

    # Check if the simulation did not start under the visibility of a ground
    # facility.
    accesses[1][1] != dt₀ && push!(gaps, (dt₀, accesses[1][1]))

    @inbounds for k in 1:(length(accesses) - 1)
        push!(gaps, (accesses[k][2], accesses[k + 1][1]))
    end

    # Check if the simulation did not end under the visibility of a ground
    # facility.
    accesses[end][2] != dt₁ && push!(gaps, (accesses[end][2], dt₁))

    return gaps
end
