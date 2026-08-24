## Description #############################################################################
#
# Main function to compute the decay analysis of a satellite in orbit.
#
############################################################################################

# Flag indicating whether the space index set required by the default `F107` keyword was
# already initialized in this session, avoiding re-fetching and re-parsing the remote file
# at every call. A benign race that initializes the set twice is acceptable.
const _PREDICTED_F107_INITIALIZED = Ref(false)

# Cache of the default gravity model (EGM2008), avoiding re-fetching and re-parsing the
# ICGEM file at every call. A benign race that loads the model twice is acceptable.
const _DEFAULT_GRAVITY_MODEL = Ref{Union{Nothing, AbstractGravityModel}}(nothing)

function SatelliteAnalysis.decay_analysis(
    orb::KeplerianElements;
    # Required keywords.
    satellite_mass::Number,
    satellite_mean_area::Number,
    # Optional keywords. The keywords `atmospheric_model` and `F107` accept any callable
    # object, hence they are not annotated with `::Function`.
    atmospheric_model = nothing,
    gravity_model::Union{AbstractGravityModel, Nothing} = nothing,
    num_sampling_points_per_orbit::Int = 17,
    abstol::Number = 1e-6,
    C_d::Number = 2.2,
    C_r::Number = 1.25,
    distance_unit::Symbol = :km,
    F107 = nothing,
    reltol::Number = 1e-6,
    return_solution::Bool = false,
    solver = VCABM(),
    terminate_altitude::Number = 120e3,
    tf::Number = 30 * 365.25 * 86400.0,
    time_unit::Symbol = :y,
)
    gm = if isnothing(gravity_model)
        if isnothing(_DEFAULT_GRAVITY_MODEL[])
            _DEFAULT_GRAVITY_MODEL[] =
                GravityModels.load(IcgemFile, fetch_icgem_file(:EGM2008))
        end

        _DEFAULT_GRAVITY_MODEL[]
    else
        gravity_model
    end

    # The default atmospheric model carries a mutable buffer, so a fresh instance is built
    # per call to keep the public API thread-safe.
    atmospheric_model′ = isnothing(atmospheric_model) ?
        _Nrlmsise00AtmosphericModel() :
        atmospheric_model

    F107′ = if isnothing(F107)
        # Notice that if the user calls `SpaceIndices.destroy()` after this initialization,
        # `space_index` raises a clear error asking to initialize the space indices again.
        if !_PREDICTED_F107_INITIALIZED[]
            SpaceIndices.init(SpaceIndices.SatelliteToolboxSpaceIndexSets)
            _PREDICTED_F107_INITIALIZED[] = true
        end

        jd -> space_index(Val(:F10predicted), jd)
    elseif F107 isa Number
        _ -> Float64(F107)
    else
        F107
    end

    # The keyword `gravity_model` is abstractly typed and the callables in
    # `atmospheric_model` and `F107′` have call-site dependent types. This function
    # barrier ensures they are concretely typed inside the numerical integration, avoiding
    # dynamic dispatch at every right-hand-side evaluation. The callables are passed as
    # positional arguments since keyword arguments cannot bind the type parameters that
    # force the specialization.
    return _decay_analysis(
        orb,
        gm,
        atmospheric_model′,
        F107′;
        satellite_mass                = satellite_mass,
        satellite_mean_area           = satellite_mean_area,
        num_sampling_points_per_orbit = num_sampling_points_per_orbit,
        abstol                        = abstol,
        C_d                           = C_d,
        C_r                           = C_r,
        distance_unit                 = distance_unit,
        reltol                        = reltol,
        return_solution               = return_solution,
        solver                        = solver,
        terminate_altitude            = terminate_altitude,
        tf                            = tf,
        time_unit                     = time_unit
    )
end

############################################################################################
#                                    Private Functions                                     #
############################################################################################

function _decay_analysis(
    orb::KeplerianElements,
    gm::AbstractGravityModel,
    atmospheric_model::AM,
    F107::FF;
    satellite_mass::Number,
    satellite_mean_area::Number,
    num_sampling_points_per_orbit::Int,
    abstol::Number,
    C_d::Number,
    C_r::Number,
    distance_unit::Symbol,
    reltol::Number,
    return_solution::Bool,
    solver,
    terminate_altitude::Number,
    tf::Number,
    time_unit::Symbol,
) where {AM, FF}
    M = true_to_mean_anomaly(orb.e, orb.f)
    ā, ē, ī, Ω̄, ω̄, M̄ = _osculating_to_mean_elements(orb.a, orb.e, orb.i, orb.Ω, orb.ω, M)
    u = Vector(_classical_to_equinoctial(ā, ē, ī, Ω̄, ω̄, M̄))

    # Force a homogeneous time span even if the user passes an integer `tf`.
    tspan = (0.0, Float64(tf))

    # Pre-allocate the J2 osculating propagator used for the mean-to-osculating
    # conversion inside the drag quadrature, avoiding one propagator allocation per
    # sampling point.
    j2osc_prop = Propagators.init(
        Val(:J2osc),
        KeplerianElements(0.0, ā, ē, ī, Ω̄, ω̄, mean_to_true_anomaly(ē, M̄))
    )

    # Pre-allocate the Legendre buffers used by the gravity model, avoiding two matrix
    # allocations per acceleration evaluation. The size supports the maximum degree 7 used
    # by `_perturbational_gravity_acceleration`.
    gravity_P  = Matrix{Float64}(undef, 8, 8)
    gravity_dP = Matrix{Float64}(undef, 8, 8)

    # NOTE: `params` carries per-call mutable workspaces: the propagator `j2osc_prop` and
    # the gravity model buffers `gravity_P` and `gravity_dP`. Hence, the assembled ODE
    # problem must not be shared across concurrent solves. Every call to `decay_analysis`
    # builds fresh workspaces, keeping the public API thread-safe.
    params = (
        satellite_mean_area           = satellite_mean_area,
        satellite_mass                = satellite_mass,
        num_sampling_points_per_orbit = num_sampling_points_per_orbit,
        atmospheric_model             = atmospheric_model,
        C_d                           = C_d,
        C_r                           = C_r,
        F107                          = F107,
        gm                            = gm,
        gravity_P                     = gravity_P,
        gravity_dP                    = gravity_dP,
        jd₀_utc                       = orb.t,
        j2osc_prop                    = j2osc_prop,
        terminate_altitude            = terminate_altitude,
    )

    cbset = CallbackSet(ContinuousCallback(_cb_altitude_condition, _cb_altitude_affect!))

    prob = ODEProblem(_dynamics!, u, tspan, params)
    sol = solve(
        prob,
        solver;
        dt       = 24.0 * 60 * 60,
        reltol   = reltol,
        abstol   = abstol,
        callback = cbset,
        maxiters = 1e8
    )

    # == Assemble the Output ===============================================================

    num_points = length(sol.t)

    date             = Vector{DateTime}(undef, num_points)
    time             = Vector{Float64}(undef, num_points)
    f107             = Vector{Float64}(undef, num_points)
    mean_elements    = Vector{KeplerianElements{Float64, Float64}}(undef, num_points)
    apogee_altitude  = Vector{Float64}(undef, num_points)
    perigee_altitude = Vector{Float64}(undef, num_points)

    @inbounds for k in 1:num_points
        aₖ, eₖ, iₖ, Ωₖ, ωₖ, Mₖ = _equinoctial_to_classical(sol.u[k])

        jdₖ = orb.t + sol.t[k] / 86400

        date[k] = julian2datetime(jdₖ)
        time[k] = sol.t[k]

        # Record the space indices used by the dynamics at each instant.
        f107[k] = F107(jdₖ)

        mean_elements[k] = KeplerianElements(
            jdₖ, aₖ, eₖ, iₖ, Ωₖ, ωₖ, mean_to_true_anomaly(eₖ, Mₖ)
        )

        apogee_altitude[k]  = aₖ * (1 + eₖ) - EARTH_EQUATORIAL_RADIUS
        perigee_altitude[k] = aₖ * (1 - eₖ) - EARTH_EQUATORIAL_RADIUS
    end

    # Convert the time and altitude columns to the selected units.
    if time_unit == :m
        time ./= 60
    elseif time_unit == :h
        time ./= 3600
    elseif time_unit == :d
        time ./= 86400
    elseif time_unit != :s
        # Julian year, consistent with the default `tf` of 30 years. If the symbol is not
        # known, we must use the default unit (years).
        time ./= 365.25 * 86400
        time_unit = :y
    end

    if distance_unit != :m
        # If the symbol is not known, we must use the default unit (kilometers).
        apogee_altitude  ./= 1000
        perigee_altitude ./= 1000
        distance_unit = :km
    end

    df = DataFrame(;
        date             = date,
        time             = time,
        f107             = f107,
        mean_elements    = mean_elements,
        apogee_altitude  = apogee_altitude,
        perigee_altitude = perigee_altitude,
    )

    # The style `:note` makes the metadata propagate through DataFrame transformations.
    metadata!(
        df,
        "Description",
        "Mean orbital element evolution during the orbital decay.";
        style = :note
    )

    metadata!(df, "Satellite Mass",      satellite_mass;      style = :note)
    metadata!(df, "Satellite Mean Area", satellite_mean_area; style = :note)
    metadata!(df, "Terminate Altitude",  terminate_altitude;  style = :note)

    colmetadata!(df, :date,             "Unit", :UTC)
    colmetadata!(df, :time,             "Unit", time_unit)
    colmetadata!(df, :f107,             "Unit", :sfu)
    colmetadata!(df, :mean_elements,    "Unit", :SI)
    colmetadata!(df, :apogee_altitude,  "Unit", distance_unit)
    colmetadata!(df, :perigee_altitude, "Unit", distance_unit)

    return_solution && return df, sol

    return df
end

"""
    struct _Nrlmsise00AtmosphericModel

Default atmospheric model of the decay analysis, wrapping the NRLMSISE-00 model provided
by **AtmosphericModels.jl** with a constant geomagnetic index Ap = 9, as in STELA.

# Fields

- `P::Matrix{Float64}`: Pre-allocated buffer for the Legendre matrix used by the model,
    avoiding one matrix allocation per density evaluation. Since it is mutated at every
    evaluation, an instance must not be shared across concurrent computations.
"""
struct _Nrlmsise00AtmosphericModel
    P::Matrix{Float64}

    _Nrlmsise00AtmosphericModel() = new(Matrix{Float64}(undef, 8, 4))
end

"""
    (m::_Nrlmsise00AtmosphericModel)(jd_utc::Number, lat::Number, lon::Number, h::Number, F107::Number) -> Float64

Compute the atmospheric density [kg/m³] using the NRLMSISE-00 model at the Julian date
`jd_utc` [UTC], geodetic latitude `lat` [rad], longitude `lon` [rad], and altitude `h` [m],
considering the 10.7 cm solar flux index `F107` [sfu] as both the daily value and the
81-day centered average.
"""
function (m::_Nrlmsise00AtmosphericModel)(
    jd_utc::Number, lat::Number, lon::Number, h::Number, F107::Number
)
    # We use the default Ap value of 9 as in STELA. Notice that the influence of Ap on the
    # atmospheric density is much smaller than that of F10.7, and the geomagnetic index is
    # not known in advance for a decay analysis.
    Ap    = 9
    atmos = AtmosphericModels.nrlmsise00(jd_utc, h, lat, lon, F107, F107, Ap; P = m.P)
    ρ     = atmos.total_density

    return ρ
end

function _cb_altitude_condition(u, t, integrator)
    a_mean, e_mean, _, _, _, _ = _equinoctial_to_classical(u)

    # Terminate on the perigee altitude of the mean orbit. This condition is a smooth,
    # monotonically decreasing function of the mean elements, whereas the instantaneous
    # altitude oscillates between the perigee and apogee altitudes within a single
    # integrator step, which can make the callback root finder miss the crossing.
    perigee_altitude = a_mean * (1 - e_mean) - EARTH_EQUATORIAL_RADIUS

    return perigee_altitude - integrator.p.terminate_altitude
end

function _cb_altitude_affect!(integrator)
    terminate!(integrator)
    return nothing
end
