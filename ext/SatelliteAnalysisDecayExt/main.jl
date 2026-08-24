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

# Maximum degree kept in the truncated default gravity model file. The decay analysis
# evaluates the gravity field up to degree 7, so the truncation is transparent while
# making the model parsing about two orders of magnitude faster than with the full
# EGM2008 file.
const _DEFAULT_GRAVITY_MODEL_MAX_DEGREE = 36

"""
    _default_gravity_model_file() -> String

Return the path of the truncated EGM2008 file used by the default gravity model, building
it in the package scratch space at the first call ever. The truncation keeps the
coefficients up to the degree `_DEFAULT_GRAVITY_MODEL_MAX_DEGREE` and is transparent for
the decay analysis, which evaluates the gravity field up to degree 7.
"""
function _default_gravity_model_file()
    dir  = SatelliteAnalysis.Scratch.get_scratch!(SatelliteAnalysis, "decay_analysis")
    path = joinpath(dir, "EGM2008_max_degree_$(_DEFAULT_GRAVITY_MODEL_MAX_DEGREE).gfc")

    isfile(path) && return path

    source = fetch_icgem_file(:EGM2008)
    tmp    = path * ".tmp"

    n    = _DEFAULT_GRAVITY_MODEL_MAX_DEGREE
    seen = falses(n + 1, n + 1)

    open(tmp, "w") do out
        for line in eachline(source)
            if startswith(line, "gfc")
                tokens = split(line; limit = 4)
                degree = parse(Int, tokens[2])

                if degree <= n
                    order = parse(Int, tokens[3])
                    seen[degree + 1, order + 1] = true
                    println(out, line)
                end
            elseif startswith(line, "max_degree")
                println(out, "max_degree ", n)
            else
                println(out, line)
            end
        end

        # The source file omits the null coefficients (e.g. degree 1 in EGM2008), which
        # the parser would leave uninitialized. Hence, we explicitly write the missing
        # entries as zeros.
        for degree in 0:n, order in 0:degree
            seen[degree + 1, order + 1] && continue
            println(out, "gfc ", degree, " ", order, " 0.0 0.0 0.0 0.0")
        end
    end

    # Move the finished file atomically so that an interrupted build does not leave a
    # partial file behind.
    mv(tmp, path; force = true)

    return path
end

function SatelliteAnalysis.decay_analysis(
    orb::KeplerianElements;
    # Required keywords.
    satellite_mass::Number,
    satellite_mean_area::Number,
    # Optional keywords. The keywords `atmospheric_model` and `F107` accept any callable
    # object, hence they are not annotated with `::Function`.
    atmospheric_model = nothing,
    atmospheric_model_name::Union{Nothing, String} = nothing,
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
    verbose::Bool = false,
)
    gm = if isnothing(gravity_model)
        if isnothing(_DEFAULT_GRAVITY_MODEL[])
            _DEFAULT_GRAVITY_MODEL[] =
                GravityModels.load(IcgemFile, _default_gravity_model_file())
        end

        _DEFAULT_GRAVITY_MODEL[]
    else
        gravity_model
    end

    # The default atmospheric model carries a mutable buffer, so a fresh instance is built
    # per call to keep the public API thread-safe.
    atmospheric_model′ = isnothing(atmospheric_model) ?
        Nrlmsise00AtmosphericModel() :
        atmospheric_model

    # Descriptions of the atmospheric model and the F10.7 source, recorded as metadata in
    # the output so that, for example, `plot_decay_analysis` can show the assumptions. The
    # name provided by the user has precedence over the derived one.
    atmospheric_model_name′ = if !isnothing(atmospheric_model_name)
        atmospheric_model_name
    elseif isnothing(atmospheric_model)
        "NRLMSISE-00"
    elseif atmospheric_model isa Function
        "Custom (" * String(nameof(atmospheric_model)) * ")"
    else
        "Custom (" * String(nameof(typeof(atmospheric_model))) * ")"
    end

    f107_source = isnothing(F107) ? "Predicted" :
        F107 isa Number ? "Constant ($(F107) sfu)" : "User function"

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
        atmospheric_model_name        = atmospheric_model_name′,
        f107_source                   = f107_source,
        C_d                           = C_d,
        C_r                           = C_r,
        distance_unit                 = distance_unit,
        reltol                        = reltol,
        return_solution               = return_solution,
        solver                        = solver,
        terminate_altitude            = terminate_altitude,
        tf                            = tf,
        time_unit                     = time_unit,
        verbose                       = verbose
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
    atmospheric_model_name::String,
    f107_source::String,
    C_d::Number,
    C_r::Number,
    distance_unit::Symbol,
    reltol::Number,
    return_solution::Bool,
    solver,
    terminate_altitude::Number,
    tf::Number,
    time_unit::Symbol,
    verbose::Bool,
) where {AM, FF}
    M = true_to_mean_anomaly(orb.e, orb.f)
    ā, ē, ī, Ω̄, ω̄, M̄ = _osculating_to_mean_elements(orb.a, orb.e, orb.i, orb.Ω, orb.ω, M)

    # The state is a static vector with an out-of-place right-hand side, removing the
    # per-step state allocations of the numerical integration.
    u₀ = _classical_to_equinoctial(ā, ē, ī, Ω̄, ω̄, M̄)

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

    # Hoist the gravity model constants used by the hot loop out of the right-hand side.
    μ  = GravityModels.gravity_constant(gm)
    Re = GravityModels.radius(gm)
    J₂ = -first(GravityModels.coefficients(gm, 2, 0)) * √5

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
        μ                             = μ,
        Re                            = Re,
        J₂                            = J₂,
        gravity_P                     = gravity_P,
        gravity_dP                    = gravity_dP,
        jd₀_utc                       = orb.t,
        j2osc_prop                    = j2osc_prop,
        terminate_altitude            = terminate_altitude,
    )

    # Progress interface shown during the integration when `verbose` is enabled.
    progress = verbose ?
        DecayProgress(
            stderr,
            tspan[2],
            ā * (1 - ē) - EARTH_EQUATORIAL_RADIUS,
            Float64(terminate_altitude)
        ) :
        nothing

    # The progress callback is always installed with a stable type: when `verbose` is
    # disabled, its condition is constantly `false`. Hence, the solver specialization is
    # shared by the verbose and silent paths, avoiding a second compilation on the first
    # verbose call.
    cbset = CallbackSet(
        ContinuousCallback(_cb_altitude_condition, _cb_altitude_affect!),
        _decay_progress_callback(progress)
    )

    isnothing(progress) ||
        _start_decay_progress!(progress, ā * (1 + ē) - EARTH_EQUATORIAL_RADIUS)

    prob = ODEProblem(_dynamics, u₀, tspan, params)

    sol = try
        solve(
            prob,
            solver;
            dt       = 24.0 * 60 * 60,
            reltol   = reltol,
            abstol   = abstol,
            callback = cbset,
            maxiters = 1e8
        )
    finally
        # Restore the terminal cursor even if the solver throws.
        isnothing(progress) || _cleanup_decay_progress!(progress)
    end

    if !isnothing(progress)
        aₑ, eₑ, _, _, _, _ = _equinoctial_to_classical(sol.u[end])

        perigee_end = aₑ * (1 - eₑ) - EARTH_EQUATORIAL_RADIUS
        apogee_end  = aₑ * (1 + eₑ) - EARTH_EQUATORIAL_RADIUS
        reentered   = perigee_end <= terminate_altitude * (1 + 1e-6)

        _finish_decay_progress!(progress, sol.t[end], perigee_end, apogee_end, reentered)
    end

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

    metadata!(df, "Atmospheric Model",   atmospheric_model_name; style = :note)
    metadata!(df, "Drag Coefficient",    C_d;                    style = :note)
    metadata!(df, "F10.7 Source",        f107_source;            style = :note)
    metadata!(df, "Satellite Mass",      satellite_mass;         style = :note)
    metadata!(df, "Satellite Mean Area", satellite_mean_area;    style = :note)
    metadata!(df, "SRP Coefficient",     C_r;                    style = :note)
    metadata!(df, "Terminate Altitude",  terminate_altitude;     style = :note)

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
    struct Nrlmsise00AtmosphericModel

Default atmospheric model of the decay analysis, wrapping the NRLMSISE-00 model provided
by **AtmosphericModels.jl** with a constant geomagnetic index Ap = 9, as in STELA.

# Fields

- `P::Matrix{Float64}`: Pre-allocated buffer for the Legendre matrix used by the model,
    avoiding one matrix allocation per density evaluation. Since it is mutated at every
    evaluation, an instance must not be shared across concurrent computations.
"""
struct Nrlmsise00AtmosphericModel
    P::Matrix{Float64}

    Nrlmsise00AtmosphericModel() = new(Matrix{Float64}(undef, 8, 4))
end

"""
    (m::Nrlmsise00AtmosphericModel)(jd_utc::Number, lat::Number, lon::Number, h::Number, F107::Number) -> Float64

Compute the atmospheric density [kg/m³] using the NRLMSISE-00 model at the Julian date
`jd_utc` [UTC], geodetic latitude `lat` [rad], longitude `lon` [rad], and altitude `h` [m],
considering the 10.7 cm solar flux index `F107` [sfu] as both the daily value and the
81-day centered average.
"""
function (m::Nrlmsise00AtmosphericModel)(
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
