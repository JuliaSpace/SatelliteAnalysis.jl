## Description #############################################################################
#
# Main function to compute the decay analysis of a satellite in orbit.
#
############################################################################################

# Flag indicating whether the space index sets required by the default `space_indices`
# keyword were already initialized in this session, avoiding re-fetching and re-parsing the
# remote files at every call. A benign race that initializes the sets twice is acceptable.
const _SPACE_INDICES_INITIALIZED = Ref(false)

# Cache of the default gravity model (EGM96), avoiding re-fetching and re-parsing the ICGEM
# file at every call. A benign race that loads the model twice is acceptable.
const _DEFAULT_GRAVITY_MODEL = Ref{Union{Nothing, AbstractGravityModel}}(nothing)

function SatelliteAnalysis.decay_analysis(
    orb::KeplerianElements;
    # Required keywords.
    satellite_mass::Number,
    satellite_mean_area::Number,
    # Optional keywords. The keywords `atmospheric_model` and `space_indices` accept any
    # callable object, hence they are not annotated with `::Function`.
    atmospheric_model = nothing,
    atmospheric_model_name::Union{Nothing, String} = nothing,
    gravity_model::Union{AbstractGravityModel, Nothing} = nothing,
    num_sampling_points_per_orbit::Int = 17,
    abstol::Number = 1e-6,
    C_d::Number = 2.2,
    C_r::Number = 1.25,
    distance_unit::Symbol = :km,
    space_indices = nothing,
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
                GravityModels.load(IcgemFile, fetch_icgem_file(:EGM96))
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

    # Descriptions of the atmospheric model and the space indices source, recorded as
    # metadata in the output so that, for example, `plot_decay_analysis` can show the
    # assumptions. The name provided by the user has precedence over the derived one.
    atmospheric_model_name′ = if !isnothing(atmospheric_model_name)
        atmospheric_model_name
    elseif isnothing(atmospheric_model)
        "NRLMSISE-00"
    elseif atmospheric_model isa Function
        "Custom (" * String(nameof(atmospheric_model)) * ")"
    else
        "Custom (" * String(nameof(typeof(atmospheric_model))) * ")"
    end

    space_indices_source = isnothing(space_indices) ? "Default (Obs. + Pred.)" :
        space_indices isa NamedTuple ? "Constant $(space_indices)" : "User function"

    space_indices′ = if isnothing(space_indices)
        # Notice that if the user calls `SpaceIndices.destroy()` after this initialization,
        # `space_index` raises a clear error asking to initialize the space indices again.
        if !_SPACE_INDICES_INITIALIZED[]
            SpaceIndices.init(SpaceIndices.Celestrak)
            SpaceIndices.init(SpaceIndices.SatelliteToolboxSpaceIndexSets)
            _SPACE_INDICES_INITIALIZED[] = true
        end

        _decay_analysis__default_space_indices
    elseif space_indices isa NamedTuple
        _ -> space_indices
    else
        space_indices
    end

    # The keyword `gravity_model` is abstractly typed and the callables in
    # `atmospheric_model` and `space_indices′` have call-site dependent types. This
    # function barrier ensures they are concretely typed inside the numerical integration,
    # avoiding dynamic dispatch at every right-hand-side evaluation. The callables are
    # passed as positional arguments since keyword arguments cannot bind the type
    # parameters that force the specialization.
    return _decay_analysis(
        orb,
        gm,
        atmospheric_model′,
        space_indices′;
        satellite_mass                = satellite_mass,
        satellite_mean_area           = satellite_mean_area,
        num_sampling_points_per_orbit = num_sampling_points_per_orbit,
        abstol                        = abstol,
        atmospheric_model_name        = atmospheric_model_name′,
        space_indices_source          = space_indices_source,
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
    space_indices::SF;
    satellite_mass::Number,
    satellite_mean_area::Number,
    num_sampling_points_per_orbit::Int,
    abstol::Number,
    atmospheric_model_name::String,
    space_indices_source::String,
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
) where {AM, SF}
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
        space_indices                 = space_indices,
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
    mean_elements    = Vector{KeplerianElements{Float64, Float64}}(undef, num_points)
    apogee_altitude  = Vector{Float64}(undef, num_points)
    perigee_altitude = Vector{Float64}(undef, num_points)

    @inbounds for k in 1:num_points
        aₖ, eₖ, iₖ, Ωₖ, ωₖ, Mₖ = _equinoctial_to_classical(sol.u[k])

        jdₖ = orb.t + sol.t[k] / 86400

        date[k] = julian2datetime(jdₖ)
        time[k] = sol.t[k]

        mean_elements[k] = KeplerianElements(
            jdₖ, aₖ, eₖ, iₖ, Ωₖ, ωₖ, mean_to_true_anomaly(eₖ, Mₖ)
        )

        apogee_altitude[k]  = aₖ * (1 + eₖ) - EARTH_EQUATORIAL_RADIUS
        perigee_altitude[k] = aₖ * (1 - eₖ) - EARTH_EQUATORIAL_RADIUS
    end

    # Record the space indices used by the dynamics at each instant. A comprehension is
    # used so that the column eltype is the concrete named tuple type of the source.
    space_indices_column = [space_indices(orb.t + tₖ / 86400) for tₖ in sol.t]

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
        space_indices    = space_indices_column,
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

    metadata!(df, "Atmospheric Model",    atmospheric_model_name; style = :note)
    metadata!(df, "Drag Coefficient",     C_d;                    style = :note)
    metadata!(df, "Satellite Mass",       satellite_mass;         style = :note)
    metadata!(df, "Satellite Mean Area",  satellite_mean_area;    style = :note)
    metadata!(df, "Space Indices Source", space_indices_source;   style = :note)
    metadata!(df, "SRP Coefficient",      C_r;                    style = :note)
    metadata!(df, "Terminate Altitude",   terminate_altitude;     style = :note)

    # Notice that the column `space_indices` has no `Unit` metadata since its fields have
    # heterogeneous units.
    colmetadata!(df, :date,             "Unit", :UTC)
    colmetadata!(df, :time,             "Unit", time_unit)
    colmetadata!(df, :mean_elements,    "Unit", :SI)
    colmetadata!(df, :apogee_altitude,  "Unit", distance_unit)
    colmetadata!(df, :perigee_altitude, "Unit", distance_unit)

    return_solution && return df, sol

    return df
end

"""
    _decay_analysis__default_space_indices(jd_utc::Number) -> NamedTuple

Default space indices computed at the Julian Day [UTC] `jd_utc` used by the decay analysis
when the user does not provide a named tuple or a callable. It returns a named tuple with
the fields:

- `f107`: Observed daily 10.7 cm solar flux [sfu], falling back to the predicted F10.7
    outside the observed timespan.
- `f107_avg`: Observed last-81-day average of the 10.7 cm solar flux [sfu], falling back to
    the predicted F10.7 outside the observed timespan.
- `ap`: Observed daily geomagnetic index [-], falling back to 9, as in STELA, outside the
    observed timespan. Notice that the influence of Ap on the atmospheric density is much
    smaller than that of F10.7, and the geomagnetic index is not known in advance for a
    decay analysis.
"""
function _decay_analysis__default_space_indices(jd_utc::Number)
    jd₀, jd₁ = SpaceIndices.timespan(Val(:F10obs))
    f107 = jd₀ <= jd_utc <= jd₁ ?
        space_index(Val(:F10obs), jd_utc) :
        space_index(Val(:F10predicted), jd_utc)

    jd₀, jd₁ = SpaceIndices.timespan(Val(:F10obs_avg_last81))
    f107_avg = jd₀ <= jd_utc <= jd₁ ?
        space_index(Val(:F10obs_avg_last81), jd_utc) :
        space_index(Val(:F10predicted), jd_utc)

    jd₀, jd₁ = SpaceIndices.timespan(Val(:Ap_daily))
    ap = jd₀ <= jd_utc <= jd₁ ? space_index(Val(:Ap_daily), jd_utc) : 9.0

    return (f107 = f107, f107_avg = f107_avg, ap = ap)
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
