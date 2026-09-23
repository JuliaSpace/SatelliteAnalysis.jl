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
    input_type::Symbol = :mean,
    num_sampling_points_per_orbit::Union{Nothing, Int} = nothing,
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
    input_type in (:mean, :osculating) ||
        throw(ArgumentError("The keyword `input_type` must be `:mean` or `:osculating`."))

    # Validate the physical inputs. Otherwise, the analysis can silently produce `NaN`s or
    # infinite values.
    satellite_mass > 0 || throw(ArgumentError("The satellite mass must be greater than 0."))

    satellite_mean_area >= 0 ||
        throw(ArgumentError("The satellite mean area must not be negative."))

    (isnothing(num_sampling_points_per_orbit) || (num_sampling_points_per_orbit >= 1)) ||
        throw(
            ArgumentError("The number of sampling points per orbit must be greater than 0.")
        )

    C_d >= 0 || throw(ArgumentError("The drag coefficient must not be negative."))

    C_r >= 0 || throw(
        ArgumentError("The solar radiation pressure coefficient must not be negative.")
    )

    abstol > 0 || throw(ArgumentError("The absolute tolerance must be greater than 0."))
    reltol > 0 || throw(ArgumentError("The relative tolerance must be greater than 0."))

    terminate_altitude >= 0 ||
        throw(ArgumentError("The terminate altitude must not be negative."))

    tf > 0 || throw(ArgumentError("The maximum propagation time must be greater than 0."))

    # Validate the units before performing the analysis.
    _decay_analysis__time_unit_factor(time_unit)
    SatelliteAnalysis._distance_unit_factor(distance_unit)

    # The propagation uses mean elements with respect to the averaged dynamics. Hence, if
    # the input elements are osculating, they must be converted to mean elements first.
    orb′ = input_type == :osculating ? _osculating_to_mean_elements(orb) : orb

    # If the user did not select the number of sampling points per orbit, we must obtain it
    # from the mean eccentricity.
    num_sampling_points_per_orbit′ =
        isnothing(num_sampling_points_per_orbit) ?
        _decay_analysis__default_num_sampling_points(orb′.eccentricity) :
        num_sampling_points_per_orbit

    gm = if isnothing(gravity_model)
        if isnothing(_DEFAULT_GRAVITY_MODEL[])
            _DEFAULT_GRAVITY_MODEL[] = GravityModels.load(IcgemFile, fetch_icgem_file(:EGM96))
        end

        _DEFAULT_GRAVITY_MODEL[]
    else
        gravity_model
    end

    # The default atmospheric model carries a mutable buffer, so a fresh instance is built
    # per call to keep the public API thread-safe.
    atmospheric_model′ =
        isnothing(atmospheric_model) ? Nrlmsise00AtmosphericModel() : atmospheric_model

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

    space_indices′ = if isnothing(space_indices)
        _decay_analysis__default_space_indices
    elseif space_indices isa NamedTuple
        _ -> space_indices
    else
        space_indices
    end

    # Description of the space indices source and flag indicating whether it is one of the
    # default sources, which require the initialization of the space indices.
    space_indices_source, is_default_space_indices =
        if space_indices′ === _decay_analysis__default_space_indices
            "Default (Obs. + Pred.)", true
        elseif space_indices′ === _decay_analysis__default_space_indices_kp
            "Default (Adj. + Pred.)", true
        elseif space_indices isa NamedTuple
            "Constant $(space_indices)", false
        else
            "User function", false
        end

    # The default space index functions require the remote data sets provided by
    # SpaceIndices.jl. Notice that if the user calls `SpaceIndices.destroy()` after this
    # initialization, `space_index` raises a clear error asking to initialize the space
    # indices again.
    if is_default_space_indices && !_SPACE_INDICES_INITIALIZED[]
        SpaceIndices.init(SpaceIndices.Celestrak)
        SpaceIndices.init(SpaceIndices.SatelliteToolboxSpaceIndexSets)
        _SPACE_INDICES_INITIALIZED[] = true
    end

    # The keyword `gravity_model` is abstractly typed and the callables in
    # `atmospheric_model` and `space_indices′` have call-site dependent types. This
    # function barrier ensures they are concretely typed inside the numerical integration,
    # avoiding dynamic dispatch at every right-hand-side evaluation. The callables are
    # passed as positional arguments since keyword arguments cannot bind the type
    # parameters that force the specialization.
    return _decay_analysis(
        orb′,
        gm,
        atmospheric_model′,
        space_indices′;
        satellite_mass                = satellite_mass,
        satellite_mean_area           = satellite_mean_area,
        num_sampling_points_per_orbit = num_sampling_points_per_orbit′,
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
        verbose                       = verbose,
    )
end

# Setups selected by the macros `@decay_analysis__jacchia77`,
# `@decay_analysis__jacchia77_stela`, and `@decay_analysis__jr1971`: the atmospheric model
# wrapper and the default space indices source of each model. The `Nothing` argument
# overrides the `::Any` fallback in the main package with a more specific method instead
# of overwriting it, which is forbidden during precompilation.
function SatelliteAnalysis._decay_analysis__jacchia77_setup(::Nothing)
    return (
        atmospheric_model      = Jacchia77AtmosphericModel(),
        atmospheric_model_name = "Jacchia 1977",
        space_indices          = _decay_analysis__default_space_indices_kp,
    )
end

function SatelliteAnalysis._decay_analysis__jacchia77_stela_setup(::Nothing)
    return (
        atmospheric_model      = Jacchia77AtmosphericModel(Val(:stela)),
        atmospheric_model_name = "Jacchia 1977 (STELA)",
        space_indices          = _decay_analysis__default_space_indices_kp,
    )
end

function SatelliteAnalysis._decay_analysis__jr1971_setup(::Nothing)
    return (
        atmospheric_model      = Jr1971AtmosphericModel(),
        atmospheric_model_name = "Jacchia-Roberts 1971",
        space_indices          = _decay_analysis__default_space_indices_kp,
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
    # The input elements are treated as mean elements with respect to the averaged
    # dynamics, following the same convention of semi-analytical tools such as STELA. The
    # analysis is computed in `Float64`, and the state of the numerical integration stores
    # the mean anomaly. Hence, we convert the elements once so that the initial state, the
    # pre-allocated propagator, and the output share the same types.
    orb = convert(KeplerianElements{MeanAnomaly, Float64, Float64}, orb)

    # The state is a static vector with an out-of-place right-hand side, removing the
    # per-step state allocations of the numerical integration.
    u₀ = _keplerian_to_state(orb)

    # The analysis is meaningless if the initial mean perigee altitude is already at or
    # below the terminate altitude, since the termination condition could never be
    # triggered.
    perigee₀, apogee₀ = _state_apsis_altitudes(u₀)

    perigee₀ > terminate_altitude || throw(
        ArgumentError(
            "The initial mean perigee altitude ($(perigee₀) m) must be greater than the terminate altitude ($(terminate_altitude) m).",
        ),
    )

    # Force a homogeneous time span even if the user passes an integer `tf`.
    tspan = (0.0, Float64(tf))

    # Pre-allocate the J2 osculating propagator used for the mean-to-osculating
    # conversion inside the drag quadrature, avoiding one propagator allocation per
    # sampling point.
    j2osc_prop = Propagators.init(Val(:J2osc), orb)

    # Pre-allocate the gravity model workspace, avoiding the buffer allocations at every
    # acceleration evaluation. It supports the maximum degree 7 and order 0 used by
    # `_perturbational_gravity_acceleration`. The element type must be the one obtained by
    # promoting the model coefficients with the state type (`Float64`).
    gravity_workspace = GravityModels.Workspace(
        gm; max_degree = 7, max_order = 0, T = _gravity_workspace_type(gm)
    )

    # Hoist the gravity model constants used by the hot loop out of the right-hand side.
    μ  = GravityModels.gravity_constant(gm)
    Re = GravityModels.radius(gm)
    J₂ = -first(GravityModels.coefficients(gm, 2, 0)) * √5

    # NOTE: `params` carries per-call mutable workspaces: the propagator `j2osc_prop` and
    # the gravity model workspace `gravity_workspace`. Hence, the assembled ODE problem
    # must not be shared across concurrent solves. Every call to `decay_analysis` builds
    # fresh workspaces, keeping the public API thread-safe.
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
        gravity_workspace             = gravity_workspace,
        jd₀_utc                       = orb.epoch,
        j2osc_prop                    = j2osc_prop,
        terminate_altitude            = terminate_altitude,
    )

    # Progress interface shown during the integration when `verbose` is enabled.
    progress =
        verbose ? DecayProgress(stderr, tspan[2], perigee₀, Float64(terminate_altitude)) :
        nothing

    # The progress callback is always installed with a stable type: when `verbose` is
    # disabled, its condition is constantly `false`. Hence, the solver specialization is
    # shared by the verbose and silent paths, avoiding a second compilation on the first
    # verbose call.
    cbset = CallbackSet(
        # Only a downward crossing of the terminate altitude ends the integration.
        ContinuousCallback(_cb_altitude_condition, nothing, _cb_altitude_affect!),
        _decay_progress_callback(progress),
    )

    isnothing(progress) || _start_decay_progress!(progress, apogee₀)

    prob = ODEProblem(_dynamics, u₀, tspan, params)

    sol = try
        solve(
            prob,
            solver;
            dt       = 24.0 * 60 * 60,
            reltol   = reltol,
            abstol   = abstol,
            callback = cbset,
            maxiters = 1e8,
        )
    finally
        # Restore the terminal cursor even if the solver throws.
        isnothing(progress) || _cleanup_decay_progress!(progress)
    end

    if !isnothing(progress)
        perigee_end, apogee_end = _state_apsis_altitudes(sol.u[end])

        # The integration is only terminated by the callback that detects when the perigee
        # altitude reaches the terminate altitude.
        reentered = sol.retcode == ReturnCode.Terminated

        # Any other return code, except `Success`, means that the solver stopped before
        # reaching the maximum propagation time.
        failure =
            (reentered || sol.retcode == ReturnCode.Success) ? nothing : string(sol.retcode)

        _finish_decay_progress!(
            progress, sol.t[end], perigee_end, apogee_end, reentered, failure
        )
    end

    # == Assemble the Output ===============================================================

    num_points = length(sol.t)

    date             = Vector{DateTime}(undef, num_points)
    time             = Vector{Float64}(undef, num_points)
    mean_elements    = Vector{KeplerianElements{MeanAnomaly, Float64, Float64}}(undef, num_points)
    apogee_altitude  = Vector{Float64}(undef, num_points)
    perigee_altitude = Vector{Float64}(undef, num_points)

    @inbounds for k in 1:num_points
        jdₖ = orb.epoch + sol.t[k] / 86400

        date[k] = julian2datetime(jdₖ)
        time[k] = sol.t[k]

        # The state already contains the mean anomaly. Hence, the conversion does not solve
        # Kepler's equation.
        mean_elements[k] = _state_to_keplerian(sol.u[k], jdₖ)

        perigee_altitude[k], apogee_altitude[k] = _state_apsis_altitudes(sol.u[k])
    end

    # Record the space indices used by the dynamics at each instant. A comprehension is
    # used so that the column eltype is the concrete named tuple type of the source.
    space_indices_column = [space_indices(orb.epoch + tₖ / 86400) for tₖ in sol.t]

    # Convert the time and altitude columns to the selected units. Notice that the year is
    # the Julian year, consistent with the default `tf` of 30 years.
    time_factor     = _decay_analysis__time_unit_factor(time_unit)
    distance_factor = SatelliteAnalysis._distance_unit_factor(distance_unit)

    time             .*= time_factor
    apogee_altitude  .*= distance_factor
    perigee_altitude .*= distance_factor

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
        style = :note,
    )

    metadata!(df, "Atmospheric Model", atmospheric_model_name; style = :note)
    metadata!(df, "Drag Coefficient", C_d; style = :note)
    metadata!(df, "Satellite Mass", satellite_mass; style = :note)
    metadata!(df, "Satellite Mean Area", satellite_mean_area; style = :note)
    metadata!(df, "Space Indices Source", space_indices_source; style = :note)
    metadata!(df, "SRP Coefficient", C_r; style = :note)
    metadata!(df, "Terminate Altitude", terminate_altitude; style = :note)

    # Notice that the column `space_indices` has no `Unit` metadata since its fields have
    # heterogeneous units.
    colmetadata!(df, :date, "Unit", :UTC; style = :note)
    colmetadata!(df, :time, "Unit", time_unit; style = :note)
    colmetadata!(df, :mean_elements, "Unit", :SI; style = :note)
    colmetadata!(df, :apogee_altitude, "Unit", distance_unit; style = :note)
    colmetadata!(df, :perigee_altitude, "Unit", distance_unit; style = :note)

    # The raw solution is attached to the `DataFrame`, keeping the return type stable.
    # Notice that this metadata has the style `:default`. Hence, it is not propagated by
    # DataFrame transformations because it would not describe the new table.
    return_solution && metadata!(df, "Solution", sol; style = :default)

    return df
end

"""
    _decay_analysis__default_space_indices(jd_utc::Number) -> NamedTuple

Default space indices computed at the Julian Day [UTC] `jd_utc` used by the decay analysis
when the user does not provide a named tuple or a callable. It returns a named tuple with
the fields:

- `f107`: Observed daily 10.7 cm solar flux of the previous day [sfu], falling back to
    the predicted F10.7 outside the observed timespan. The NRLMSISE-00 documentation
    requires the observed flux (measured at the actual Earth-Sun distance) instead of the
    flux adjusted to 1 AU, and prescribes the value of the previous day for the daily
    index.
- `f107_avg`: Observed centered 81-day average of the 10.7 cm solar flux [sfu], falling
    back to the predicted F10.7 outside the observed timespan.
- `ap`: Observed daily geomagnetic index [-], falling back to 9, as in STELA, outside the
    observed timespan. Notice that the influence of Ap on the atmospheric density is much
    smaller than that of F10.7, and the geomagnetic index is not known in advance for a
    decay analysis.
"""
function _decay_analysis__default_space_indices(jd_utc::Number)
    jd₀, jd₁ = SpaceIndices.timespan(Val(:Ap_daily))
    ap = jd₀ <= jd_utc <= jd₁ ? space_index(Val(:Ap_daily), jd_utc) : 9.0

    # The NRLMSISE-00 documentation prescribes the daily F10.7 of the previous day.
    return (
        f107     = _default_f107(jd_utc - 1, Val(:F10obs), Val(:F10obs_predicted)),
        f107_avg = _default_f107(jd_utc, Val(:F10obs_avg_center81), Val(:F10obs_predicted)),
        ap       = ap,
    )
end

"""
    _decay_analysis__default_space_indices_kp(jd_utc::Number) -> NamedTuple

Default space indices computed at the Julian Day [UTC] `jd_utc` used by the decay analysis
with the atmospheric models that require the geomagnetic index Kp (Jacchia 1977 and
Jacchia-Roberts 1971), selected by the macros `@decay_analysis__jacchia77` and
`@decay_analysis__jr1971`. It returns a named tuple with the fields:

- `f107`: Daily 10.7 cm solar flux adjusted to 1 AU [sfu], falling back to the predicted
    F10.7 outside the adjusted flux timespan. The adjusted flux is used because the
    Jacchia models were derived using the flux normalized to 1 AU, unlike NRLMSISE-00,
    which uses the observed flux at the actual Earth-Sun distance.
- `f107_avg`: Centered 81-day average of the 10.7 cm solar flux adjusted to 1 AU [sfu],
    falling back to the predicted F10.7 outside the adjusted flux timespan.
- `kp`: Observed daily geomagnetic index Kp [-], falling back to 7 / 3 (equivalent to
    Ap = 9, as in STELA) outside the observed timespan. Notice that the influence of the
    geomagnetic index on the atmospheric density is much smaller than that of F10.7, and
    it is not known in advance for a decay analysis.
"""
function _decay_analysis__default_space_indices_kp(jd_utc::Number)
    jd₀, jd₁ = SpaceIndices.timespan(Val(:Kp_daily))
    kp = jd₀ <= jd_utc <= jd₁ ? space_index(Val(:Kp_daily), jd_utc) : 7 / 3

    return (
        f107     = _default_f107(jd_utc, Val(:F10adj), Val(:F10adj_predicted)),
        f107_avg = _default_f107(jd_utc, Val(:F10adj_avg_center81), Val(:F10adj_predicted)),
        kp       = kp,
    )
end

# F10.7 [sfu] from the space index `index` at the Julian Day [UTC] `jd_utc`, falling back
# to the space index `predicted` outside the `index` timespan. The default sources use
# the observed indices (`F10obs`, `F10obs_avg_center81`, and `F10obs_predicted`) for
# NRLMSISE-00, whose documentation requires the flux at the actual Earth-Sun distance,
# and the adjusted indices (`F10adj`, `F10adj_avg_center81`, and `F10adj_predicted`) for
# the Jacchia models, which were derived using the flux normalized to 1 AU. Both use a
# centered 81-day average, the convention shared by the supported atmospheric models (the
# Jacchia 1977 report uses a centered Gaussian-weighted mean, which the centered 81-day
# average approximates).
function _default_f107(jd_utc::Number, index::Val, predicted::Val)
    jd₀, jd₁ = SpaceIndices.timespan(index)

    jd₀ <= jd_utc <= jd₁ && return space_index(index, jd_utc)

    return space_index(predicted, jd_utc)
end

# Default number of sampling points per orbit used to average the perturbations given the
# mean eccentricity `e`. The perturbations are concentrated near the perigee in eccentric
# orbits, requiring more points. Those values keep the quadrature error of the lifetime
# below 0.03% in the following reference cases compared to the result with 257 points: a
# 300 km circular orbit (17 points), a 250 km x 2000 km orbit (e = 0.12, 33 points), and a
# 250 km x 20000 km orbit (e = 0.60, 65 points). With 17 points, the errors in the last two
# cases are 0.02% and 7.6%, respectively.
function _decay_analysis__default_num_sampling_points(e::Number)
    e < 0.05 && return 17
    e < 0.30 && return 33
    return 65
end

# Factor to convert the time in seconds to the unit `time_unit`, which can be `:s`, `:min`,
# `:h`, `:d`, or `:y` (Julian year). It throws an `ArgumentError` if the unit is not valid.
function _decay_analysis__time_unit_factor(time_unit::Symbol)
    return SatelliteAnalysis._time_unit_factor(time_unit, (:s, :min, :h, :d, :y))
end

# Element type of the gravity model workspace: the promotion between the type of the
# coefficients in the gravity model `gm` and the state type (`Float64`).
function _gravity_workspace_type(::AbstractGravityModel{Tm}) where {Tm <: Number}
    return promote_type(float(Tm), Float64)
end

function _cb_altitude_condition(u, t, integrator)
    # Terminate on the perigee altitude of the mean orbit. This condition is a smooth,
    # monotonically decreasing function of the mean elements, whereas the instantaneous
    # altitude oscillates between the perigee and apogee altitudes within a single
    # integrator step, which can make the callback root finder miss the crossing.
    perigee_altitude, ~ = _state_apsis_altitudes(u)

    return perigee_altitude - integrator.p.terminate_altitude
end

function _cb_altitude_affect!(integrator)
    terminate!(integrator)
    return nothing
end
