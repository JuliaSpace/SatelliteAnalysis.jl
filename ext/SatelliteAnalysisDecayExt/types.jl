## Description #############################################################################
#
# Definition of types and structures for the decay analysis.
#
############################################################################################

"""
    struct Nrlmsise00AtmosphericModel

Default atmospheric model of the decay analysis, wrapping the NRLMSISE-00 model provided
by **AtmosphericModels.jl**. It consumes the space indices `f107` (daily 10.7 cm solar
flux) [sfu], `f107_avg` (81-day average of the 10.7 cm solar flux) [sfu], and `ap` (daily
geomagnetic index) [-] from the named tuple passed to the callable.

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
    (m::Nrlmsise00AtmosphericModel)(
        jd_utc::Number,
        lat::Number,
        lon::Number,
        h::Number,
        space_indices::NamedTuple
    ) -> Float64

Compute the atmospheric density [kg/m³] using the NRLMSISE-00 model at the Julian date
`jd_utc` [UTC], geodetic latitude `lat` [rad], longitude `lon` [rad], and altitude `h` [m],
considering the space indices in the named tuple `space_indices`, which must contain the
fields:

- `f107`: Daily 10.7 cm solar flux [sfu].
- `f107_avg`: 81-day average of the 10.7 cm solar flux [sfu].
- `ap`: Daily geomagnetic index [-].
"""
function (m::Nrlmsise00AtmosphericModel)(
    jd_utc::Number, lat::Number, lon::Number, h::Number, space_indices::NamedTuple
)
    atmos = AtmosphericModels.nrlmsise00(
        jd_utc,
        h,
        lat,
        lon,
        space_indices.f107_avg,
        space_indices.f107,
        space_indices.ap;
        P = m.P
    )

    return atmos.total_density
end

"""
    struct Jacchia77AtmosphericModel

Atmospheric model of the decay analysis wrapping the Jacchia 1977 model provided by
**AtmosphericModels.jl**, selected by the macro
[`@decay_analysis__jacchia77`](@ref SatelliteAnalysis.@decay_analysis__jacchia77). It
consumes the space indices `f107` (daily 10.7 cm solar flux) [sfu], `f107_avg` (81-day
average of the 10.7 cm solar flux) [sfu], and `kp` (daily geomagnetic index Kp) [-] from
the named tuple passed to the callable.
"""
struct Jacchia77AtmosphericModel end

"""
    (m::Jacchia77AtmosphericModel)(
        jd_utc::Number,
        lat::Number,
        lon::Number,
        h::Number,
        space_indices::NamedTuple
    ) -> Float64

Compute the atmospheric density [kg/m³] using the Jacchia 1977 model at the Julian date
`jd_utc` [UTC], geodetic latitude `lat` [rad], longitude `lon` [rad], and altitude `h` [m],
considering the space indices in the named tuple `space_indices`, which must contain the
fields:

- `f107`: Daily 10.7 cm solar flux [sfu].
- `f107_avg`: 81-day average of the 10.7 cm solar flux [sfu].
- `kp`: Daily geomagnetic index Kp [-].
"""
function (::Jacchia77AtmosphericModel)(
    jd_utc::Number, lat::Number, lon::Number, h::Number, space_indices::NamedTuple
)
    # The Jacchia 1977 model is only valid between 90 km and 2000 km, whereas the
    # integrator can evaluate trial states outside this range (unphysical states near the
    # decay end or apogees above 2000 km). Hence, we clamp the altitude to keep the model
    # valid: the analysis terminates well above 90 km and the drag is negligible above
    # 2000 km, so the clamping does not change the result.
    h′ = clamp(h, 90.0e3, 2000.0e3)

    atmos = AtmosphericModels.jacchia1977(
        jd_utc,
        lat,
        lon,
        h′,
        space_indices.f107,
        space_indices.f107_avg,
        space_indices.kp
    )

    return atmos.total_density
end

"""
    struct Jr1971AtmosphericModel

Atmospheric model of the decay analysis wrapping the Jacchia-Roberts 1971 model provided
by **AtmosphericModels.jl**, selected by the macro
[`@decay_analysis__jr1971`](@ref SatelliteAnalysis.@decay_analysis__jr1971). It consumes
the space indices `f107` (daily 10.7 cm solar flux) [sfu], `f107_avg` (81-day average of
the 10.7 cm solar flux) [sfu], and `kp` (daily geomagnetic index Kp) [-] from the named
tuple passed to the callable.
"""
struct Jr1971AtmosphericModel end

"""
    (m::Jr1971AtmosphericModel)(
        jd_utc::Number,
        lat::Number,
        lon::Number,
        h::Number,
        space_indices::NamedTuple
    ) -> Float64

Compute the atmospheric density [kg/m³] using the Jacchia-Roberts 1971 model at the Julian
date `jd_utc` [UTC], geodetic latitude `lat` [rad], longitude `lon` [rad], and altitude
`h` [m], considering the space indices in the named tuple `space_indices`, which must
contain the fields:

- `f107`: Daily 10.7 cm solar flux [sfu].
- `f107_avg`: 81-day average of the 10.7 cm solar flux [sfu].
- `kp`: Daily geomagnetic index Kp [-].
"""
function (::Jr1971AtmosphericModel)(
    jd_utc::Number, lat::Number, lon::Number, h::Number, space_indices::NamedTuple
)
    # The Jacchia-Roberts 1971 model is only valid above 90 km, whereas the integrator can
    # evaluate trial states with unphysical altitudes near the decay end. Hence, we clamp
    # the altitude to keep the model valid: the analysis terminates well above 90 km, so
    # the clamping does not change the result.
    h′ = max(h, 90.0e3)

    atmos = AtmosphericModels.jr1971(
        jd_utc,
        lat,
        lon,
        h′,
        space_indices.f107,
        space_indices.f107_avg,
        space_indices.kp
    )

    return atmos.total_density
end
