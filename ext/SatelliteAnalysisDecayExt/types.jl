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

