# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#   Compute the satellite lighting condition considering the Earth eclipse.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
# ==========================================================================================
#
# [1] Longo, C. R. O., Rickman, S. L (1995). Method for the Calculation of Spacecraft Umbra
#     and Penumbra Shadow Terminator Points. NASA Technical Paper 3547.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export lighting_condition

"""
    lighting_condition(r_i::AbstractVector, s_i::AbstractVector)

Compute the lighting condition at the position `r_i` [m] considering the Sun position vector
`s_i` [m]. The possible return values are:

- `:sunlight`: The point is under direct sunlight.
- `:penumbra`: The point is in penumbra region.
- `:umbra`: The point is in umbra region.

The algorithm used in this function was based on **[1]**.

!!! note
    The vectors `r_i` and `s_i` must be represented in the same reference frame.

# References

- **[1]** Longo, C. R. O., Rickman, S. L (1995). Method for the Calculation of Spacecraft
    Umbra and Penumbra Shadow Terminator Points. NASA Technical Paper 3547.
"""
function lighting_condition(r_i::AbstractVector, s_i::AbstractVector)
    # Constants.
    R₀ = EARTH_EQUATORIAL_RADIUS
    Rs = SUN_RADIUS

    # Distance from the Earth to the Sun.
    norm_s_i = norm(s_i)

    # Projection of the satellite position vector on the Sun direction [1].
    rs_i = (r_i ⋅ s_i) * s_i / (norm_s_i * norm_s_i)

    # Check if the satellite is under sunlight, umbra or penumbra.
    if (rs_i ⋅ s_i / norm_s_i) < 0
        # Distance of the umbral cone and the spacecraft [1].
        Δ_i = r_i - rs_i
        norm_Δ_i = norm(Δ_i)

        # Penumbra section [1].
        xp = R₀ * norm_s_i / (Rs + R₀)
        αp = asin(R₀ / xp)

        # Location of the penumbral cone terminator at the projected spacecraft
        # location [1].
        kp_i = (xp + norm(rs_i)) * tan(αp)

        if norm_Δ_i > kp_i
            return :sunlight
        else
            # Umbra section [1].
            xu = R₀ * norm_s_i / (Rs - R₀)
            αu = asin(R₀ / xu)

            # Location of the umbral cone terminator at the projected spacecraft
            # location [1].
            ep_i = (xu - norm(rs_i)) * tan(αu)

            if norm_Δ_i < ep_i
                return :umbra
            else
                return :penumbra
            end
        end
    else
        return :sunlight
    end
end
