# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==============================================================================
#
#   Compute the satellite lighting condition considering the Earth eclipse.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
# ==============================================================================
#
# [1] Longo, C. R. O., Rickman, S. L (1995). Method for the Calculation of
#     Spacecraft Umbra and Penumbra Shadow Terminator Points. NASA Technical
#     Paper 3547.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export lighting_condition

"""
    lighting_condition(r_i::AbstractVector, s_i::AbstractVector)

Compute the lighting condition at the position `r_i` [m] considering the Sun
position vector `s_i` [m]. The possible return values are:

- `:sunlight`: The point is under direct sunlight.
- `:penumbra`: The point is in penumbra region.
- `:umbra`: The point is in umbra region.

!!! note
    The vectors `r_i` and `s_i` must be represented in the same reference frame.
"""
function lighting_condition(r_i::AbstractVector, s_i::AbstractVector)
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
        xp = R0 * norm_s_i / (Rs + R0)
        αp = asin(R0 / xp)

        # Location of the penumbral cone terminator at the projected spacecraft
        # location [1].
        kp_i = (xp + norm(rs_i)) * tan(αp)

        if norm_Δ_i > kp_i
            return :sunlight
        else
            # Umbra section [1].
            xu = R0 * norm_s_i / (Rs - R0)
            αu = asin(R0 / xu)

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
