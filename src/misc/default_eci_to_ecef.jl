## Description #############################################################################
#
# Default function to convert vectors from the ECI to the ECEF reference frame.
#
############################################################################################

############################################################################################
#                                    Private Functions                                     #
############################################################################################

"""
    _default_eci_to_ecef(r_eci::AbstractVector, jd::Number) -> AbstractVector

Convert the vector `r_eci` from the Earth-centered inertial (ECI) reference frame to the
Earth-centered, Earth-fixed (ECEF) reference frame at the instant `jd` [Julian Day, UTC].
This function considers TEME as the ECI reference frame and PEF as the ECEF reference frame.
It is the default value of the keyword `f_eci_to_ecef` in the analyses, and it is
thread-safe.
"""
function _default_eci_to_ecef(r_eci::AbstractVector, jd::Number)
    D_ecef_eci = r_eci_to_ecef(TEME(), PEF(), jd)
    return D_ecef_eci * r_eci
end
