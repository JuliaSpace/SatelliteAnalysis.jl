module SatelliteAnalysisDecayExt

using LinearAlgebra
using OrdinaryDiffEqAdamsBashforthMoulton
using PrecompileTools
using SatelliteAnalysis
using StaticArrays

############################################################################################
#                                        Constants                                         #
############################################################################################

# Solar radiation pressure at 1 AU [N / m²].
const _SOLAR_PRESSURE_1AU = 1367 / 299792458

# Sun gravitational parameter [m³ / s²].
const _μ_SUN  = 1.32712440018e20

# Moon gravitational parameter [m³ / s²].
const _μ_MOON = 4.9027988e12

############################################################################################
#                                          Types                                           #
############################################################################################

include("types.jl")

############################################################################################
#                                         Includes                                         #
############################################################################################

include("accelerations.jl")
include("dynamics.jl")
include("gauss_matrices.jl")
include("misc.jl")
include("progress.jl")
include("variational_rates.jl")
include("main.jl")
include("precompile.jl")

end

