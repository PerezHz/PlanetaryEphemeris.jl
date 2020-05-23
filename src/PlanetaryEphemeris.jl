module PlanetaryEphemeris

__precompile__(false)

export au, yr, sundofs, earthdofs,
    c_au_per_day, μ, NBP_pN_A_J23E_J23M_J2S!,
    NBP_pN_A_J23E_J23M_J2S_threads!,
    semimajoraxis, eccentricity, inclination,
    longascnode, argperi, longperi,
    trueanomaly, ecanomaly, meananomaly,
    timeperipass, lrlvec, eccentricanomaly,
    meanan2truean, meanmotion, time2truean,
    su, ea, mo, au, yr, daysec, clightkms,
    c_au_per_day, c_au_per_sec, c_cm_per_sec,
    J2000, R_sun, α_p_sun, δ_p_sun, au,
    UJ_interaction, de430_343ast_ids, Rx, Ry, Rz,
    ITM_und, ITM1, ITM2, R_moon, τ_M, k_2M,
    JSEM, CM, SM, n1SEM, n2M, J2E, J2EDOT, RE

using TaylorIntegration, LinearAlgebra
using Printf
using Dates: DateTime, julian2datetime, datetime2julian
using DelimitedFiles
using Test
using JLD
using Quadmath

# integration parameters
const order = 30
const abstol = 1.0E-30

include("constants.jl")
include("jpl-de-430-431-earth-orientation-model.jl")
include("initial_conditions.jl")
include("dynamical_model.jl")
include("propagation.jl")
include("osculating.jl")

function __init__()
    @show length(methods(TaylorIntegration.jetcoeffs!))
    @show methods(NBP_pN_A_J23E_J23M_J2S!)
    @show methods(NBP_pN_A_J23E_J23M_J2S_threads!)
    @show methods(TaylorIntegration.jetcoeffs!)
end

end # module
