module PlanetaryEphemeris

# __precompile__(false)

export PE, au, yr, sundofs, earthdofs,
    c_au_per_day, μ, NBP_pN_A_J23E_J23M_J2S!,
    NBP_pN_A_J23E_J23M_J2S_threads!, DE430!,
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
    JSEM, CM, SM, n1SEM, n2M, J2E, J2EDOT, RE,
    k_20E, k_21E, k_22E, τ_0p, τ_1p, τ_2p, τ_0, τ_1, τ_2, ω_E, EMRAT,
    TaylorInterpolant, selecteph2jld, ssb_posvel_pN, nbodyind, propagate, propagate_dense

using AutoHashEquals
using TaylorIntegration, LinearAlgebra
using Printf
using Dates: DateTime, julian2datetime, datetime2julian, year
using DelimitedFiles
using JLD
using Quadmath

import Base.reverse

@doc raw"""
    nbodyind(N::Int, i::Int)
    nbodyind(N::Int, ivec::AbstractVector{Int})

Returns the indexes of the positions and velocities of the `i`-th body (or the 
`ivec`-th bodies) in a vector with `N` bodies. The function assumes that the vector has 
the form: `3N` positions + `3N` velocities (+ Lunar physical librations + TT-TDB). 
"""
nbodyind(N::Int, i::Int) = union(3i-2:3i, 3*(N+i)-2:3*(N+i))

function nbodyind(N::Int, ivec::AbstractVector{Int})
    a = Int[]
    for i in ivec
        i > N && continue
        a = union(a, nbodyind(N,i))
    end
    return sort(a)
end

include("constants.jl")
include("jpl-de-430-431-earth-orientation-model.jl")
include("initial_conditions.jl")
include("dynamical_model.jl")
include("jetcoeffs.jl")
include("interpolation.jl")
include("plephinteg.jl")
include("propagation.jl")
include("osculating.jl")
include("barycenter.jl")

end # module
