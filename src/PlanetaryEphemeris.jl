module PlanetaryEphemeris

__precompile__(false)

export propagate, au, yr, sundofs, earthdofs,
    ssdofs, c_au_per_day, μ, NBP_pN_A_J23E_J23M_J2S!,
    semimajoraxis, eccentricity, inclination

using Reexport
@reexport using TaylorIntegration, LinearAlgebra # so that JLD may interpret previously saved Taylor1 objects saved in .jld files
@reexport using Printf
using Dates: DateTime, julian2datetime, datetime2julian
using DelimitedFiles
using Test
using JLD
using SPICE

# integration parameters
const varorder = 10
const order = 30
const abstol = 1.0E-30

const su = 1 #Sun's index
const ea = 4 #Earth's index
const mo = 5 #Moon's index

# vector of G*m values

const μ = [0.0002959122082855911, 4.91248045036476e-11, 7.24345233264412e-10, # Sun, Mercury, Venus
8.887692445125634e-10, 1.093189450742374e-11, 9.54954869555077e-11, # Earth, Moon, Mars
2.82534584083387e-7, 8.45970607324503e-8, 1.29202482578296e-8, # Jupiter, Saturn, Uranus
1.52435734788511e-8, # Neptune

2.17844105197418e-12, # Pluto

1.4004765561723440E-13, # 1 Ceres
3.8547501878088100E-14, # 4 Vesta
3.1044481989387130E-14, # 2 Pallas
1.2358007872941250E-14, # 10 Hygiea
6.3432804736486020E-15, # 31 Euphrosyne
5.2561686784936620E-15, # 704 Interamnia
5.1981269794574980E-15, # 511 Davida
4.6783074183509050E-15, # 15 Eunomia
3.6175383171479370E-15, # 3 Juno
3.4115868261938120E-15, # 16 Psyche
3.1806592826525410E-15, # 65 Cybele
2.5771141273110470E-15, # 88 Thisbe
2.5310917260150680E-15, # 48 Doris
2.4767881012558670E-15, # 52 Europa
2.2955593906374620E-15, # 451 Patientia
2.1992951735740730E-15, # 87 Sylvia

#TODO: ADD THE REST OF 343 MAIN-BELT ASTEROIDS
]

const N = length(μ)

# vector of J2*R^2 values
const Λ2 = zeros(N)
Λ2[su] = 4.5685187392703475e-12
Λ2[ea] = 1.9679542578489185e-12
Λ2[mo] = 2.7428745500623694e-14

# vector of J3*R^3 values
const Λ3 = zeros(N)
Λ3[ea] = -1.962633335678878e-19
Λ3[mo] = 1.3265639193531515e-20

# Matrix of J2 interactions included in DE430 ephemeris, according to Folkner et al., 2014
const UJ_interaction = fill(false, N, N)
UJ_interaction[2:end, su] .= true
UJ_interaction[union(1:ea-1,ea+1:7,N), ea] .= true
UJ_interaction[union(1:mo-1,mo+1:7), mo] .= true

const au = 1.495978707E8 # astronomical unit value in km
const yr = 365.25 # days in a Julian year
const daysec = 86_400 # number of seconds in a day
const c_au_per_day = daysec*(299_792.458/au) # speed of light in au per day
const c_au_per_sec = 299_792.458/au # speed of light in au per sec
const c_cm_per_sec = 100_000*299_792.458 # speed of light in cm per sec

const apophisdofs = union(3N-2:3N, 6N-2:6N) #union(34:36, 70:72)
const sundofs = union(1:3, 3(N+su)-2:3(N+su))
const earthdofs = union(3ea-2:3ea, 3(N+ea)-2:3(N+ea))
const ssdofs = setdiff(1:6N, apophisdofs)

const J2000 = 2.451545e6

const R_sun = 696000.0/au # Solar radius in au, value taken from DE430 docs

const α_p_sun = 268.13 # Sun's rotation pole right ascension (degrees)
const δ_p_sun = 63.87 # Sun's rotation pole declination (degrees)

function __init__()
    @show length(methods(NBP_pN_A_J23E_J23M_J2S!))
    @show length(methods(TaylorIntegration.jetcoeffs!))
    @show methods(NBP_pN_A_J23E_J23M_J2S!)
end

include("jpl-de-430-431-earth-orientation-model.jl")
include("initial_conditions.jl")
include("dynamical_model.jl")
include("propagation.jl")
include("osculating.jl")

end # module
