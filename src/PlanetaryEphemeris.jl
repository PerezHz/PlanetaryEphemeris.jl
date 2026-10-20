module PlanetaryEphemeris

# __precompile__(false)

export PlanetaryEphemerisProblem
export PE, au, yr, su, ea, mo, daysec, clightkms, c_au_per_day, c_au_per_sec,
       c_cm_per_sec, J2000, ω_E, kmsec2auday, auday2kmsec, numberofbodies,
       ssb_posvel_pN, nbodyind, read_initial_conditions, read_parameters
export Rx, Ry, Rz, ITM1, ITM2, pole_rotation
export DE430Constants, DE430Params, NBP_pN_A_J23E_J23M_J2S_threads!, DE430!,
       t2c_jpl_de430, c2t_jpl_de430
export DE440Constants, DE440Params, t2c_jpl_de440, c2t_jpl_de440
export semimajoraxis, eccentricity, inclination, longascnode, argperi, longperi, ecanomaly,
       trueanomaly,  meananomaly, timeperipass, lrlvec, eccentricanomaly, meanan2truean,
       meanmotion, time2truean
export propagate, selecteph2jld2, save2jld2andcheck, selecteph

using TaylorIntegration, LinearAlgebra, Printf, DelimitedFiles, JLD2, Quadmath

using Dates: DateTime, datetime2julian, year
using Parameters: @unpack
using TaylorIntegration: RetAlloc, _determine_parsing!, init_expansions
using TaylorSeries: numtype, matmul!

import Base: show
import Dates: julian2datetime
import TaylorSeries

const TS = TaylorSeries

include("abstractproblem.jl")
include("constants.jl")
include("orientation.jl")
include("initcond.jl")
include("parameters.jl")
include("interpolation/TaylorSolution.jl")
include("propagation.jl")
include("osculating.jl")
include("barycenter.jl")
#include("precompile.jl")

include("trivial.jl")
include("DE430/params.jl")
include("DE430/orientation.jl")
include("DE430/dynamicalmodel.jl")
include("DE430/jetcoeffs.jl")
include("DE440/params.jl")
include("DE440/orientation.jl")

end # module
