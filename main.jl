#Multi-threaded:
#JULIA_NUM_THREADS=<number-of-threads> julia --project=@. main.jl
#Single-threaded:
#julia --project=@. main.jl
using PlanetaryEphemeris
using Dates, JLD

#script parameters (TODO: use ArgParse.jl instead)
const maxsteps = 10000
const nyears = 5.0 #24.0
const t0 = datetime2julian(DateTime(2008,9,24,0,0,0)) #starting time of integration
@show t0 == 2454733.5
const dense = true#false
# const dynamics = NBP_pN_A_J23E_J23M_J2S!
const dynamics = NBP_pN_A_J23E_J23M_J2S_threads!
@show dynamics
const nast = 16 # number of asteroid perturbers

# path of lunar Euler angles ephemeris file
eulangfile = joinpath(dirname(pathof(PlanetaryEphemeris)), "../eph/moon_pa_de430_tx_2008_2013.jld")

#integrator warmup
PlanetaryEphemeris.propagate(1, t0, nyears, eulangfile, output=false, dense=dense, dynamics=dynamics, nast=nast)
println("*** Finished warmup")

PlanetaryEphemeris.propagate(2, t0, nyears, eulangfile, output=false, dense=dense, dynamics=dynamics, nast=nast)
println("*** Finished 2 steps")

PlanetaryEphemeris.propagate(3, t0, nyears, eulangfile, output=false, dense=dense, dynamics=dynamics, nast=nast)
println("*** Finished 3 steps")

PlanetaryEphemeris.propagate(3, t0, nyears, eulangfile, dense=dense, dynamics=dynamics, nast=nast)
println("*** Finished 3 steps with output file")

#PlanetaryEphemeris.propagate(maxsteps, t0, nyears, eulangfile, dense=dense, dynamics=dynamics, nast=nast)
#println("*** Finished full integration")
