#Multi-threaded:
#JULIA_NUM_THREADS=<number-of-threads> julia --project=@. main.jl
#Single-threaded:
#julia --project=@. main.jl

@show Threads.nthreads()

using PlanetaryEphemeris
using Dates, JLD

#script parameters (TODO: use ArgParse.jl instead)
const maxsteps = 10000
const nyears = 6.0 #21.0 #-5.0
const jd0 = datetime2julian(DateTime(2008,9,24,0,0,0)) #starting time of integration
@show jd0 == 2454733.5
const dense = true #false
const dynamics = DE430!
@show dynamics
const nast = 16 #343 # number of asteroid perturbers
const quadmath = false #true # use quadruple precision
const bodyind = 1:(11+16) # body indices in output

# path of lunar Euler angles ephemeris file
eulangfile = joinpath(pkgdir(PlanetaryEphemeris), "eph/moon_pa_de430_2003_2036_et.jld")

#integrator warmup
PlanetaryEphemeris.propagate(1, jd0, nyears, eulangfile, output=false, dense=dense, dynamics=dynamics, nast=nast, quadmath=quadmath, bodyind=bodyind)
println("*** Finished warmup")

PlanetaryEphemeris.propagate(2, jd0, nyears, eulangfile, output=false, dense=dense, dynamics=dynamics, nast=nast, quadmath=quadmath, bodyind=bodyind)
println("*** Finished 2 steps")

PlanetaryEphemeris.propagate(3, jd0, nyears, eulangfile, output=false, dense=dense, dynamics=dynamics, nast=nast, quadmath=quadmath, bodyind=bodyind)
println("*** Finished 3 steps")

PlanetaryEphemeris.propagate(3, jd0, nyears, eulangfile, dense=dense, dynamics=dynamics, nast=nast, quadmath=quadmath, bodyind=bodyind)
println("*** Finished 3 steps with output file")

#PlanetaryEphemeris.propagate(maxsteps, jd0, nyears, eulangfile, dense=dense, dynamics=dynamics, nast=nast, quadmath=quadmath, bodyind=bodyind)
#println("*** Finished full integration")
