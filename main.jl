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
const jd0 = datetime2julian(DateTime(1969,6,28,0,0,0)) #starting time of integration
@show jd0, J2000, jd0-J2000
const dense = true #false
const dynamics = DE430!
@show dynamics
const nast = 16 #343 #25 # number of asteroid perturbers
const quadmath = false #true # use quadruple precision
###const bodyind = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 17, 18, 19, 20, 21, 22, 23, 25, 27, 28, 30, 40, 41, 46, 55, 62, 73, 113, 115, 277, 322] # SS + 25 ast perturbers
const bodyind = 1:(11+nast) # body indices in output

# integration parameters
const order = 25
const abstol = 1.0E-20

#integrator warmup
PlanetaryEphemeris.propagate(1, jd0, nyears, output=false, dense=dense, dynamics=dynamics, nast=nast, quadmath=quadmath, bodyind=bodyind, order=order, abstol=abstol)
println("*** Finished warmup")

# PlanetaryEphemeris.propagate(2, jd0, nyears, output=false, dense=dense, dynamics=dynamics, nast=nast, quadmath=quadmath, bodyind=bodyind, order=order, abstol=abstol)
# println("*** Finished 2 steps")

# PlanetaryEphemeris.propagate(3, jd0, nyears, output=false, dense=dense, dynamics=dynamics, nast=nast, quadmath=quadmath, bodyind=bodyind, order=order, abstol=abstol)
# println("*** Finished 3 steps")

PlanetaryEphemeris.propagate(3, jd0, nyears, dense=dense, dynamics=dynamics, nast=nast, quadmath=quadmath, bodyind=bodyind, order=order, abstol=abstol)
println("*** Finished 3 steps with output file")

#PlanetaryEphemeris.propagate(maxsteps, jd0, nyears, dense=dense, dynamics=dynamics, nast=nast, quadmath=quadmath, bodyind=bodyind, order=order, abstol=abstol)
#println("*** Finished full integration")
