#Multi-threaded:
#JULIA_NUM_THREADS=<number-of-threads> julia --project=@. integrate_ephemeris.jl
#Single-threaded:
#julia --project=@. integrate_ephemeris.jl

@show Threads.nthreads()

using PlanetaryEphemeris
using Dates, JLD

#script parameters (TODO: use ArgParse.jl instead)
maxsteps = 1000000
nyears = 100.0
# jd0 = datetime2julian(DateTime(1969,6,28,0,0,0)) #starting time of integration
jd0 = datetime2julian(DateTime(2000,1,1,12)) #starting time of integration
@show jd0, J2000, jd0-J2000
dense = true #false
dynamics = DE430!
@show dynamics
nast = 343 #16 # number of asteroid perturbers
quadmath = false #true # use quadruple precision
###bodyind = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 17, 18, 19, 20, 21, 22, 23, 25, 27, 28, 30, 40, 41, 46, 55, 62, 73, 113, 115, 277, 322] # SS + 25 ast perturbers
bodyind = 1:(11+16) #1:(11+nast) # body indices in output

# integration parameters
order = 25
abstol = 1.0E-20

#integrator warmup
PlanetaryEphemeris.propagate(1, jd0, nyears, output=false, dense=dense, dynamics=dynamics, nast=nast, quadmath=quadmath, bodyind=bodyind, order=order, abstol=abstol)
println("*** Finished warmup")

# perform full integration
PlanetaryEphemeris.propagate(maxsteps, jd0, nyears, dense=dense, dynamics=dynamics, nast=nast, quadmath=quadmath, bodyind=bodyind, order=order, abstol=abstol)
println("*** Finished full integration")
