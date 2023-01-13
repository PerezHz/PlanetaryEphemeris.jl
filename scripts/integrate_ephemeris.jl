#Multi-threaded:
# julia -t <number-of-threads> --project=@. integrate_ephemeris.jl
# Single-threaded:
# julia --project=@. integrate_ephemeris.jl

println("Number of threads: ", Threads.nthreads())

using PlanetaryEphemeris
using Dates

#script parameters (TODO: use ArgParse.jl instead)
const maxsteps = 100 # 1000000
# jd0 = datetime2julian(DateTime(1969,6,28,0,0,0)) #starting time of integration
const jd0 = datetime2julian(DateTime(2000,1,1,12)) #starting time of integration
const nyears = 2031.0 - year(julian2datetime(jd0))
const dynamics = DE430!
println("Dynamical function: ", dynamics)
const nast = 343 #16 # number of asteroid perturbers
const quadmath = false #true # use quadruple precision
###bodyind = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 17, 18, 19, 20, 21, 22, 23, 25, 27, 28, 30, 40, 41, 46, 55, 62, 73, 113, 115, 277, 322] # SS + 25 ast perturbers
const bodyind = 1:(11+16) #1:(11+nast) # body indices in output

# Integration parameters
const order = 25
const abstol = 1.0E-20

# Integrator warmup
propagate_dense(1, jd0, nyears, Val(quadmath); output = false, dynamics=dynamics, nast=nast, bodyind=bodyind, order=order, abstol=abstol)
println("*** Finished warmup")

# Perform full integration
propagate_dense(maxsteps, jd0, nyears, Val(quadmath); dynamics=dynamics, nast=nast, bodyind=bodyind, order=order, abstol=abstol)
println("*** Finished full integration")
