#julia --project=@. main.jl
using PlanetaryEphemeris
using Dates

#script parameters (TODO: use ArgParse.jl instead)
const maxsteps = 10000
const nyears = 5.0 # since t0 is 2008-9-24, this ends the integration on 2032-9-24; NOTE: this value is overriden when evaluating solution at JPL radar observation times
const t0 = datetime2julian(DateTime(2008,9,24,0,0,0)) #starting time of integration
@show t0 == 2454733.5
const dense = true#false

#integrator warmup
PlanetaryEphemeris.propagate(1, t0, nyears, output=false, dense=dense)
println("*** Finished warmup")

PlanetaryEphemeris.propagate(2, t0, nyears, dense=dense)
println("*** Finished 2nd warmup")

#PlanetaryEphemeris.propagate(maxsteps, t0, nyears, dense=dense)
#println("*** Finished full integration")
