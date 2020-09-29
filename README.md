# PlanetaryEphemeris.jl

A `TaylorIntegration.jl`-based Taylor integrator for the JPL DE430 planetary
ephemeris dynamical model (Folkner et al., 2014).

# Usage

`PlanetaryEphemeris.propagate` is a high-level function which performs the
numerical integration. The file `main.jl` in this package root directory
contains an example script. This script may be called from any subfolder of this
package simply as:

`julia --project=@. main.jl`

`PlanetaryEphemeris.propagate` also supports multi-threading:

`JULIA_NUM_THREADS=<number-of-threads> julia --project=@. main.jl`.