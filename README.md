# PlanetaryEphemeris.jl

`PlanetaryEphemeris.jl` is a Taylor integrator for the JPL DE430 planetary
ephemeris dynamical model (Folkner et al., 2014), based on
[TaylorIntegration.jl](https://github.com/PerezHzTaylorIntegration.jl).

## Authors

- [Jorge A. Pérez](https://www.linkedin.com/in/perezhz),
Instituto de Ciencias Físicas, Universidad Nacional Autónoma de México (UNAM)
- [Luis Benet](http://www.cicc.unam.mx/~benet/),
Instituto de Ciencias Físicas, Universidad Nacional Autónoma de México (UNAM)

## Installation

The current development version of this package may be installed in Julia via:
```
import Pkg
Pkg.add(Pkg.PackageSpec(url="https://github.com/PerezHz/PlanetaryEphemeris.jl.git", rev="jp/threads"))
```

## Usage

`PlanetaryEphemeris.propagate` is a high-level function which performs the
numerical integration. The file `main.jl` in this package root directory
contains an example script. This script may be called from any subfolder of this
package simply as:

`julia --project=@. main.jl`

`PlanetaryEphemeris.propagate` also supports multi-threading:

`JULIA_NUM_THREADS=<number-of-threads> julia --project=@. main.jl`

## Acknowledgments

We acknowledge financial support from UNAM-PAPIIT grant IG100819 and computing
resources provided by LANCAD-UNAM-DGTIC-284.

## References

- Pérez-Hernández, Jorge A., & Benet, Luis. (2020, October 13).
    [PerezHz/TaylorIntegration.jl](https://github.com/PerezHzTaylorIntegration.jl):
    v0.8.4 (Version v0.8.4). Zenodo. http://doi.org/10.5281/zenodo.4086166
- Folkner, W. M., Williams, J. G., Boggs, D. H., Park, R. S., & Kuchynka, P.
  (2014). The planetary and lunar ephemerides DE430 and DE431. Interplanetary
  Network Progress Report, 196(1).
