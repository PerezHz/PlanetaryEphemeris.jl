# PlanetaryEphemeris.jl

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5152451.svg)](https://doi.org/10.5281/zenodo.5152451)
[![CI](https://github.com/PerezHz/PlanetaryEphemeris.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/PerezHz/PlanetaryEphemeris.jl/actions/workflows/CI.yml)
[![Coverage Status](https://coveralls.io/repos/github/PerezHz/PlanetaryEphemeris.jl/badge.svg?branch=main)](https://coveralls.io/github/PerezHz/PlanetaryEphemeris.jl?branch=main)
[![codecov](https://codecov.io/gh/PerezHz/PlanetaryEphemeris.jl/branch/main/graph/badge.svg?token=CZE0SONYYX)](https://codecov.io/gh/PerezHz/PlanetaryEphemeris.jl)

`PlanetaryEphemeris.jl` is a Taylor integrator for the JPL DE430 planetary
ephemeris dynamical model (Folkner et al., 2014), based on
[TaylorIntegration.jl](https://github.com/PerezHzTaylorIntegration.jl).

## Authors

- [Jorge A. Pérez Hernández](https://github.com/PerezHz),
Instituto de Ciencias Físicas, Universidad Nacional Autónoma de México (UNAM)
- [Luis Benet](http://www.cicc.unam.mx/~benet/),
Instituto de Ciencias Físicas, Universidad Nacional Autónoma de México (UNAM)
- [Luis Eduardo Ramírez Montoya](https://github.com/LuEdRaMo),
Facultad de Ciencias, Universidad Nacional Autónoma de México (UNAM)

## Installation

The current version of this package may be installed in Julia pkg manager via:
```
] add PlanetaryEphemeris
```

## Usage

`PlanetaryEphemeris.propagate` is a high-level function which performs the
numerical integration. The file `integrate_ephemeris.jl` in the `scripts` directory
contains an example script. This script may be called as:

`julia --project integrate_ephemeris.jl --help`

`PlanetaryEphemeris.propagate` also supports multi-threading:

`julia -t <number-of-threads> --project integrate_ephemeris.jl --help`

## Acknowledgments

We acknowledge financial support from UNAM-PAPIIT grants IG-100819 and IG-101122, as well as
computing resources provided by LANCAD-UNAM-DGTIC-284.

## References

- Pérez-Hernández, Jorge A., & Benet, Luis. (2020, October 13).
    [PerezHz/TaylorIntegration.jl](https://github.com/PerezHzTaylorIntegration.jl):
    v0.8.4 (Version v0.8.4). Zenodo. http://doi.org/10.5281/zenodo.4086166
- Folkner, W. M., Williams, J. G., Boggs, D. H., Park, R. S., & Kuchynka, P.
  (2014). The planetary and lunar ephemerides DE430 and DE431. Interplanetary
  Network Progress Report, 196(1).
