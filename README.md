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
Minor Planet Center, Center for Astrophysics | Harvard & Smithsonian
- [Luis Benet](http://www.cicc.unam.mx/~benet/),
Instituto de Ciencias Físicas, Universidad Nacional Autónoma de México (UNAM)
- [Luis Eduardo Ramírez Montoya](https://github.com/LuEdRaMo),
Instituto de Ciencias Físicas, Universidad Nacional Autónoma de México (UNAM)

## Installation

The current version of this package may be installed in Julia pkg manager via:
```
] add PlanetaryEphemeris
```

## Usage

`PlanetaryEphemeris.propagate` is a high-level function which performs the
numerical integration. The file `de430.jl` in the `scripts` directory
contains an example script. This script may be called as:

`julia --project de430.jl --help`

`PlanetaryEphemeris.propagate` also supports multi-threading:

`julia -t <number-of-threads> --project de430.jl --help`

## Earth orientation models

`c2t_jpl_de430`/`t2c_jpl_de430` and `c2t_jpl_de440`/`t2c_jpl_de440` implement the Earth
orientation models of the JPL DE430/DE431 (Folkner et al., 2014) and DE440/DE441 (Park et
al., 2021) ephemerides. Only the direction of the Earth's pole is modeled (precession,
18.6-year nutation and the offsets of the EME2000 pole from the ICRF pole), not the
Earth's rotation. The DE440/DE441 model uses the long-term precession of Vondrák et al.
(2011); the cosine coefficient of Q_A with period 882 cy is 198.296701, as in SOFA,
instead of the 198.296071 printed in the paper's code (the pole is not affected).

## Known typos in Folkner et al. (2014)

The following values in Folkner et al. (2014) are typos, confirmed by the JPL team.
PlanetaryEphemeris.jl uses the correct values.

| Quantity | Location in the paper | Printed | Correct | Code |
|---|---|---|---|---|
| Rate of the lunar node longitude Ω | Eq. (21) | 1934°08′10″.549 | 1934°08′10″.539 | `coeffs_Ω_A` (`src/orientation.jl`) |
| Amplitude of the nutation in longitude Δψ | Eq. (22) | −17″.206262 | −17″.1996 | `coeff_Δψ` (`src/orientation.jl`) |
| Amplitude of the nutation in obliquity Δε | Eq. (22) | 9″.205348 | 9″.2025 | `coeff_Δϵ` (`src/orientation.jl`) |
| Right ascension of the Sun's pole | Table 9 | 268.13° | 286.13° | `α_p_sun` (`src/constants.jl`) |
| Earth's J4 | Table 10 | −1.616×10⁻⁶ | −1.619898×10⁻⁶ | `J4E` (`src/constants.jl`) |
| Degree of the Earth's zonal harmonics | Section III.B, Table 10 | 4 (no J5) | 5, J5 = −2.277345×10⁻⁷ | `J5E`, `n1SEM` (`src/constants.jl`) |

The correct values of Ω and of the nutation amplitudes are those of the IAU 1980 nutation
model, also printed in Park et al. (2021) for DE440/DE441.

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
- Park, R. S., Folkner, W. M., Williams, J. G., & Boggs, D. H. (2021). The JPL
  planetary and lunar ephemerides DE440 and DE441. The Astronomical Journal,
  161(3), 105. https://doi.org/10.3847/1538-3881/abd414
- Vondrák, J., Capitaine, N., & Wallace, P. (2011). New precession expressions,
  valid for long time intervals. Astronomy & Astrophysics, 534, A22.
  https://doi.org/10.1051/0004-6361/201117274
