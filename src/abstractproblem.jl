"""
    AbstractPlanetaryEphemerisProblem

Supertype for the planetary ephemeris problems interface.
"""
abstract type AbstractPlanetaryEphemerisProblem end

"""
    PlanetaryEphemerisProblem{D, T, P} <: AbstractPlanetaryEphemerisProblem

A planetary ephemeris problem.

# Fields

- `dynamics::D`: dynamical model.
- `epoch::T`: reference epoch [JDTDB].
- `initcond::Vector{T}`: initial condition.
- `params::P`: parameters.
"""
struct PlanetaryEphemerisProblem{D, T, P} <: AbstractPlanetaryEphemerisProblem
    dynamics::D
    epoch::T
    initcond::Vector{T}
    params::P
end

# Print method for PlanetaryEphemerisProblem
function show(io::IO, x::PlanetaryEphemerisProblem)
    t = repeat(' ', 4)
    print(io,
        "Planetary ephemeris problem\n",
        t, rpad("Dynamical model:", 21), x.dynamics, "\n",
        t, rpad("Reference epoch:", 21), julian2datetime(x.epoch), "\n",
        t, rpad("Initial condition:", 21), length(x.initcond), " degrees of freedom", "\n",
        t, rpad("Parameters:", 21), typeof(x.params), "\n"
    )
end