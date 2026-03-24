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
- `tspan::Tuple{T, T}`: time span [JDTDB].
- `initcond::Vector{T}`: initial condition.
- `params::P`: parameters.
"""
struct PlanetaryEphemerisProblem{D, T, P} <: AbstractPlanetaryEphemerisProblem
    dynamics::D
    tspan::Tuple{T, T}
    initcond::Vector{T}
    params::P
end

# Print method for PlanetaryEphemerisProblem
function show(io::IO, x::PlanetaryEphemerisProblem)
    tab = repeat(' ', 4)
    t0, tf = @. string(julian2datetime(x.tspan))
    print(io,
        "Planetary ephemeris problem\n",
        tab, rpad("Dynamical model:", 21), x.dynamics, "\n",
        tab, rpad("Time span:", 21), (t0, tf), "\n",
        tab, rpad("Initial condition:", 21), length(x.initcond), " degrees of freedom", "\n",
        tab, rpad("Parameters:", 21), typeof(x.params), "\n"
    )
end