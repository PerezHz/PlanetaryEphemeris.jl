# This file is part of the PlanetaryEphemeris.jl package; MIT licensed

"""
    read_initial_conditions(filename)

Read a vector of initial conditions from a file.

# Extended help

For a planetary ephemeris model with `N` bodies, the vector of initial conditions
consists of `3N` positions [au], `3N` velocities [au/day], 3 lunar mantle Euler
angles [rad], 3 lunar mantle angular velocities [rad/day], 3 lunar core Euler
angles [rad], 3 lunar core angular velocities [rad/day] and the initial value
of TT-TDB [sec].
"""
function read_initial_conditions(filename::AbstractString)
    # Read lines except header
    lines = readlines(filename)
    popfirst!(lines)
    # Number of bodies
    N = length(lines) - 3
    # Allocate the initial conditions vector
    q = Vector{Float64}(undef, 6N+13)
    # Planets + asteroids
    for i in 1:N
        q[nbodyind(N, i)] .= parse.(Float64, Iterators.drop(eachsplit(lines[i]), 2))
    end
    # Lunar mantle
    q[6N+1:6N+6] .= parse.(Float64, Iterators.drop(eachsplit(lines[end-2]), 2))
    # Lunar core
    q[6N+7:6N+12] .= parse.(Float64, Iterators.drop(eachsplit(lines[end-1]), 2))
    # TT-TDB
    q[6N+13] = parse(Float64, split(lines[end])[3])

    return q
end