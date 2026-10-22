# This file is part of the PlanetaryEphemeris.jl package; MIT licensed

@doc """
    ordpres_differentiate(a::Taylor1)

Returns the derivative of `a`, but preserving the order/degree of `a`. In
comparison, `TaylorSeries.differentiate` returns a `Taylor1` object with one
order/degree less than the one of `a`.

See also [`TaylorSeries.differentiate`](@ref).
"""
function ordpres_differentiate(a::Taylor1{T}) where {T}
    res = zero(a)
    for ord in eachindex(res)
        TaylorSeries.differentiate!(res, a, ord)
    end
    return res
end

@doc raw"""
    accpoteph(sseph::TaylorSolution, GM::Vector{<:Real})

Return the accelerations and the Newtonian N-body potentials of all the bodies in the
Solar System ephemeris `sseph`, as two dense `TaylorSolution`s, `acceph` and `poteph`.
The accelerations are obtained by differentiating the velocities in `sseph`.

# Arguments

- `sseph::TaylorSolution`: Solar System ephemeris.
- `GM::Vector{<:Real}`: mass parameters of the bodies in `sseph` [au^3/day^2].
"""
function accpoteph(sseph::TaylorSolution, GM::Vector{<:Real})

    # Compute the accelerations and Newtonian potentials of all bodies in the ephemeris,
    # which are needed to compute the post-Newtonian acceleration of e.g. Solar System
    # minor bodies
    # Number of bodies that contribute to the minor body's acceleration
    Nm1 = numberofbodies(sseph)
    # Initialize dense polynomial arrays for the accelerations
    acceph = Matrix{eltype(sseph.p)}(undef, length(sseph.t)-1, 3Nm1)
    # Initialize dense polynomial arrays for the Newtonian N-body potential
    poteph = Matrix{eltype(sseph.p)}(undef, length(sseph.t)-1, Nm1)
    fill!(acceph, zero(sseph.p[1]))
    fill!(poteph, zero(sseph.p[1]))

    # Iterator over all bodies in the ephemeris
    for j in 1:Nm1
        for i in 1:Nm1
            if i == j
                #
            else
                # Difference between two positions (\mathbf{r}_i - \mathbf{r}_j)
                X_ij = sseph.p[:, 3i-2] .- sseph.p[:, 3j-2]  # X-axis component
                Y_ij = sseph.p[:, 3i-1] .- sseph.p[:, 3j-1]  # Y-axis component
                Z_ij = sseph.p[:, 3i  ] .- sseph.p[:, 3j  ]  # Z-axis component
                # Distance between two bodies squared ||\mathbf{r}_i - \mathbf{r}_j||^2
                r_p2_ij = ( (X_ij.^2) .+ (Y_ij.^2) ) .+ (Z_ij.^2)
                # Distance between two bodies ||\mathbf{r}_i - \mathbf{r}_j||
                r_ij = sqrt.(r_p2_ij)
                # Newtonian potential
                poteph[:, j] .+= (GM[i]./r_ij)
            end
        end

        # Fill accelerations by differentiating velocities (X, Y and Z-axis components)
        acceph[:, 3j-2] .= ordpres_differentiate.(sseph.p[:, 3(Nm1+j)-2])
        acceph[:, 3j-1] .= ordpres_differentiate.(sseph.p[:, 3(Nm1+j)-1])
        acceph[:, 3j  ] .= ordpres_differentiate.(sseph.p[:, 3(Nm1+j)  ])
    end

    return TaylorSolution(sseph.t, acceph), TaylorSolution(sseph.t, poteph)
end

"""
    propagate(PE; kwargs...)

Propagate a planetary ephemeris problem `PE` using the Taylor method
implemented in `TaylorIntegration`.

# Keyword arguments

- `maxsteps::Int`: maximum number of steps for the integration (default: `500`).
- `order::Int`: order of Taylor expansions wrt time (default: 25).
- `abstol::T`: absolute tolerance used to compute the propagation timestep
    (default: `1E-20`).
- `parse_eqs::Bool`: whether to use the specialized method of `jetcoeffs` or not
    (default: `true`).
"""
function propagate(PE::PlanetaryEphemerisProblem{D, T, P};
                   maxsteps::Int = 500, order::Int = _order,
                   abstol::T = abstol, parse_eqs::Bool = true) where {D, T, P}
    # Unpack
    @unpack dynamics, tspan, initcond, params = PE
    # Integration
    sol = taylorinteg(dynamics, initcond, zero(T), tspan[2] - tspan[1], order,
                      abstol, params; maxsteps, parse_eqs)
    return TaylorSolution(collect((tspan[1] - J2000) .+ sol.t), collect(sol.p))
end
