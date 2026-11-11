# This file is part of the PlanetaryEphemeris.jl package; MIT licensed

@doc """
    ordpres_differentiate(a::Taylor1)

Returns the derivative of `a`, but preserving the order/degree of `a`. In comparison,
`TaylorSeries.differentiate` returns the returns a `Taylor1` object with one order/degree
less than the one of `a`.

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
    loadeph(ss16asteph::TaylorSolution, μ::Vector{<:Real})

Taking Solar System ephemeris `ss16asteph` and their gravitational parameters `μ` as input, returns for all bodies the point-mass Newtonian acceleration and the Newtonian N body potential.

# Arguments

- `ss16asteph`: Solar System ephemeris.
- `μ::Vector{<:Real}`: vector of mass parameters.
"""
function loadeph(ss16asteph::TaylorSolution, μ::Vector{<:Real})

    # Compute point-mass Newtonian accelerations from ephemeris
    # accelerations of all bodies are needed to compute the post-Newtonian acceleration of e.g. Solar System minor bodies
    # Number of bodies that contibute to the asteroid's acceleration
    Nm1 = numberofbodies(ss16asteph)
    # Initialize dense polynomial arrays for the point-mass Newtonian accelerations
    acc_p = Matrix{eltype(ss16asteph.p)}(undef, length(ss16asteph.t)-1, 3Nm1)
    # Initialize dense polynomial arrays for the newtonian N body potential
    pot_p = Matrix{eltype(ss16asteph.p)}(undef, length(ss16asteph.t)-1, Nm1)
    fill!(acc_p, zero(ss16asteph.p[1]))
    fill!(pot_p, zero(ss16asteph.p[1]))

    # Iterator over all bodies except asteroid
    for j in 1:Nm1
        for i in 1:Nm1
            if i == j
                #
            else
                # Difference between two positions (\mathbf{r}_i - \mathbf{r}_j)
                X_ij = ss16asteph.p[:, 3i-2] .- ss16asteph.p[:, 3j-2]  # X-axis component
                Y_ij = ss16asteph.p[:, 3i-1] .- ss16asteph.p[:, 3j-1]  # Y-axis component
                Z_ij = ss16asteph.p[:, 3i  ] .- ss16asteph.p[:, 3j  ]  # Z-axis component
                # Distance between two bodies squared ||\mathbf{r}_i - \mathbf{r}_j||^2
                r_p2_ij = ( (X_ij.^2) .+ (Y_ij.^2) ) .+ (Z_ij.^2)
                # Distance between two bodies ||\mathbf{r}_i - \mathbf{r}_j||
                r_ij = sqrt.(r_p2_ij)
                # Newtonian potential
                pot_p[:, j] .+= (μ[i]./r_ij)
            end
        end

        # Fill acelerations by differentiating velocities
        acc_p[:, 3j-2] .= ordpres_differentiate.(ss16asteph.p[:, 3(Nm1+j)-2])  # X-axis component
        acc_p[:, 3j-1] .= ordpres_differentiate.(ss16asteph.p[:, 3(Nm1+j)-1])  # Y-axis component
        acc_p[:, 3j  ] .= ordpres_differentiate.(ss16asteph.p[:, 3(Nm1+j)  ])  # Z-axis component
    end

    return TaylorSolution(ss16asteph.t, acc_p), TaylorSolution(ss16asteph.t, pot_p)
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
