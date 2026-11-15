# This file is part of the PlanetaryEphemeris.jl package; MIT licensed

@doc """
    ordpres_differentiate(a::Taylor1)

Returns the derivative of `a`, but preserving the order/degree of `a`. In comparison,
`TaylorSeries.differentiate` returns a `Taylor1` object with one order/degree
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
    accpoteph(sseph::TaylorSolution, GM::AbstractVector)

Return the point-mass Newtonian accelerations and the Newtonian ``N``-body potentials of
all the bodies in the Solar System ephemeris `sseph`, given their gravitational parameters
`GM`.

# Arguments

- `sseph::TaylorSolution`: Solar System ephemeris (with real Taylor coefficients).
- `GM::AbstractVector`: vector of gravitational parameters.
"""
function accpoteph(sseph::TaylorSolution{T, U},
                   GM::AbstractVector{T}) where {T <: Real, U <: Real}
    # Number of bodies
    Nm1 = numberofbodies(sseph)
    # Number of time steps and order of the Taylor polynomials
    nsteps = size(sseph.p, 1)
    order = TS.order(sseph.p[1])
    # Point-mass Newtonian accelerations (X, Y and Z components of each body) and
    # Newtonian N body potentials
    acceph = [Taylor1(zeros(U, order+1)) for _ in 1:nsteps, _ in 1:3Nm1]
    poteph = [Taylor1(zeros(U, order+1)) for _ in 1:nsteps, _ in 1:Nm1]
    # Auxiliaries, reused across the whole computation
    X, Y, Z = Taylor1(zeros(U, order+1)), Taylor1(zeros(U, order+1)),
        Taylor1(zeros(U, order+1))          # Difference between two positions
    X2, Y2, Z2 = Taylor1(zeros(U, order+1)), Taylor1(zeros(U, order+1)),
        Taylor1(zeros(U, order+1))          # Squares of the components
    X2pY2 = Taylor1(zeros(U, order+1))      # X^2 + Y^2
    r2 = Taylor1(zeros(U, order+1))         # Distance squared
    r = Taylor1(zeros(U, order+1))          # Distance
    GMi_r = Taylor1(zeros(U, order+1))      # GM[i] / r
    GMj_r = Taylor1(zeros(U, order+1))      # GM[j] / r
    aux = Taylor1(zeros(U, order+1))        # Auxiliary for sqrt!
    aux0 = zero(U)                          # Auxiliary for sqr!

    for k in 1:nsteps
        for j in 1:Nm1
            # Accelerations: derivatives of the velocities (the last coefficient stays 0)
            for c in 1:3
                vj = sseph.p[k, 3(Nm1+j)-3+c]
                for ord in 0:order-1
                    TS.differentiate!(acceph[k, 3j-3+c], vj, ord)
                end
            end
            # Newtonian potentials: each pair (i, j) with i > j contributes GM[i] / r to
            # the potential of j and GM[j] / r to the potential of i. Terms are added in
            # increasing order of the other body's index
            for i in j+1:Nm1
                for ord in 0:order
                    # Difference between two positions (\mathbf{r}_i - \mathbf{r}_j)
                    TS.subst!(X, sseph.p[k, 3i-2], sseph.p[k, 3j-2], ord)
                    TS.subst!(Y, sseph.p[k, 3i-1], sseph.p[k, 3j-1], ord)
                    TS.subst!(Z, sseph.p[k, 3i  ], sseph.p[k, 3j  ], ord)
                    # Distance between two bodies squared ||\mathbf{r}_i - \mathbf{r}_j||^2
                    TS.sqr!(X2, X, aux0, ord)
                    TS.sqr!(Y2, Y, aux0, ord)
                    TS.sqr!(Z2, Z, aux0, ord)
                    TS.add!(X2pY2, X2, Y2, ord)
                    TS.add!(r2, X2pY2, Z2, ord)
                    # Distance between two bodies ||\mathbf{r}_i - \mathbf{r}_j||
                    TS.sqrt!(r, r2, aux, ord)
                    # Newtonian potentials
                    TS.div!(GMi_r, GM[i], r, ord)
                    TS.add!(poteph[k, j], poteph[k, j], GMi_r, ord)
                    TS.div!(GMj_r, GM[j], r, ord)
                    TS.add!(poteph[k, i], poteph[k, i], GMj_r, ord)
                end
            end
        end
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
