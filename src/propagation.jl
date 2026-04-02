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
    loadeph(ss16asteph::TaylorInterpolant, μ::Vector{<:Real})

Taking Solar System ephemeris `ss16asteph` and their gravitational parameters `μ` as input, returns for all bodies the point-mass Newtonian acceleration and the Newtonian N body potential.

# Arguments

- `ss16asteph`: Solar System ephemeris.
- `μ::Vector{<:Real}`: vector of mass parameters.
"""
function loadeph(ss16asteph::TaylorInterpolant, μ::Vector{<:Real})

    # Compute point-mass Newtonian accelerations from ephemeris
    # accelerations of all bodies are needed to compute the post-Newtonian acceleration of e.g. Solar System minor bodies
    # Number of bodies that contibute to the asteroid's acceleration
    Nm1 = numberofbodies(ss16asteph)
    # Initialize a TaylorInterpolant for the point-mass Newtonian accelerations
    acc_eph = TaylorInterpolant(ss16asteph.t0, ss16asteph.t, Matrix{eltype(ss16asteph.x)}(undef, length(ss16asteph.t)-1, 3Nm1))
    # Initialize a TaylorInterpolant for the newtonian N body potential
    pot_eph = TaylorInterpolant(ss16asteph.t0, ss16asteph.t, Matrix{eltype(ss16asteph.x)}(undef, length(ss16asteph.t)-1, Nm1))
    # Fill TaylorInterpolant.x with zero polynomials
    fill!(acc_eph.x, zero(ss16asteph.x[1]))
    fill!(pot_eph.x, zero(ss16asteph.x[1]))

    # Iterator over all bodies except asteroid
    for j in 1:Nm1
        for i in 1:Nm1
            if i == j
                #
            else
                # Difference between two positions (\mathbf{r}_i - \mathbf{r}_j)
                X_ij = ss16asteph.x[:, 3i-2] .- ss16asteph.x[:, 3j-2]  # X-axis component
                Y_ij = ss16asteph.x[:, 3i-1] .- ss16asteph.x[:, 3j-1]  # Y-axis component
                Z_ij = ss16asteph.x[:, 3i  ] .- ss16asteph.x[:, 3j  ]  # Z-axis component
                # Distance between two bodies squared ||\mathbf{r}_i - \mathbf{r}_j||^2
                r_p2_ij = ( (X_ij.^2) .+ (Y_ij.^2) ) .+ (Z_ij.^2)
                # Distance between two bodies ||\mathbf{r}_i - \mathbf{r}_j||
                r_ij = sqrt.(r_p2_ij)
                # Newtonian potential
                pot_eph.x[:, j] .+= (μ[i]./r_ij)
            end
        end

        # Fill acelerations by differentiating velocities
        acc_eph.x[:, 3j-2] .= ordpres_differentiate.(ss16asteph.x[:, 3(Nm1+j)-2])  # X-axis component
        acc_eph.x[:, 3j-1] .= ordpres_differentiate.(ss16asteph.x[:, 3(Nm1+j)-1])  # Y-axis component
        acc_eph.x[:, 3j  ] .= ordpres_differentiate.(ss16asteph.x[:, 3(Nm1+j)  ])  # Z-axis component
    end

    return acc_eph, pot_eph
end

"""
    selecteph2jld2(sseph::TaylorInterpolant, bodyind, nyears)

Save a subset of `sseph` containing only the ephemeris of the
`bodyind`-th bodies in a `.jld2` file named as follows

    sseph{N}ast{n}_{p/m}{nyears}y_et.jld2

where `N` is the number of asteroids in sseph, `n` is the number
of asteroids to be saved in the file, `p/m` indicates a forward
or backward integration and `nyears` is the number of years.
"""
function selecteph2jld2(sseph::TaylorInterpolant, bodyind::AbstractVector{Int},
                        nyears::Number)
    # Total number of bodies
    N = numberofbodies(sseph)
    # Number of asteroids in sseph
    nast = N - 11
    # Number of asteroids to be saved
    nastout = length(bodyind) - 11
    # Check nastout <= nast
    @assert nastout <= nast "Cannot save $nastout asteroids from ephemeris \
        with $nast asteroids"
    # Prefix to distinguish between forward (p) / backward (m) integration
    sgn_yrs = signbit(nyears) ? "m" : "p"
    # Number of years
    nyrs_int = floor(Int, abs(nyears))
    # Name of the file
    ss16ast_fname = "sseph$(lpad(nast,3,'0'))ast$(lpad(nastout,3,'0'))_" *
        sgn_yrs * "$(nyrs_int)y_et.jld2"
    # Select bodies to be saved + Lunar orientation + TT-TDB
    t0, tf = timespan(sseph)
    ss16ast_eph = selecteph(sseph, bodyind, t0, tf, euler = true, ttmtdb = true)
    # Open file
    println("Saving solution to file: ", ss16ast_fname)
    jldsave(ss16ast_fname; ss16ast_eph)
    # Check that written output is equal to original variable ss16ast_eph
    recovered_sol_i = JLD2.load(ss16ast_fname, "ss16ast_eph")
    if recovered_sol_i == ss16ast_eph
        println("Solution saved correctly")
    else
        println("Saved and recovered solution are not equal")
    end

    return ss16ast_fname
end

@doc raw"""
    save2jld2andcheck(outfilename::String, sol)

Save `sol` in `outfilename` (.jld2) and check that recovered solution equals `sol`.
"""
function save2jld2andcheck(outfilename::String, sol)

    println("Saving solution to file: ", outfilename)

    # Open file
    JLD2.jldopen(outfilename, "w") do file
        # Loop over solution variables
        for ind in eachindex(sol)
            # Name of the variable
            varname = string(ind)
            println("Saving variable: ", varname)
            # Write the variable
            write(file, varname, sol[ind])
        end
    end

    # Check that saved solution is equal to the original
    println("Checking that all variables were saved correctly...")

    # Loop over solution variables
    for ind in eachindex(sol)
        # Name of the variable
        varname = string(ind)
        # Read varname from files and assign recovered variable to recovered_sol_i
        recovered_sol_i = JLD2.load(outfilename, varname)
        # Check that varname was recovered succesfully
        if recovered_sol_i == sol[ind]
            println("Variable ", varname, " saved correctly" )
        else
            println("Recovered variable ", varname, " is not equal to the original" )
        end
    end

    println("Saved solution")

    return nothing
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
                   maxsteps::Int = 500, order::Int = order,
                   abstol::T = abstol, parse_eqs::Bool = true) where {D, T, P}
    # Unpack
    @unpack dynamics, tspan, initcond, params = PE
    # Integration
    sol = taylorinteg(dynamics, initcond, zero(T), tspan[2] - tspan[1], order,
                      abstol, params; maxsteps, parse_eqs)
    # Convert from TaylorSolution to TaylorInterpolant
    return TaylorInterpolant{T, T, 2}(tspan[1] - J2000, sol.t, sol.p)
end
