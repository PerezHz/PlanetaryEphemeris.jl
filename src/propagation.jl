# This file is part of the PlanetaryEphemeris.jl package; MIT licensed

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
        acc_eph.x[:, 3j-2] .= NEOs.ordpres_differentiate.(ss16asteph.x[:, 3(Nm1+j)-2])  # X-axis component
        acc_eph.x[:, 3j-1] .= NEOs.ordpres_differentiate.(ss16asteph.x[:, 3(Nm1+j)-1])  # Y-axis component
        acc_eph.x[:, 3j  ] .= NEOs.ordpres_differentiate.(ss16asteph.x[:, 3(Nm1+j)  ])  # Z-axis component
    end

    return acc_eph, pot_eph
end

@doc raw"""
    selecteph2jld2(sseph::TaylorInterpolant, bodyind::T, tspan::S) where {T <: AbstractVector{Int}, S <: Number}

Save the ephemeris, contained in `sseph`, of the bodies with indices `bodyind`, in a `.jld2` file named as follows

    "sseph" * number of asteroids in sseph * "ast" * number of asteroids to be saved in file * "_"
    * "p" / "m" (forward / backward integration) * number of years in sseph *  "y_et.jld2"

# Arguments

- `sseph::TaylorInterpolant`: ephemeris of all the bodies.
- `bodyind::T`: indices of the bodies to be saved.
- `tspan::S`: time span of the integration (positive -> forward integration / negative -> backward integration).
"""
function selecteph2jld2(sseph::TaylorInterpolant, bodyind::T, tspan::S) where {T <: AbstractVector{Int}, S <: Number}

    # Total number of bodies
    N = numberofbodies(sseph)
    # Number of asteroids in sseph
    nast = N - 11
    # Number of asteroids to be saved
    nastout = length(bodyind) - 11
    # Check nastout <= nast
    @assert nastout <= nast "Cannot save $nastout asteroids from ephemeris with $nast asteroids"
    # Prefix to distinguish between forward (p) / backward (m) integration
    sgn_yrs = signbit(tspan) ? "m" : "p"
    # Number of years
    nyrs_int = floor(Int, abs(tspan))

    # Write output to .jld2 file

    # Name of the file
    ss16ast_fname = "sseph$(lpad(nast,3,'0'))ast$(lpad(nastout,3,'0'))_" * sgn_yrs * "$(nyrs_int)y_et.jld2"

    # TaylorInterpolant with only the information of the bodies to be saved + Lunar orientation + TT-TDB
    ss16ast_eph = selecteph(sseph, bodyind, euler = true, ttmtdb = true)

    println("Saving solution to file: ", ss16ast_fname)

    # Open file
    JLD2.jldopen(ss16ast_fname, "w") do file
        # Write the ephemeris to file
        write(file, "ss16ast_eph", ss16ast_eph)
    end
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

@doc raw"""
    propagate(maxsteps::Int, jd0::T, tspan::T; dynamics::Function = NBP_pN_A_J23E_J23M_J2S_threads!,
              nast::Int = 343, order::Int = order, abstol::T = abstol, parse_eqs::Bool = true) where {T <: Real}

Integrate the Solar System via the Taylor method.

# Arguments

- `maxsteps::Int`: maximum number of steps for the integration.
- `jd0::T`: initial Julian date.
- `tspan::T`: time span of the integration (in Julian days).
- `dynamics::Function`: dynamical model function.
- `nast::Int`: number of asteroids to be considered in the integration.
- `order::Int`: order of the Taylor expansions to be used in the integration.
- `abstol::T`: absolute tolerance.
- `parse_eqs::Bool`: whether to use the specialized method of `jetcoeffs!` (`true`) created with `@taylorize` or not.
""" propagate

function propagate(maxsteps::Int, jd0::T, tspan::T; dynamics::Function = NBP_pN_A_J23E_J23M_J2S_threads!,
                    nast::Int = 343, order::Int = order, abstol::T = abstol, parse_eqs::Bool = true) where {T <: Real}

    # Total number of bodies (Sun + 8 planets + Moon + Pluto + Asteroid)
    N = 11 + nast

    # Get 6N + 13 initial conditions (3N positions + 3N velocities + 6 lunar mantle angles + 6 lunar core angles + TT-TDB)
    q0 = initialcond(N, jd0)

    # Set initial time equal to zero (improves accuracy in data reductions)
    t0 = zero(T)

    # Parameters for dynamical function
    params = (N, jd0)

    # Final time of integration (days)
    tmax = t0 + tspan*yr

    # Integration
    sol = @time taylorinteg(dynamics, q0, t0, tmax, order, abstol, params;
                            maxsteps, parse_eqs)

    return TaylorInterpolant{T, T, 2}(jd0 - J2000, sol.t, sol.p)
end

function propagate(maxsteps::Int, jd0::T1, tspan::T2; dynamics::Function = NBP_pN_A_J23E_J23M_J2S_threads!,
                    nast::Int = 343, order::Int = order, abstol::T3 = abstol, parse_eqs::Bool = true) where {T1, T2, T3 <: Real}

    _jd0, _tspan, _abstol = promote(jd0, tspan, abstol)

    return propagate(maxsteps, _jd0, _tspan; dynamics = dynamics, nast = nast, order = order,
                        abstol = _abstol, parse_eqs = parse_eqs)

end
