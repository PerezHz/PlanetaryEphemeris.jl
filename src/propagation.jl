# This file is part of the PlanetaryEphemeris.jl package; MIT licensed

@doc raw"""
    numberofbodies(::TaylorSolution)

Return the number of bodies saved in a `TaylorSolution` produced by `PlanetaryEphemeris`.
"""
numberofbodies(L::Int) = (L - 13) ÷ 6
numberofbodies(v::AbstractVector{T}) where {T} = numberofbodies(length(v))
numberofbodies(m::AbstractMatrix{T}) where {T} = numberofbodies(size(m, 2))
numberofbodies(eph::TaylorSolution) = numberofbodies(size(eph.p, 2))

@doc raw"""
    loadeph(sseph::TaylorSolution, μ::Vector{T}) where {T <: Real}

Taking Solar System ephemeris `sseph` and their gravitational parameters `μ` as input,
returns for all bodies the point-mass Newtonian acceleration and the Newtonian N-body
potential.
"""
function loadeph(sseph::TaylorSolution, μ::Vector{T}) where {T <: Real}
    # Number of bodies and timesteps
    Nb = numberofbodies(sseph)
    Nt = length(sseph.t)
    # Initialize memory for the Newtonian point-mass accelerations and N-body potentials
    acc_eph = TaylorSolution(sseph.t, [zero(sseph.x[1]) for i in 1:Nt, j in 1:3Nb],
        [zero(sseph.p[1]) for i in 1:Nt-1, j in 1:3Nb])
    pot_eph = TaylorSolution(sseph.t, [zero(sseph.x[1]) for i in 1:Nt, j in 1:Nb],
        [zero(sseph.p[1]) for i in 1:Nt-1, j in 1:Nb])
    # Iterate over all bodies
    for j in 1:Nb
        for i in 1:Nb
            if i == j
                #
            else
                # Difference between two positions (\mathbf{r}_i - \mathbf{r}_j)
                X_ij = sseph.p[:, 3i-2] - sseph.p[:, 3j-2]
                Y_ij = sseph.p[:, 3i-1] - sseph.p[:, 3j-1]
                Z_ij = sseph.p[:, 3i  ] - sseph.p[:, 3j  ]
                # Distance between two bodies squared ||\mathbf{r}_i - \mathbf{r}_j||^2
                r_p2_ij = @. ( (X_ij^2) + (Y_ij^2) ) + (Z_ij^2)
                # Distance between two bodies ||\mathbf{r}_i - \mathbf{r}_j||
                r_ij = @. sqrt(r_p2_ij)
                # Newtonian potential
                @. pot_eph.p[:, j] += (μ[i]/r_ij)
            end
        end
        # Fill acelerations by differentiating velocities
        @. acc_eph.p[:, 3j-2] = ordpres_differentiate(sseph.p[:, 3(Nb+j)-2])
        @. acc_eph.p[:, 3j-1] = ordpres_differentiate(sseph.p[:, 3(Nb+j)-1])
        @. acc_eph.p[:, 3j  ] = ordpres_differentiate(sseph.p[:, 3(Nb+j)  ])
    end

    return acc_eph, pot_eph
end

@doc raw"""
    selecteph(eph::TaylorSolution, bodyind::Union{Int, AbstractVector{Int}},
        t0::T = eph.t[1], tf::T = eph.t[end]; kwargs...) where {T <: Real}

Return a subset of `eph` containing bodies `bodyind` in timerange `[t0, tf]`.

## Keyword arguments

- `euler::Bool`: whether to include lunar euler angles (default: `false`).
- `ttmtdb::Bool`: whether to include TT-TDB (default: `false`).
"""
function selecteph(eph::TaylorSolution, bodyind::Union{Int, AbstractVector{Int}},
    t0::T = eph.t[1], tf::T = eph.t[end]; euler::Bool = false,
    ttmtdb::Bool = false) where {T <: Real}
    # Check whether eph spans [t0, tf]
    tmin, tmax = minmax(eph.t[1], eph.t[end])
    @assert tmin ≤ t0 ≤ tf ≤ tmax "[t0, tf] beyond ephemeris timerange [$tmin, $tmax]"
    # Subset of eph.t spanning [t0, tf]
    if issorted(eph.t)
        j0 = searchsortedlast(eph.t, t0)
        jf = searchsortedfirst(eph.t, tf)
    else
        j0 = searchsortedlast(eph.t, tf, rev = true)
        jf = searchsortedfirst(eph.t, t0, rev = true)
    end
    t = eph.t[j0:jf]
    # Indices of bodies in bodyind (+ lunar euler angles + TT-TDB)
    N = numberofbodies(eph)
    @assert all(bodyind .< N) "bodyind is beyond ephemeris number of bodies ($N)"
    idxs = nbodyind(N, bodyind)
    if euler
        idxs = vcat(idxs, 6N+1:6N+12)
    end
    if ttmtdb
        idxs = vcat(idxs, 6N+13)
    end
    # Subsets of eph.x and eph.p containing bodies in bodyind and spanning [t0, tf]
    x = eph.x[j0:jf, idxs]
    p = eph.p[j0:jf-1, idxs]

    return TaylorSolution(t, x, p)
end

@doc raw"""
    selecteph2jld2(sseph::TaylorSolution, bodyind::AbstractVector{Int},
        tspan::Number)

Save the ephemeris, contained in `sseph`, of the bodies with indices `bodyind`,
in a `.jld2` file named as follows:

    "sseph" * number of asteroids in sseph * "ast" * number of asteroids to be
    saved in file * "_" * signbit(tspan) ? "m" : "p" * number of years in sseph *
    "y_et.jld2"
"""
function selecteph2jld2(sseph::TaylorSolution, bodyind::AbstractVector{Int},
    tspan::Number)
    # Total number of bodies
    N = numberofbodies(sseph)
    # Number of asteroids in sseph
    nast = N - 11
    # Number of asteroids to be saved
    nastout = length(bodyind) - 11
    # Check nastout <= nast
    @assert nastout <= nast "Cannot save $nastout asteroids from ephemeris with \
        $nast asteroids"
    # Prefix to distinguish between forward (p) / backward (m) integration
    sgn_yrs = signbit(tspan) ? "m" : "p"
    # Number of years
    nyrs_int = floor(Int, abs(tspan))
    # Name of the (.jld2) file
    ssfname = "sseph$(lpad(nast,3,'0'))ast$(lpad(nastout,3,'0'))_" * sgn_yrs *
        "$(nyrs_int)y_et.jld2"
    # TaylorSolution with only the information of the bodies to be saved +
    # Lunar orientation + TT-TDB
    sseph = selecteph(sseph, bodyind, euler = true, ttmtdb = true)

    println("Saving solution to file: ", ssfname)

    # Open file
    JLD2.jldopen(ssfname, "w") do file
        # Write the ephemeris to file
        write(file, "sseph", sseph)
    end
    # Check that written output is equal to original variable sseph
    recovered_sol = JLD2.load(ssfname, "sseph")
    if recovered_sol == sseph
        println("Solution saved correctly")
    else
        println("Saved and recovered solution are not equal")
    end

    return ssfname
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
    propagate(maxsteps::Int, jd0::T, tspan::T; kwargs...) where {T <: Real}

Integrate the Solar System via the Taylor method.

## Arguments

- `maxsteps::Int`: maximum number of steps for the integration.
- `jd0::T`: initial Julian date.
- `tspan::T`: time span of the integration (in Julian days).

## Keyword arguments

- `dynamics::Function`: dynamical model (default: `NBP_pN_A_J23E_J23M_J2S_threads!`).
- `nast::Int`: number of asteroid perturbers (default: `343`).
- `order::Int`: order of the Taylor expansions (default: `25`).
- `abstol::T`: absolute tolerance (default: `1.0E-20`).
- `parse_eqs::Bool`: whether to use the specialized method of `jetcoeffs!` or not
    (default: `true`).
"""
function propagate(maxsteps::Int, jd0::T, tspan::T;
    dynamics::Function = NBP_pN_A_J23E_J23M_J2S_threads!,
    nast::Int = 343, order::Int = order, abstol::T = abstol,
    parse_eqs::Bool = true) where {T <: Real}
    # Total number of bodies (Sun + 8 planets + Moon + Pluto + Asteroid)
    N = 11 + nast
    # Get 6N + 13 initial conditions (3N positions + 3N velocities +
    # 6 lunar mantle angles + 6 lunar core angles + TT-TDB)
    q0 = initialcond(N, jd0)
    # Set initial time equal to zero (improves accuracy in data reductions)
    t0 = zero(T)
    # Parameters for dynamical function
    params = (N, jd0)
    # Final time of integration (days)
    tmax = t0 + tspan*yr
    # Integration
    @time sol = taylorinteg(dynamics, q0, t0, tmax, order, abstol, params;
        maxsteps, parse_eqs)
    # Shift initial time
    @. sol.t += jd0 - J2000

    return sol
end

function propagate(maxsteps::Int, jd0::T1, tspan::T2;
    dynamics::Function = NBP_pN_A_J23E_J23M_J2S_threads!,
    nast::Int = 343, order::Int = order, abstol::T3 = abstol,
    parse_eqs::Bool = true) where {T1, T2, T3 <: Real}
    # Type promotion
    _jd0, _tspan, _abstol = promote(jd0, tspan, abstol)
    # Integration
    return propagate(maxsteps, _jd0, _tspan; dynamics, nast, order,
        abstol = _abstol, parse_eqs = parse_eqs)
end
