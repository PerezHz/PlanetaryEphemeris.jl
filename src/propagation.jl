@doc raw"""
    selecteph2jld2(sseph::TaylorInterpolant, bodyind::AbstractVector{Int}, tspan::Number, N::Int)

Save the ephemeris, contained in `sseph`, of the bodies with indexes `bodyind`, in a `jld2` file
named as follows

    "sseph" * number of asteroids in sseph * "ast" * number of asteroids to be saved in file 
    * "p" (forward integration) or "m" (backward integration) * "y_et.jld2"

# Arguments

- `sseph::TaylorInterpolant`: ephemeris of all the bodies. 
- `bodyind::AbstractVector{Int}`: indexes of the bodies to be saved.
- `tspan::Number`: time span of the integration (positive -> forward integration / negative -> backward integration).
- `N::Int`: total number of bodies.
"""
function selecteph2jld2(sseph::TaylorInterpolant, bodyind::AbstractVector{Int}, tspan::Number, N::Int)
    # Number of asteroids in sseph
    nast = N - 11                 
    # Indexes of the positions and velocities of the bodies to be saved
    indvec = nbodyind(N, bodyind)      
    # Number of asteroids to be saved
    nastout = length(bodyind) - 11     
    @assert nastout <= nast "Cannot save $nastout asteroids from ephemeris with $nast asteroids"
    # Prefix to distinguish between forward (p) / backward (m) integration
    sgn_yrs = signbit(tspan) ? "m" : "p"     
    # Number of years
    nyrs_int = Int(abs(tspan))                
    
    # Write output to jld2 file

    # Name of the file
    ss16ast_fname = "sseph$(lpad(nast,3,'0'))ast$(lpad(nastout,3,'0'))_"*sgn_yrs*"$(nyrs_int)y_et.jld2"
    # TaylorInterpolant with only the information of the bodies to be saved 
    # + Lunar orientation + TT-TDB
    ss16ast_eph = TaylorInterpolant(sseph.t0, sseph.t, sseph.x[:, union(indvec, 6N+1:6N+13)])
    # Open file 
    JLD2.jldopen(ss16ast_fname, "w") do file
        # Write the ephemeris to file 
        write(file, "ss16ast_eph", ss16ast_eph) 
    end
    # Check that written output is equal to original variable ss16ast_eph
    recovered_sol_i = load(ss16ast_fname, "ss16ast_eph")
    if recovered_sol_i == ss16ast_eph
        println("ss16ast_eph variable saved correctly")
    end 

    println("Saved solution")

    return nothing
end

@doc raw"""
    save2jld2andcheck(outfilename, sol)

Save `sol` in `outfilename` (.jld2). 
"""
function save2jld2andcheck(outfilename, sol)
    println("Saving solution to file: $outfilename")
    # Open file 
    JLD2.jldopen(outfilename, "w") do file
        # Loop over solution variables
        for ind in eachindex(sol)
            # Name of the variable 
            varname = string(ind)
            println("Saving variable: ", varname)
            # Write the varaible 
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
        end 
    end
    println("Saved solution")
    return outfilename
end

@doc raw"""
    propagate_params(jd0::T, nast::Int = 343) where {T <: Real}

Return initial time, initial conditions and the parameters necessary for the integration of the Solar System. 
"""
function propagate_params(jd0::T, nast::Int = 343) where {T <: Real}

    # Total number of bodies
    N = 11 + nast

    # Get initial conditions (6N translational + 6 lunar mantle physical librations + 6 lunar core + TT-TDB)
    q0 = initialcond(N, jd0) # length(_q0) == 6N + 13

    # Set initial time equal to zero (improves accuracy in data reductions)
    t0 = zero(T)

    # N: Total number of bodies
    # jd0: Initial Julian date
    params = (N, jd0)

    return t0, q0, params 

end 

@doc raw"""
    propagate(maxsteps::Int, jd0::T, tspan::T, ::Val{false/true}; output::Bool = true, ephfile::String = "sseph.jld2", 
              dynamics::Function = NBP_pN_A_J23E_J23M_J2S!, nast::Int = 343, ss16ast::Bool = true, 
              bodyind::AbstractVector{Int} = 1:(11+nast), order::Int = order, abstol::T = abstol,
              parse_eqs::Bool = true) where {T <: Real}

Integrate the Solar System via the Taylor method. 

# Arguments 

- `maxsteps::Int`: maximum number of steps for the integration.
- `jd0::T`: initial Julian date.
- `tspan::T`: time span of the integration (in Julian days). 
- `::Val{false/true}`: whether to save the Taylor polynomials at each step (`true`) or not (`false`).
- `output::Bool`: whether to write the output to a file (`true`) or not.
- `ephfile::String`: name of the file where to save the solution if `ss16ast` is `false`.
- `dynamics::Function`: dynamical model function.
- `nast::Int`: number of asteroids to be considered in the integration.
- `ss16ast::Bool`: whether to save the solution using `selecteph2jld2` (`true`) or not.
- `bodyind::AbstractVector{Int}`: indexes of the bodies to be saved.
- `order::Int=order`: order of the Taylor expansions to be used in the integration. 
- `abstol::T`: absolute tolerance.
- `parse_eqs::Bool`: whether to use the specialized method of `jetcoeffs!` (`true`) created with `@taylorize` or not.
""" propagate

for V_dense in (true, false)
    @eval begin

        function propagate(maxsteps::Int, jd0::T, tspan::T, ::Val{$V_dense}; output::Bool = true, ephfile::String = "sseph.jld2", 
                           dynamics::Function = NBP_pN_A_J23E_J23M_J2S!, nast::Int = 343, ss16ast::Bool = true, 
                           bodyind::AbstractVector{Int} = 1:(11+nast), order::Int = order, abstol::T = abstol,
                           parse_eqs::Bool = true) where {T <: Real}

            t0, q0, params = propagate_params(jd0, nast)
        
            # Final time (julian days)
            println("Initial time of integration: ", julian2datetime(jd0))
            # Final time of integration (days)
            tmax = t0 + tspan*yr
            println("Final time of integration: ", julian2datetime(jd0 + tmax))
        
            # Integration 
            @time sol_ = taylorinteg_threads(dynamics, q0, t0, tmax, order, abstol, Val($V_dense), params, maxsteps = maxsteps,
                                            parse_eqs = parse_eqs)

            if $V_dense 

                # Parameters for TaylorInterpolant

                # Initial time (seconds)
                et0 = T( (jd0 - J2000) * daysec )

                # Vector of times (seconds)
                etv = T.( sol_.t[:]*daysec )
                
                # Vector of Taylor polynomials
                sseph_x_et = map( x -> x(Taylor1(order)/daysec), map(x -> Taylor1(T.(x.coeffs)), sol_.x[:, :]) )

                # Element type of initial conditions 
                U = eltype(q0)

                # Save ephemeris in TaylorInterpolant object
                sseph = TaylorInterpolant{T, U, 2}(et0, etv, sseph_x_et)

                sol = (sseph = sseph,)

            else
                
                sol = (t = sol_[1][:], x = sol_[2][:, :])

            end 

            # Write solution to .jld2 files
            if output
                if $V_dense 
                    if ss16ast
                        selecteph2jld2(sseph, bodyind, tspan, params[1])
                    else
                        _ = save2jldandcheck(ephfile, sol)
                    end
                else
                    _ = save2jldandcheck(ephfile, sol)
                end
            end

            return sol 

        end 

        function propagate(maxsteps::Int, jd0::T1, tspan::T2, ::Val{$V_dense}; output::Bool = true, ephfile::String = "sseph.jld", 
                           dynamics::Function = NBP_pN_A_J23E_J23M_J2S!, nast::Int = 343, ss16ast::Bool = true, 
                           bodyind::AbstractVector{Int} = 1:(11+nast), order::Int = order, abstol::T3 = abstol,
                           parse_eqs::Bool = true) where {T1, T2, T3 <: Real}
            
            _jd0, _tspan, _abstol = promote(jd0, tspan, abstol)                           

            return propagate(maxsteps, _jd0, _tspan, Val($V_dense); output = output, ephfile = ephfile, dynamics = dynamics, 
                             nast = nast, ss16ast = ss16ast, bodyind = bodyind, order = order, abstol = abstol, parse_eqs = parse_eqs)
        end 

    end
end 