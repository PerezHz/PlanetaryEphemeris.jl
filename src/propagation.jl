
@doc raw"""
    selecteph2jld(sseph::TaylorInterpolant, bodyind::AbstractVector{Int}, tspan::Number, N::Int)

Saves the ephemeris, contained in `sseph`, of the bodies with indexes `bodyind`, in a jld file
named as follows

    "sseph" * number of asteroids in sseph * "ast" * number of asteroids to be saved in file 
    * "p" (forward integration) or "m" (backward integration) * "y_et.jld"

# Arguments

- `sseph::TaylorInterpolant`: Ephemeris of all the bodies. 
- `bodyind::AbstractVector{Int}`: Indexes of the bodies to be saved.
- `tspan::Number`: Time span of the integration (positive -> forward integration / negative -> backward integration).
- `N::Int`: Total number of bodies.
"""
function selecteph2jld(sseph::TaylorInterpolant, bodyind::AbstractVector{Int}, tspan::Number, N::Int)
    nast = N - 11                      # Number of asteroids in sseph
    indvec = nbodyind(N, bodyind)      # Indexes of the positions and velocities of the bodies to be saved
    nastout = length(bodyind) - 11     # Number of asteroids to be saved
    @assert nastout <= nast
    sgn_yrs = signbit(tspan) ? "m" : "p"     # Prefix to distinguish between forward (p) / backward (m) integration
    nyrs_int = Int(abs(tspan))               # Number of years 
    
    # Write output to jld file

    # Name of the file
    ss16ast_fname = "sseph$(lpad(nast,3,'0'))ast$(lpad(nastout,3,'0'))_"*sgn_yrs*"$(nyrs_int)y_et.jld"
    # TaylorInterpolant with only the information of the bodies to be saved 
    # + Lunar orientation + TT-TDB
    ss16ast_eph = TaylorInterpolant(sseph.t0, sseph.t, sseph.x[:, union(indvec, 6N+1:6N+13)])
    # Open file 
    jldopen(ss16ast_fname, "w") do file
        addrequire(file, TaylorSeries)          # Require TaylorSeries
        addrequire(file, PlanetaryEphemeris)    # Require PlanetaryEphemeris
        write(file, "ss16ast_eph", ss16ast_eph) # Write the ephemeris to file 
    end
    # Check that written output is equal to original variable ss16ast_eph
    recovered_sol_i = load(ss16ast_fname, "ss16ast_eph")
    @show recovered_sol_i == ss16ast_eph
    return nothing
end

@doc raw"""
    propagate(maxsteps::Int, jd0::T, tspan::T; output::Bool=true, dense::Bool=false, 
              ephfile::String="sseph.jld", dynamics::Function=NBP_pN_A_J23E_J23M_J2S!, 
              nast::Int=343, quadmath::Bool=false, ss16ast::Bool=true, bodyind::AbstractVector{Int}=1:(11+nast),
              order::Int=order, abstol::T=abstol, parse_eqs::Bool=true) where {T<:Real}

Integrates the Solar System via the Taylor method. 

# Arguments 

- `maxsteps::Int`: maximum number of steps for the integration.
- `jd0::T`: initial Julian date.
- `tspan::T`: time span of the integration (in Julian days). 
- `output::Bool`: whether to write the output to a file (`true`) or not.
- `dense::Bool`: whether to save the Taylor polynomials at each step (`true`) or not.
- `ephfile::String`: name of the file where to save the solution if `output` is `true` but one or both of `dense` and `ss16ast` is `false`.
- `dynamics::Function`: dynamical model function.
- `nast::Int`: number of asteroids to be considered in the integration.
- `quadmath::Bool`: whether to use quadruple precision (`true`) or not. 
- `ss16ast::Bool`: wheter to save the solution using `selecteph2jld` (`true`) or not.
- `bodyind::AbstractVector{Int}`: indexes of the bodies to be saved.
- `order::Int=order`: order of the Taylor expansions to be used in the integration. 
- `abstol::T`: absolute tolerance.
- `parse_eqs::Bool`: whether to use the specialized method of `jetcoeffs!` (`true`) created with `@taylorize` or not.
"""
function propagate(maxsteps::Int, jd0::T, tspan::T; output::Bool=true, dense::Bool=false, 
                   ephfile::String="sseph.jld", dynamics::Function=NBP_pN_A_J23E_J23M_J2S!,
                   nast::Int=343, quadmath::Bool=false, ss16ast::Bool=true, 
                   bodyind::AbstractVector{Int}=1:(11+nast), order::Int=order,
                   abstol::T=abstol, parse_eqs::Bool=true) where {T<:Real}

    # Total number of bodies
    N = 11+nast
    # Get initial conditions (6N translational + 6 lunar mantle physical librations + 6 lunar core + TT-TDB)
    _q0 = initialcond(N, jd0) # <--- length(_q0) == 6N+13
    # Set initial time equal to zero (improves accuracy in data reductions)
    _t0 = zero(jd0)
    # Final time (julian days)
    @show _tmax = zero(_t0)+tspan*yr

    if quadmath
        # Use quadruple precision
        q0 = Float128.( _q0 )
        t0 = Float128(_t0)
        tmax = Float128(_tmax)
        _abstol = Float128(abstol)
        _jd0 = Float128(jd0)
    else
        q0 = _q0
        t0 = _t0
        tmax = _tmax
        _abstol = abstol
        _jd0 = jd0
    end

    # N: Total number of bodies
    # jd0: Initial Julian date
    params = (N, _jd0)

    # Do integration
    if dense
        # @time sol_ = taylorinteg(dynamics, q0, t0, tmax, order, _abstol, params, maxsteps=maxsteps, dense=dense)
        @time sol_ = taylorinteg_threads(dynamics, q0, t0, tmax, order, _abstol, params, maxsteps=maxsteps, dense=dense, parse_eqs=parse_eqs)
        # Parameters for TaylorInterpolant
        if quadmath  # with quadruple precision
            # Initial time (seconds)
            et0 = (jd0-J2000)*daysec
            # Vector of times (seconds)
            etv = Float64.( sol_.t[:]*daysec )
            # Vector of Taylor polynomials
            sseph_x_et = map( x->x(Taylor1(order)/daysec), map(x->Taylor1(Float64.(x.coeffs)), sol_.x[:,:]) )
        else
            # Initial time (seconds)
            et0 = (jd0-J2000)*daysec
            # Vector of times (seconds)
            etv = sol_.t[:]*daysec
            # Vector of Taylor polynomials
            sseph_x_et = map(x->x(Taylor1(order)/daysec), sol_.x[:,:])
        end
        # Save ephemeris in TaylorInterpolant object
        sseph = TaylorInterpolant(et0, etv, sseph_x_et)
        sol = (sseph=sseph,)
    else
        # @time sol_ = taylorinteg(dynamics, q0, t0, tmax, order, _abstol, params, maxsteps=maxsteps, dense=dense)
        @time sol_ = taylorinteg_threads(dynamics, q0, t0, tmax, order, _abstol, params, maxsteps=maxsteps, dense=dense, parse_eqs=parse_eqs)
        sol = (t=sol_[1][:], x=sol_[2][:,:])
    end

    # Write solution to .jld files
    if output
        if dense && ss16ast
            selecteph2jld(sseph, bodyind, tspan, N)
        else
            println("Saving solution to file: $ephfile")
            # Open file
            jldopen(ephfile, "w") do file
                addrequire(file, TaylorSeries)        # Require TaylorSeries
                addrequire(file, PlanetaryEphemeris)  # Require PlanetaryEphemeris
                # Write variables to jld file
                for ind in eachindex(sol)
                    varname = string(ind)
                    println("Saving variable: ", varname)
                    write(file, varname, sol[ind])
                end
            end
            # Check that recovered variables are equal to original variables
            for ind in eachindex(sol)
                varname = string(ind)
                # Read varname from jld file and assign recovered variable to recovered_sol_i
                recovered_sol_i = load(ephfile, varname)
                # Check that recovered variable is equal to original variable
                @show recovered_sol_i == sol[ind]
            end
        end
        println("Saved solution")
        return nothing
    else
        return sol
    end
end
