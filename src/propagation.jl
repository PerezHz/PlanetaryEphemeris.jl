function propagate(maxsteps::Int, jd0::T, tspan::T, eulangfile::String;
        output::Bool=true, dense::Bool=false, ephfile::String="sseph.jld",
        dynamics::Function=NBP_pN_A_J23E_J23M_J2S!, nast::Int=343,
        quadmath::Bool=false, ss16ast::Bool=true) where {T<:Real}

    # total number of bodies
    N = 11+nast
    # get initial conditions
    if quadmath
        # use quadruple precision
        _q0 = Float128.( initialcond(11+nast) )
        _t0 = Float128(zero(jd0))
        @show _tmax = zero(_t0)+tspan*yr #final time of integration
        _abstol = Float128(abstol)
        # load DE430 lunar Euler angles Taylor ephemeris
        __eulang_de430 = load(eulangfile, "eulang_de430")
        _eulang_de430 = TaylorInterpolant(Float128(__eulang_de430.t0), Float128.(__eulang_de430.t), map(x->Taylor1(Float128.(x.coeffs)), __eulang_de430.x))
    else
        _q0 = initialcond(11+nast)
        _t0 = zero(jd0)
        @show _tmax = zero(_t0)+tspan*yr #final time of integration
        _abstol = abstol
        # load DE430 lunar Euler angles Taylor ephemeris
        _eulang_de430 = load(eulangfile, "eulang_de430")
    end
    # auxiliary variable
    S = eltype(_q0)

    params = (N, S, _eulang_de430, jd0)

    # do integration
    if dense
        # @time sol_ = taylorinteg(dynamics, _q0, _t0, _tmax, order, _abstol, params, maxsteps=maxsteps, dense=dense)
        @time sol_ = taylorinteg_threads(dynamics, _q0, _t0, _tmax, order, _abstol, params, maxsteps=maxsteps, dense=dense)
        # transform Julian dates to TDB seconds since initial Julian date
        if quadmath
            et0 = (jd0-J2000)*daysec
            etv = Float64.( sol_.t[:]*daysec )
            sseph_x_et = map( x->x(Taylor1(order)/daysec), map(x->Taylor1(Float64.(x.coeffs)), sol_.x[:,:]) )
        else
            et0 = (jd0-J2000)*daysec
            etv = sol_.t[:]*daysec
            sseph_x_et = map(x->x(Taylor1(order)/daysec), sol_.x[:,:])
        end
        sseph = TaylorInterpolant(et0, etv, sseph_x_et)
        sol = (sseph=sseph,)
    else
        # @time sol_ = taylorinteg(dynamics, _q0, _t0, _tmax, order, _abstol, params, maxsteps=maxsteps, dense=dense)
        @time sol_ = taylorinteg_threads(dynamics, _q0, _t0, _tmax, order, _abstol, params, maxsteps=maxsteps, dense=dense)
        sol = (t=sol_[1][:], x=sol_[2][:,:])
    end

    #write solution to .jld files
    if output
        if ss16ast
            i1 = 1
            i2 = 27 # 11 major bodies + 16 asteroids
            i2 > N && (i2 = N)
            indvec = union(3i1-2:3i2, 3*(N+i1)-2:3*(N+i2))
            sgn_yrs = sign(tspan) == 1.0 ? "p" : "m"
            nyrs_int = Int(abs(tspan))
            # write output to jld file
            ss16ast_fname = "ss16ast343_eph_"*sgn_yrs*"$(nyrs_int)y_et.jld"
            ss16ast_eph = TaylorInterpolant(sseph.t0, sseph.t, sseph.x[:, indvec])
            jldopen(ss16ast_fname, "w") do file
                write(file, "ss16ast_eph", ss16ast_eph)
            end
            #check that written output is equal to original variable `ss16ast_eph`
            recovered_sol_i = load(ss16ast_fname, "ss16ast_eph")
            @show recovered_sol_i == ss16ast_eph
        else
            println("Saving solution to file: $ephfile")
            jldopen(ephfile, "w") do file
                # write variables to jld file
                for ind in eachindex(sol)
                    varname = string(ind)
                    println("Saving variable: ", varname)
                    write(file, varname, sol[ind])
                end
            end
            #check that recovered variables are equal to original variables
            for ind in eachindex(sol)
                varname = string(ind)
                #read varname from jld file and assign recovered variable to recovered_sol_i
                recovered_sol_i = load(ephfile, varname)
                #check that recovered variable is equal to original variable
                @show recovered_sol_i == sol[ind]
            end
        end
        println("Saved solution")
        return nothing
    else
        return sol
    end
end
