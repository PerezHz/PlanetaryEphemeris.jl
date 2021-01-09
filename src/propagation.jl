nbodyind(N::Int, i::Int) = union(3i-2:3i, 3*(N+i)-2:3*(N+i))

function nbodyind(N::Int, ivec::AbstractVector{Int})
    a = Int[]
    for i in ivec
        i > N && continue
        a = union(a, nbodyind(N,i))
    end
    return sort(a)
end

function selecteph2jld(sseph::TaylorInterpolant, bodyind::AbstractVector{Int}, tspan::Number)
    N = size(sseph.x)[2] ÷ 6 # total number of bodies (Sun+planets+Moon+Pluto+asts)
    nast = N - 11 # number of asteroids in sseph
    indvec = nbodyind(N, bodyind)
    nastout = length(bodyind) - 11 # number of asteroids to be saved in file
    @assert nastout <= nast
    sgn_yrs = sign(tspan) == 1.0 ? "p" : "m"
    nyrs_int = Int(abs(tspan))
    # write output to jld file
    ss16ast_fname = "sseph$(lpad(nast,3,'0'))ast$(lpad(nastout,3,'0'))_"*sgn_yrs*"$(nyrs_int)y_et.jld"
    ss16ast_eph = TaylorInterpolant(sseph.t0, sseph.t, sseph.x[:, indvec])
    jldopen(ss16ast_fname, "w") do file
        addrequire(file, TaylorSeries)
        addrequire(file, PlanetaryEphemeris)
        write(file, "ss16ast_eph", ss16ast_eph)
    end
    #check that written output is equal to original variable `ss16ast_eph`
    recovered_sol_i = load(ss16ast_fname, "ss16ast_eph")
    @show recovered_sol_i == ss16ast_eph
    return nothing
end

function propagate(maxsteps::Int, jd0::T, tspan::T, eulangfile::String;
        output::Bool=true, dense::Bool=false, ephfile::String="sseph.jld",
        dynamics::Function=NBP_pN_A_J23E_J23M_J2S!, nast::Int=343,
        quadmath::Bool=false, ss16ast::Bool=true,
        bodyind::AbstractVector{Int}=1:(11+nast)) where {T<:Real}

    # total number of bodies
    N = 11+nast
    # get initial conditions
    _q0 = initialcond(N)
    # initial time (Julian date)
    _t0 = zero(jd0)
    # final time (Julian date)
    @show _tmax = zero(_t0)+tspan*yr
    # load DE430 lunar Euler angles ephemeris (TaylorInterpolant)
    _eulang_de430 = load(eulangfile, "eulang_de430")

    if quadmath
        # use quadruple precision
        q0 = Float128.( _q0 )
        t0 = Float128(_t0)
        tmax = Float128(_tmax)
        _abstol = Float128(abstol)
        eulang_de430 = TaylorInterpolant(Float128(_eulang_de430.t0), Float128.(_eulang_de430.t), map(x->Taylor1(Float128.(x.coeffs)), _eulang_de430.x))
    else
        q0 = _q0
        t0 = _t0
        tmax = _tmax
        _abstol = abstol
        eulang_de430 = _eulang_de430
    end

    params = (N, eulang_de430, jd0)

    # do integration
    if dense
        # @time sol_ = taylorinteg(dynamics, q0, t0, tmax, order, _abstol, params, maxsteps=maxsteps, dense=dense)
        @time sol_ = taylorinteg_threads(dynamics, q0, t0, tmax, order, _abstol, params, maxsteps=maxsteps, dense=dense)
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
        # @time sol_ = taylorinteg(dynamics, q0, t0, tmax, order, _abstol, params, maxsteps=maxsteps, dense=dense)
        @time sol_ = taylorinteg_threads(dynamics, q0, t0, tmax, order, _abstol, params, maxsteps=maxsteps, dense=dense)
        sol = (t=sol_[1][:], x=sol_[2][:,:])
    end

    #write solution to .jld files
    if output
        if ss16ast
            selecteph2jld(sseph, bodyind, tspan)
        else
            println("Saving solution to file: $ephfile")
            jldopen(ephfile, "w") do file
                addrequire(file, TaylorSeries)
                addrequire(file, PlanetaryEphemeris)
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
