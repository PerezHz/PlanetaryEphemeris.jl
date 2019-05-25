function propagate(maxsteps::Int, t0::T, tspan::T;
        output::Bool=true, dense::Bool=false) where {T<:Real}

    # get initial conditions
    q0 = initialcond(length(μ))

    @show tmax = t0+tspan*yr #final time of integration

    # do integration
    if dense
        @time interp = taylorinteg(NBP_pN_A_J23E_J23M_J2S!, q0, t0, tmax, order, abstol, maxsteps=maxsteps, dense=dense)
        sol = (t=interp.t, x=interp.x)
    else
        @time t, x = taylorinteg(NBP_pN_A_J23E_J23M_J2S!, q0, t0, tmax, order, abstol, maxsteps=maxsteps, dense=dense)
        sol = (t=t, x=x)
    end

    #write solution to .jld files
    if output
        filename = string("sseph.jld")
        println("Saving solution to file: $filename")
        jldopen(filename, "w") do file
            addrequire(file, TaylorSeries)
            addrequire(file, TaylorIntegration)
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
            recovered_sol_i = load(filename, varname)
            #check that recovered variable is equal to original variable
            @show recovered_sol_i == sol[ind]
        end
        println("Saved solution")
    end
    return nothing
end
