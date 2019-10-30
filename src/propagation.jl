function propagate(maxsteps::Int, t0::T, tspan::T;
        output::Bool=true, dense::Bool=false, ephfile::String="sseph.jld",
        dynamics::Function=NBP_pN_A_J23E_J23M_J2S!) where {T<:Real}

    # get initial conditions
    q0 = initialcond(length(μ))

    @show tmax = t0+tspan*yr #final time of integration

    # do integration
    if dense
        @time interp = taylorinteg(dynamics, q0, t0, tmax, order, abstol, maxsteps=maxsteps, dense=dense)
        sol = (t=interp.t[:], x=interp.x[:,:])
    else
        @time t, x = taylorinteg(dynamics, q0, t0, tmax, order, abstol, maxsteps=maxsteps, dense=dense)
        sol = (t=t[:], x=x[:,:])
    end

    #write solution to .jld files
    if output
        println("Saving solution to file: $ephfile")
        jldopen(ephfile, "w") do file
            addrequire(file, TaylorSeries)
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
        println("Saved solution")
    end
    return nothing
end

# # auxiliary function; generates ephemeris file for bodies `i1` through `i2` out of `ephfile`, which stores ephemeris of a full Solar System integration
# # by default, the selected bodies are the Sun, the eight planets, the Moon, and 16 most-massive main-belt asteroids
function save_bodies_eph(ephfile::String, outfile::String, i1::Int=1, i2::Int=27)
    @assert i2 ≥ i1
    t = load(ephfile, "t")
    x = load(ephfile, "x")
    indvec = union(3i1-2:3i2, 3*(N+i1)-2:3*(N+i2))
    # ny = Int(floor( (t[end]-t[1])/yr ))
    jldopen(outfile, "w") do file
        addrequire(file, TaylorSeries)
        write(file, "t", t)
        write(file, "x", x[:, indvec])
    end
end