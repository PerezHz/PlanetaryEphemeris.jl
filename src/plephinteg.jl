# Threaded version of TaylorSeries.evaluate!
function evaluate_threads!(x::Array{Taylor1{T},1}, δt::T,
        x0::Union{Array{T,1},SubArray{T,1}}) where {T<:Number}
    # @assert length(x) == length(x0)
    Threads.@threads for i in eachindex(x, x0)
        x0[i] = evaluate( x[i], δt )
    end
    nothing
end

# Threaded version of TaylorIntegration.stepsize
function stepsize_threads(q::AbstractArray{Taylor1{U},1}, epsilon::T) where
        {T<:Real, U<:Number}
    R = promote_type(typeof(norm(constant_term(q[1]), Inf)), T)
    h = convert(R, Inf)
    Threads.@threads for i in eachindex(q)
        @inbounds hi = TaylorIntegration.stepsize( q[i], epsilon )
        h = min( h, hi )
    end

    # If `isinf(h)==true`, we use the maximum (finite)
    # step-size obtained from all coefficients as above.
    # Note that the time step is independent from `epsilon`.
    if isinf(h)
        h = zero(R)
        Threads.@threads for i in eachindex(q)
            @inbounds hi = TaylorIntegration._second_stepsize(q[i], epsilon)
            h = max( h, hi )
        end
    end
    return h::R
    # h = TaylorIntegration.stepsize( q[1], epsilon )
    # return one(h)::R
end

function taylorstep_threads!(f!, t::Taylor1{T}, x::Vector{Taylor1{U}},
        dx::Vector{Taylor1{U}}, xaux::Vector{Taylor1{U}}, abstol::T, params,
        parse_eqs::Bool=true) where {T<:Real, U<:Number}

    # Compute the Taylor coefficients
    TaylorIntegration.__jetcoeffs!(Val(parse_eqs), f!, t, x, dx, xaux, params)

    # Compute the step-size of the integration using `abstol`
    δt = stepsize_threads(x, abstol)

    return δt
end

function taylorinteg_threads(f!, q0::Array{U,1}, t0::T, tmax::T, order::Int, abstol::T,
        params = nothing; maxsteps::Int=500, parse_eqs::Bool=true, dense::Bool=false) where {T<:Real, U<:Number}

    # Allocation
    tv = Array{T}(undef, maxsteps+1)
    dof = length(q0)
    xv = Array{U}(undef, dof, maxsteps+1)
    if dense
        xv_interp = Array{Taylor1{U}}(undef, dof, maxsteps+1)
    end

    # Initialize the vector of Taylor1 expansions
    t = Taylor1(T, order)
    x = Array{Taylor1{U}}(undef, dof)
    dx = Array{Taylor1{U}}(undef, dof)
    xaux = Array{Taylor1{U}}(undef, dof)
    dx .= Taylor1.(zeros(U), order)

    # Initial conditions
    @inbounds t[0] = t0
    x .= Taylor1.(q0, order)
    x0 = deepcopy(q0)
    @inbounds tv[1] = t0
    @inbounds xv[:,1] .= q0
    sign_tstep = copysign(1, tmax-t0)

    # Determine if specialized jetcoeffs! method exists
    parse_eqs = TaylorIntegration._determine_parsing!(parse_eqs, f!, t, x, dx, params)

    # Integration
    nsteps = 1
    while sign_tstep*t0 < sign_tstep*tmax
        δt = taylorstep_threads!(f!, t, x, dx, xaux, abstol, params, parse_eqs) # δt is positive!
        # δt = TaylorIntegration.taylorstep!(f!, t, x, dx, xaux, abstol, params, parse_eqs) # δt is positive!
        # Below, δt has the proper sign according to the direction of the integration
        # δt = sign_tstep * min(δt, sign_tstep*(tmax-t0))
        δt = sign_tstep * min(δt, sign_tstep*(tmax-t0))
        evaluate_threads!(x, δt, x0) # new initial condition
        # evaluate!(x, δt, x0) # new initial condition
        if dense
            # @inbounds xv_interp[:,nsteps] .= deepcopy(x)
            Threads.@threads for i in eachindex(x0)
                @inbounds xv_interp[i,nsteps] = deepcopy(x[i])
            end
        else
            # @inbounds xv[:,nsteps] .= x0
            Threads.@threads for i in eachindex(x0)
                @inbounds xv[i,nsteps] = x0[i]
            end
        end
        Threads.@threads for i in eachindex(x0)
            @inbounds x[i][0] = x0[i]
            @inbounds dx[i] = Taylor1( zero(x0[i]), order )
        end
        t0 += δt
        @inbounds t[0] = t0
        nsteps += 1
        @inbounds tv[nsteps] = t0
        if nsteps > maxsteps
            @warn("""
            Maximum number of integration steps reached; exiting.
            """)
            break
        end
    end

    if dense
        return TaylorInterpolant(tv[1], view(tv.-tv[1],1:nsteps), view(transpose(view(xv_interp,:,1:nsteps-1)),1:nsteps-1,:))
    else
        return view(tv,1:nsteps), view(transpose(view(xv,:,1:nsteps)),1:nsteps,:)
    end
end
