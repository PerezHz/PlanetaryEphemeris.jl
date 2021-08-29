# Threaded version of TaylorSeries.evaluate!
function evaluate_threads!(x::Array{Taylor1{T},1}, δt::T,
        x0::Union{Array{T,1},SubArray{T,1}}) where {T<:Number}
    # @assert length(x) == length(x0)
    Threads.@threads for i in eachindex(x, x0)
        x0[i] = evaluate( x[i], δt )
    end
    nothing
end

#### Threaded version of TaylorIntegration.stepsize
function stepsize_threads(q::AbstractArray{Taylor1{U},1}, epsilon::T) where
        {T<:Real, U<:Number}
    R = promote_type(typeof(norm(constant_term(q[1]), Inf)), T)
    h = convert(R, Inf)
    #= Threads.@threads =# for i in eachindex(q)
        @inbounds hi = TaylorIntegration.stepsize( q[i], epsilon )
        h = min( h, hi )
    end
    # If `isinf(h)==true`, we use the maximum (finite)
    # step-size obtained from all coefficients as above.
    # Note that the time step is independent from `epsilon`.
    if isinf(h)
        h = zero(R)
        #= Threads.@threads =# for i in eachindex(q)
            @inbounds hi = TaylorIntegration._second_stepsize(q[i], epsilon)
            h = max( h, hi )
        end
    end
    return h::R
end

# first step-size control from Jorba & Zou, 2005
function stepsize_jz05(q::AbstractArray{Taylor1{U}, N}, epsilon::T) where
        {T<:Real, U<:Number, N}
    nbodies = (length(q)-13)÷6
    q0_norminf = norm(constant_term.(q[1:6nbodies]), Inf)
    pred = epsilon*q0_norminf ≤ epsilon

    if pred
        p_jz05 = Int(ceil(-0.5log(epsilon)+1)) # atol
    else
        p_jz05 = Int(ceil(-0.5log(epsilon)+1)) # rtol
    end

    order = min(p_jz05, q[1].order) # q[1].order
    ordm1 = order-1
    invorder = 1/order
    invordm1 = 1/ordm1
    qordm1_norminf = norm(getcoeff.(q[1:6nbodies], ordm1), Inf)
    qorder_norminf = norm(getcoeff.(q[1:6nbodies], order), Inf)

    if pred
        ρ_ordm1 = ( 1/qordm1_norminf )^invordm1
        ρ_order = ( 1/qorder_norminf )^invorder
    else
        ρ_ordm1 = ( q0_norminf/qordm1_norminf )^invordm1
        ρ_order = ( q0_norminf/qorder_norminf )^invorder
    end
    ρ = min(ρ_ordm1, ρ_order)
    return ρ*exp(-2.0)
end

# ##### Constant timestep method: set timestep equal to 1 day
# function stepsize_threads(q::AbstractArray{Taylor1{U},1}, epsilon::T) where
#         {T<:Real, U<:Number}
#     R = promote_type(typeof(norm(constant_term(q[1]), Inf)), T)
#     h = convert(R, Inf)
#     h = TaylorIntegration.stepsize( q[1], epsilon )
#     return one(h)::R
# end

function taylorstep_threads!(f!, t::Taylor1{T}, x::Vector{Taylor1{U}},
        dx::Vector{Taylor1{U}}, xaux::Vector{Taylor1{U}}, abstol::T, params,
        parse_eqs::Bool=true) where {T<:Real, U<:Number}

    # Compute the Taylor coefficients
    TaylorIntegration.__jetcoeffs!(Val(parse_eqs), f!, t, x, dx, xaux, params)
    # @time TaylorIntegration.__jetcoeffs!(Val(parse_eqs), f!, t, x, dx, xaux, params)
    # @time TaylorIntegration.__jetcoeffs!(Val(false), f!, t, x, dx, xaux, params)
    # @time TaylorIntegration.__jetcoeffs!(Val(true), f!, t, x, dx, xaux, params)

    # Compute the step-size of the integration using `abstol`
    δt = stepsize_threads(x, abstol)
    # δt = stepsize_jz05(x, abstol)

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

    @show parse_eqs

    # Integration
    nsteps = 1
    while sign_tstep*t0 < sign_tstep*tmax
        δt = taylorstep_threads!(f!, t, x, dx, xaux, abstol, params, parse_eqs) # δt is positive!
        # δt = TaylorIntegration.taylorstep!(f!, t, x, dx, xaux, abstol, params, parse_eqs) # δt is positive!
        # Below, δt has the proper sign according to the direction of the integration
        δt = sign_tstep * min(δt, sign_tstep*(tmax-t0))
        evaluate_threads!(x, δt, x0) # new initial condition
        # evaluate!(x, δt, x0) # new initial condition
        t0 += δt
        @inbounds t[0] = t0
        nsteps += 1
        @inbounds tv[nsteps] = t0
        if dense
            # @inbounds xv_interp[:,nsteps-1] .= deepcopy(x)
            Threads.@threads for i in eachindex(x0)
                @inbounds xv_interp[i,nsteps-1] = deepcopy(x[i])
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
