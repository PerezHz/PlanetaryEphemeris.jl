@doc raw"""
    evaluate_threads!(x::Vector{Taylor1{T}}, δt::T, x0::Vector{T}) where { T <: Number}

Threaded version of `TaylorSeries.evaluate!`.

See also [`TaylorSeries.evaluate!`](@ref).
"""
function evaluate_threads!(x::Vector{Taylor1{T}}, δt::T, x0::Vector{T}) where { T <: Number}

    Threads.@threads for i in eachindex(x)
        x0[i] = evaluate( x[i], δt )
    end

    nothing
end

@doc raw"""
    stepsize_threads(q::Vector{Taylor1{U}}, epsilon::T) where {T <: Real, U <: Number}

Threaded version of `TaylorIntegration.stepsize`.

See also [`TaylorIntegration.stepsize`](@ref) and [`TaylorIntegration._second_stepsize`](@ref).
"""
function stepsize_threads(q::Vector{Taylor1{U}}, epsilon::T) where {T <: Real, U <: Number}
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

@doc raw"""
    stepsize_jz05(q::AbstractArray{Taylor1{U}, N}, epsilon::T) where {T<:Real, U<:Number, N}

First step-size control. See section 3.2 of https://doi.org/10.1080/10586458.2005.10128904.

See also [`stepsize_threads`](@ref) and [`TaylorIntegration.stepsize`](@ref).
"""
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

# Constant timestep method: set timestep equal to 1 day
# function stepsize_threads(q::AbstractArray{Taylor1{U},1}, epsilon::T) where
#         {T<:Real, U<:Number}
#     R = promote_type(typeof(norm(constant_term(q[1]), Inf)), T)
#     h = convert(R, Inf)
#     h = TaylorIntegration.stepsize( q[1], epsilon )
#     return one(h)::R
# end

@doc raw"""
    taylorstep_threads!(f!, t::Taylor1{T}, x::Vector{Taylor1{U}}, dx::Vector{Taylor1{U}}, xaux::Vector{Taylor1{U}},
                        abstol::T, params, parse_eqs::Bool=true) where {T<:Real, U<:Number}
    taylorstep_threads!(f!, t::Taylor1{T}, x::Vector{Taylor1{U}}, dx::Vector{Taylor1{U}}, abstol::T, params,
                        rv::TaylorIntegration.RetAlloc{Taylor1{U}}) where {T<:Real, U<:Number}

Threaded version of `TaylorIntegration.taylorstep`.

See also [`stepsize_threads`](@ref) and [`TaylorIntegration.taylorstep`](@ref).
"""
function taylorstep_threads!(f!, t::Taylor1{T}, x::Vector{Taylor1{U}}, dx::Vector{Taylor1{U}}, xaux::Vector{Taylor1{U}},
                             abstol::T, params) where {T<:Real, U<:Number}

    # Compute the Taylor coefficients
    TaylorIntegration.__jetcoeffs!(Val(false), f!, t, x, dx, xaux, params)

    # Compute the step-size of the integration using `abstol`
    δt = stepsize_threads(x, abstol)

    return δt
end

function taylorstep_threads!(f!, t::Taylor1{T}, x::Vector{Taylor1{U}}, dx::Vector{Taylor1{U}}, abstol::T, params,
                             rv::TaylorIntegration.RetAlloc{Taylor1{U}}) where {T<:Real, U<:Number}

    # Compute the Taylor coefficients
    TaylorIntegration.__jetcoeffs!(Val(true), f!, t, x, dx, params, rv)

    # Compute the step-size of the integration using `abstol`
    δt = stepsize_threads(x, abstol)

    return δt
end

@doc raw"""
    __determine_parsing!(parse_eqs::Bool, f, t, x, dx, params)

Specialized method of `TaylorIntegration._determine_parsing!` to avoid invalidations.

See also [`TaylorIntegration._determine_parsing!`](@ref).
"""
function __determine_parsing!(parse_eqs::Bool, f, t, x, dx, params)

    rv = TaylorIntegration._allocate_jetcoeffs!(t, x, dx, params)

    if parse_eqs
        try
            rv = TaylorIntegration._allocate_jetcoeffs!(Val(f), t, x, dx, params)
            TaylorIntegration.jetcoeffs!(Val(f), t, x, dx, params, rv)
        catch
            @warn("""Unable to use the parsed method of `jetcoeffs!` for `$f`,
            despite of having `parse_eqs=true`, due to some internal error.
            Using `parse_eqs = false`.""")
            parse_eqs = false
        end
    end

    return parse_eqs, rv
end

@doc raw"""
    taylorinteg_threads(f!, q0::Array{U,1}, t0::T, tmax::T, order::Int, abstol::T, Val(true/false),
                        params = nothing; maxsteps::Int=500, parse_eqs::Bool=true) where {T<:Real, U<:Number}

Threaded version of `TaylorIntegration.taylorinteg`.

See also [`TaylorIntegration.taylorinteg`](@ref).
""" taylorinteg_threads

for V in (:(Val{true}), :(Val{false}))
    @eval begin

        function taylorinteg_threads(f!, q0::Array{U, 1}, t0::T, tmax::T, order::Int, abstol::T, ::$V, params = nothing;
                                     maxsteps::Int = 500, parse_eqs::Bool = true) where {T <: Real, U <: Number}

            # Initialize the vector of Taylor1 expansions
            dof = length(q0)
            t = t0 + Taylor1( T, order )
            x = Array{Taylor1{U}}(undef, dof)
            dx = Array{Taylor1{U}}(undef, dof)
            @inbounds for i in eachindex(q0)
                @inbounds x[i] = Taylor1( q0[i], order )
                @inbounds dx[i] = Taylor1( zero(q0[i]), order )
            end

            # Determine if specialized jetcoeffs! method exists
            parse_eqs, rv = __determine_parsing!(parse_eqs, f!, t, x, dx, params)

            if parse_eqs
                # Re-initialize the Taylor1 expansions
                t = t0 + Taylor1( T, order )
                x .= Taylor1.( q0, order )
                return _taylorinteg_threads!(f!, t, x, dx, q0, t0, tmax, abstol, rv, $V(), params, maxsteps = maxsteps)
            else
                return _taylorinteg_threads!(f!, t, x, dx, q0, t0, tmax, abstol, $V(), params, maxsteps = maxsteps)
            end

        end

        function _taylorinteg_threads!(f!, t::Taylor1{T}, x::Array{Taylor1{U}, 1}, dx::Array{Taylor1{U}, 1}, q0::Array{U, 1}, t0::T,
                                       tmax::T, abstol::T, ::$V, params; maxsteps::Int = 500) where {T <: Real, U <: Number}

            # Initialize the vector of Taylor1 expansions
            dof = length(q0)

            # Allocation
            tv = Array{T}(undef, maxsteps+1)
            xv = Array{U}(undef, dof, maxsteps+1)
            if $V == Val{true}
                psol = Array{Taylor1{U}}(undef, dof, maxsteps)
            end
            xaux = Array{Taylor1{U}}(undef, dof)

            # Initial conditions
            @inbounds t[0] = t0
            # x .= Taylor1.(q0, order)
            x0 = deepcopy(q0)
            @inbounds tv[1] = t0
            @inbounds xv[:,1] .= q0
            sign_tstep = copysign(1, tmax-t0)

            # Integration
            nsteps = 1
            while sign_tstep*t0 < sign_tstep*tmax
                δt = taylorstep_threads!(f!, t, x, dx, xaux, abstol, params) # δt is positive!
                # Below, δt has the proper sign according to the direction of the integration
                δt = sign_tstep * min(δt, sign_tstep*(tmax-t0))
                evaluate_threads!(x, δt, x0) # new initial condition
                if $V == Val{true}
                    # Store the Taylor polynomial solution
                    @inbounds psol[:,nsteps] .= deepcopy.(x)
                end
                @inbounds Threads.@threads for i in eachindex(x0)
                    x[i][0] = x0[i]
                    dx[i][0] = zero(x0[i])
                end
                t0 += δt
                @inbounds t[0] = t0
                nsteps += 1
                @inbounds tv[nsteps] = t0
                @inbounds xv[:,nsteps] .= x0
                if nsteps > maxsteps
                    @warn("""
                    Maximum number of integration steps reached; exiting.
                    """)
                    break
                end
            end

            if $V == Val{true}
                return TaylorInterpolant(tv[1], view(tv.-tv[1],1:nsteps), view(transpose(view(psol,:,1:nsteps-1)),1:nsteps-1,:))
            elseif $V == Val{false}
                return view(tv,1:nsteps), view(transpose(view(xv,:,1:nsteps)),1:nsteps,:)
            end
        end

        function _taylorinteg_threads!(f!, t::Taylor1{T}, x::Array{Taylor1{U}, 1}, dx::Array{Taylor1{U}, 1}, q0::Array{U, 1}, t0::T,
                                       tmax::T, abstol::T, rv::TaylorIntegration.RetAlloc{Taylor1{U}}, ::$V, params; maxsteps::Int = 500) where {T <: Real, U <: Number}

            # Initialize the vector of Taylor1 expansions
            dof = length(q0)

            # Allocation of output
            tv = Array{T}(undef, maxsteps+1)
            xv = Array{U}(undef, dof, maxsteps+1)
            if $V == Val{true}
                psol = Array{Taylor1{U}}(undef, dof, maxsteps)
            end

            # Initial conditions
            @inbounds t[0] = t0
            x0 = deepcopy(q0)
            @inbounds tv[1] = t0
            @inbounds xv[:,1] .= q0
            sign_tstep = copysign(1, tmax-t0)

            # Integration
            nsteps = 1
            while sign_tstep*t0 < sign_tstep*tmax
                δt = taylorstep_threads!(f!, t, x, dx, abstol, params, rv) # δt is positive!
                # Below, δt has the proper sign according to the direction of the integration
                δt = sign_tstep * min(δt, sign_tstep*(tmax-t0))
                evaluate_threads!(x, δt, x0) # new initial condition
                if $V == Val{true}
                    # Store the Taylor polynomial solution
                    @inbounds psol[:,nsteps] .= deepcopy.(x)
                end

                Threads.@threads for i in eachindex(x0)
                    x[i][0] = x0[i]
                    dx[i][0] = zero(x0[i])
                end
                t0 += δt
                @inbounds t[0] = t0
                nsteps += 1
                @inbounds tv[nsteps] = t0
                @inbounds xv[:,nsteps] .= x0
                if nsteps > maxsteps
                    @warn("""
                    Maximum number of integration steps reached; exiting.
                    """)
                    break
                end
            end

            if $V == Val{true}
                return TaylorInterpolant(tv[1], view(tv.-tv[1],1:nsteps), view(transpose(view(psol,:,1:nsteps-1)),1:nsteps-1,:))
            elseif $V == Val{false}
                return view(tv,1:nsteps), view(transpose(view(xv,:,1:nsteps)),1:nsteps,:)
            end
        end

    end
end
