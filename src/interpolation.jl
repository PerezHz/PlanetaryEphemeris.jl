# This file is part of the TaylorIntegration.jl package; MIT licensed

@doc raw"""
    TaylorInterpolant{T, U, N}
 
Collection of Taylor polynomials that interpolate a dependent variable as a function of 
an independent variable. For example, the ``x``-axis position of the Earth as a function of 
time ``x(t)``; or a lunar core Euler angle as a function of time ``\theta_c(t)``.

# Fields

- `t0::T`: Start time.
- `t::Vector{T}`: Vector of time instances when the timespan of the ``i``-th element of `x` ends and the ``(i+1)``-th element of `x` starts being valid. 
- `x::Array{Taylor1{U},N}`: Vector of Taylor polynomials that interpolate the dependent variable as a function of the independent variable.
"""
@auto_hash_equals struct TaylorInterpolant{T, U, N}
    t0::T
    t::Vector{T}
    x::Array{Taylor1{U}, N}
    # Inner constructor
    function TaylorInterpolant{T, U, N}(t0::T, t::Vector{T}, x::Array{Taylor1{U}, N}) where {T <: Real, U <: Number, N}
        @assert size(x)[1] == length(t)-1
        @assert issorted(t) || issorted(t, rev = true)
        return new{T, U, N}(t0, t, x)
    end
end

# Outer constructors
function TaylorInterpolant(t0::T, t::Vector{T}, x::Array{Taylor1{U}, N}) where {T <: Real, U <: Number, N}
    return TaylorInterpolant{T, U, N}(t0, t, x)
end

function TaylorInterpolant(t0::T, t::SubArray{T, 1}, x::SubArray{Taylor1{U}, N}) where {T <: Real, U <: Number, N}
    return TaylorInterpolant{T, U, N}(t0, t.parent[t.indices...], x.parent[x.indices...])
end 

# Custom print 
function show(io::IO, interp::TaylorInterpolant{T, U, 2}) where {T, U}
    t_range = minmax(interp.t0 + interp.t[1], interp.t0 + interp.t[end])
    N = size(interp.x, 2)
    S = eltype(interp.x)
    print(io, "t: ", t_range, ", x: ", N, " ", S, " variables")
end 

@doc raw"""
    convert(::Type{T}, interp::TaylorInterpolant) where {T <: Real}

Convert `inter.t0`, `inter.t` and coefficients of `interp.x` to type `T`. 
"""
function convert(::Type{T}, interp::TaylorInterpolant) where {T <: Real}
    return TaylorInterpolant(
        T(interp.t0), 
        T.(interp.t), 
        map( x -> Taylor1( T.(x.coeffs) ), interp.x)
    )
end 

@doc raw"""
    getinterpindex(tinterp::TaylorInterpolant{T,U,N}, t::V) where {T<:Real, U<:Number, V<:Number, N}

Return the index of `tinterp.t` corresponding to `t` and the time elapsed from `tinterp.t0`
to `t`.
"""
function getinterpindex(tinterp::TaylorInterpolant{T, U, N}, t::V) where {T<:Real, U<:Number, V<:Number, N}
    t00 = constant_term(constant_term(t))                # Current time
    tmin, tmax = minmax(tinterp.t[end], tinterp.t[1])    # Min and max time in tinterp
    Δt = t-tinterp.t0                                    # Time since start of tinterp
    Δt00 = t00-tinterp.t0                                # Time since start of tinterp

    @assert tmin ≤ Δt00 ≤ tmax "Evaluation time outside range of interpolation"

    if Δt00 == tinterp.t[end]        # Compute solution at final time from last step expansion
        ind = lastindex(tinterp.t)-1
    elseif issorted(tinterp.t)       # Forward integration
        ind = searchsortedlast(tinterp.t, Δt00)
    elseif issorted(tinterp.t, rev=true) # Backward integration
        ind = searchsortedlast(tinterp.t, Δt00, rev=true)
    end

    # Return index and elapsed time
    return ind, Δt
end

numberofbodies(interp::TaylorInterpolant{T, U, 2}) where {T, U} = numberofbodies(size(interp.x, 2))

# Function-like (callability) methods

@doc raw"""
    (tinterp::TaylorInterpolant{T,U,1})(t::V) where {T<:Real, U<:Number, V<:Number}
    (tinterp::TaylorInterpolant{T,U,2})(t::V) where {T<:Real, U<:Number, V<:Number}
    (tinterp::TaylorInterpolant{T,U,2})(target::Int, t::V) where {T<:Real, U<:Number, V<:Number}
    (tinterp::TaylorInterpolant{T,U,2})(target::Int, observer::Int, t::V) where {T<:Real, U<:Number, V<:Number}

Evaluate `tinterp.x` at time `t`.

See also [`getinterpindex`](@ref).
"""
function (tinterp::TaylorInterpolant{T,U,1})(t::V) where {T<:Real, U<:Number, V<:Number}
    # Get index of tinterp.x that interpolates at time t
    ind, Δt = getinterpindex(tinterp, t)
    # Time since the start of the ind-th timespan
    δt = Δt-tinterp.t[ind]
    # Evaluate tinterp.x[ind] at δt
    return tinterp.x[ind](δt)
end

function (tinterp::TaylorInterpolant{T,U,2})(t::V) where {T<:Real, U<:Number, V<:Number}
    # Get index of tinterp.x that interpolates at time t
    ind, Δt = getinterpindex(tinterp, t)
    # Time since the start of the ind-th timespan
    δt = Δt-tinterp.t[ind]
    # Evaluate tinterp.x[ind] at δt
    return tinterp.x[ind,:](δt)
end

function (tinterp::TaylorInterpolant{T,U,2})(target::Int, observer::Int, t::V) where {T<:Real, U<:Number, V<:Number}
    # Number of bodies in tinterp
    N = numberofbodies(tinterp)
    # Ephemeris at time t
    eph_t = tinterp(t)
    # Relative state vector 
    if observer == 0
        return eph_t[nbodyind(N, target)]
    else
        return eph_t[nbodyind(N, target)] - eph_t[nbodyind(N, observer)]
    end 
end

(tinterp::TaylorInterpolant{T,U,2})(target::Int, t::V) where {T<:Real, U<:Number, V<:Number} = tinterp(target, 0, t)

@doc raw"""
    reverse(tinterp::TaylorInterpolant{T,U,N}) where {T<:Real, U<:Number, N}

Return a `TaylorInterpolant` object with the same information as `tinterp` but 
the independent variable reversed.

See also [`TaylorInterpolant`](@ref).
"""
function reverse(tinterp::TaylorInterpolant{T,U,N}) where {T<:Real, U<:Number, N}
    # tinterp end time is the new start time 
    tinterp_rev_t0 = tinterp.t[end]
    # reverse independent variable vector tinterp.t
    tinterp_rev_t = tinterp.t[end:-1:1] .- tinterp_rev_t0
    # reverse dependent variable vector tinterp.x
    tinterp_rev_x = vcat(tinterp(tinterp.t[end]+Taylor1(tinterp.x[1].order))', tinterp.x[end:-1:2,:])
    # Return reversed TaylorInterpolant
    return TaylorInterpolant(tinterp_rev_t0, tinterp_rev_t, tinterp_rev_x)
end

function join(bwd::TaylorInterpolant{T, U, 2}, fwd::TaylorInterpolant{T, U, 2}) where {T, U}
    @assert bwd.t0 == fwd.t0 "Initial time must be the same for both TaylorInterpolant"
    order_bwd = get_order(bwd.x[1, 1])
    order_fwd = get_order(fwd.x[1, 1])
    @assert order_bwd == order_fwd "Expansion order must be the same for both TaylorInterpolant"

    t0 = pre.t0
    t = vcat(reverse(pre.t), post.t[2, :])
    x = vcat(reverse(pre.x, dims = 1), post.x)

    return TaylorInterpolant(t0, t, x)
end 

@doc raw"""
    kmsec2auday(pv)
Convert a `[x, y, z, v_x, v_y, v_z]` state vector from km, km/sec to au, au/day.
See also [`auday2kmsec`](@ref).
"""
function kmsec2auday(pv)
    pv /= au          # (km, km/sec) -> (au, au/sec)
    pv[4:6] *= daysec # (au, au/sec) -> (au, au/day)
    return pv
end

@doc raw"""
    auday2kmsec(pv)
Convert a `[x, y, z, v_x, v_y, v_z]` state vector from au, au/day to km, km/sec.
See also [`kmsec2auday`](@ref).
"""
function auday2kmsec(pv)
    pv *= au          # (au, au/day) -> (km, km/day)
    pv[4:6] /= daysec # (km, km/day) -> (km, km/sec)
    return pv
end
