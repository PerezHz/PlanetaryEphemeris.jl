# This file is part of the TaylorIntegration.jl package; MIT licensed

@doc raw"""
    TaylorInterpolant{T,U,N}
 
Collection of Taylor polynomials that interpolate a dependent variable as a function of 
an independent variable. For example, the ``x``-axis position of the Earth as a function of 
time ``x(t)``; or a lunar core Euler angle as a function of time ``\theta_c(t)``.

# Fields

- `t0::T`: Start time.
- `t::AbstractVector{T}`: Vector of time instances when the timespan of the ``i``-th element of `x` ends and the ``(i+1)``-th element of `x` starts being valid. 
- `x::AbstractArray{Taylor1{U},N}`: Vector of Taylor polynomials that interpolate the dependent variable as a function of the independent variable.
"""
@auto_hash_equals struct TaylorInterpolant{T,U,N}
    t0::T
    t::AbstractVector{T}
    x::AbstractArray{Taylor1{U},N}
    # Inner constructor
    function TaylorInterpolant{T,U,N}(
        t0::T,
        t::AbstractVector{T},
        x::AbstractArray{Taylor1{U},N}
    ) where {T<:Real, U<:Number, N}
        @assert size(x)[1] == length(t)-1
        @assert issorted(t) || issorted(t, rev=true)
        return new{T,U,N}(t0, t, x)
    end
end

# Outer constructor
function TaylorInterpolant(t0::T, t::AbstractVector{T},
        x::AbstractArray{Taylor1{U},N}) where {T<:Real, U<:Number, N}
    return TaylorInterpolant{T,U,N}(t0, t, x)
end

@doc raw"""
    getinterpindex(tinterp::TaylorInterpolant{T,U,N}, t::V) where {T<:Real, U<:Number, V<:Number, N}

Returns the index of `tinterp.t` corresponding to `t` and the time elapsed from `tinterp.t0`
to `t`.
"""
function getinterpindex(tinterp::TaylorInterpolant{T,U,N}, t::V) where {T<:Real, U<:Number, V<:Number, N}
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

# Function-like (callability) methods

@doc raw"""
    (tinterp::TaylorInterpolant{T,U,1})(t::V) where {T<:Real, U<:Number, V<:Number}
    (tinterp::TaylorInterpolant{T,U,2})(t::V) where {T<:Real, U<:Number, V<:Number}

Evaluates `tinterp.x` at time `t`.

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

@doc raw"""
    reverse(tinterp::TaylorInterpolant{T,U,N}) where {T<:Real, U<:Number, N}

Returns a `TaylorInterpolant` object with the same information as `tinterp` but 
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
