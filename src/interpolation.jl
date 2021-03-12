# This file is part of the TaylorIntegration.jl package; MIT licensed
@auto_hash_equals struct TaylorInterpolant{T,U,N}
    t0::T
    t::AbstractVector{T}
    x::AbstractArray{Taylor1{U},N}
    #Inner constructor
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

#outer constructor
function TaylorInterpolant(t0::T, t::AbstractVector{T},
        x::AbstractArray{Taylor1{U},N}) where {T<:Real, U<:Number, N}
    return TaylorInterpolant{T,U,N}(t0, t, x)
end

# return time vector index corresponding to interpolation range
function getinterpindex(tinterp::TaylorInterpolant{T,U,N}, t::V) where {T<:Real, U<:Number, V<:Number, N}
    t00 = constant_term(constant_term(t))
    tmin, tmax = minmax(tinterp.t[end], tinterp.t[1])
    Δt = t-tinterp.t0
    Δt00 = t00-tinterp.t0
    @assert tmin ≤ Δt00 ≤ tmax "Evaluation time outside range of interpolation"
    if Δt00 == tinterp.t[end] # compute solution at final time from last step expansion
        ind = lastindex(tinterp.t)-1
    elseif issorted(tinterp.t) # forward integration
        ind = searchsortedlast(tinterp.t, Δt00)
    elseif issorted(tinterp.t, rev=true) # backward integration
        ind = searchsortedlast(tinterp.t, Δt00, rev=true)
    end
    return ind, Δt
end

# function-like (callability) methods
function (tinterp::TaylorInterpolant{T,U,1})(t::V) where {T<:Real, U<:Number, V<:Number}
    ind, Δt = getinterpindex(tinterp, t)
    δt = Δt-tinterp.t[ind]
    return tinterp.x[ind](δt)
end

function (tinterp::TaylorInterpolant{T,U,2})(t::V) where {T<:Real, U<:Number, V<:Number}
    ind, Δt = getinterpindex(tinterp, t)
    δt = Δt-tinterp.t[ind]
    return tinterp.x[ind,:](δt)
end

# reverse independent variable function
function reverse(tinterp::TaylorInterpolant{T,U,N}) where {T<:Real, U<:Number, N}
    tinterp_rev_t0 = tinterp.t[end]
    tinterp_rev_t = tinterp.t[end:-1:1] .- tinterp_rev_t0
    tinterp_rev_x = vcat(tinterp(tinterp.t[end]+Taylor1(tinterp.x[1].order))', tinterp.x[end:-1:2,:])
    return TaylorInterpolant(tinterp_rev_t0, tinterp_rev_t, tinterp_rev_x)
end
