# This file is part of the PlanetaryEphemeris.jl package; MIT licensed

@doc raw"""
    TaylorInterpolant{T, U, N}

Collection of Taylor polynomials that interpolate a dependent variable as a function of
an independent variable. For example, the ``x``-axis position of the Earth as a function of
time ``x(t)``; or a lunar core Euler angle as a function of time ``\theta_c(t)``.

# Fields

- `t0::T`: Start time.
- `t::AbstractVector{T}`: Vector of time instances when the timespan of the ``i``-th element of `x` ends and the ``(i+1)``-th element of `x` starts being valid.
- `x::AbstractArray{Taylor1{U},N}`: Array of Taylor polynomials that interpolate the dependent variable as a function of the independent variable.
"""
@auto_hash_equals struct TaylorInterpolant{T<:Number, U<:Number, N, VT<:AbstractVector{T}, X<:AbstractArray{Taylor1{U}, N}}
    t0::T
    t::VT
    x::X
    # Inner constructor
    function TaylorInterpolant{T, U, N, VT, X}(t0::T, t::VT, x::X) where {T<:Number, U<:Number, N, VT<:AbstractVector{T}, X<:AbstractArray{Taylor1{U}, N}}
        @assert size(x)[1] == length(t)-1
        @assert issorted(t) || issorted(t, rev = true)
        return new{T, U, N, VT, X}(t0, t, x)
    end
end

const TaylorInterpCallingArgs{T,U} = Union{T, U, Taylor1{U}, TaylorN{U}, Taylor1{TaylorN{U}}} where {T,U}

# Outer constructors
function TaylorInterpolant{T, U, N}(t0::T, t::VT, x::X) where {T<:Number, U<:Number, N, VT<:AbstractVector{T}, X<:AbstractArray{Taylor1{U}, N}}
    return TaylorInterpolant{T, U, N, Vector{T}, Array{Taylor1{U},N}}(t0, convert(Vector{T},t), convert(Array{Taylor1{U},N},x))
end

function TaylorInterpolant(t0::T, t::VT, x::X) where {T<:Number, U<:Number, N, VT<:AbstractVector{T}, X<:AbstractArray{Taylor1{U}, N}}
    return TaylorInterpolant{T, U, N, VT, X}(t0, t, x)
end

function TaylorInterpolant(t0::T, t::SubArray{T, 1}, x::SubArray{Taylor1{U}, N}) where {T<:Number, U<:Number, N}
    return TaylorInterpolant{T, U, N}(t0, t, x)
end

# Definition of zero TaylorInterpolant
function zero(::Type{TaylorInterpolant{T, U, N, VT, X}}) where {T <: Number, U <: Number, N, VT <: AbstractVector{T}, X <: AbstractArray{Taylor1{U}, N}}
    return TaylorInterpolant(zero(T), zeros(T, 1), Array{Taylor1{U}, N}(undef, zeros(Int, N)...))
end

function iszero(x::TaylorInterpolant{T, U, N, VT, X}) where {T <: Number, U <: Number, N, VT <: AbstractVector{T}, X <: AbstractArray{Taylor1{U}, N}}
    return x == zero(TaylorInterpolant{T, U, N, VT, X})
end

# Custom print
function show(io::IO, interp::T) where {U, V, N, T<:TaylorInterpolant{U,V,N}}
    t_range = minmax(interp.t0 + interp.t[1], interp.t0 + interp.t[end])
    S = eltype(interp.x)
    if isone(N)
        print(io, "t: ", t_range, ", x: 1 ", S, " variable")
    else
        L = size(interp.x, 2)
        print(io, "t: ", t_range, ", x: ", L, " ", S, " variables")
    end
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
function getinterpindex(tinterp::TaylorInterpolant{T,U,N}, t::TT) where {T,U,N,TT<:TaylorInterpCallingArgs{T,U}}
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

    # Time since the start of the ind-th timestep
    δt = Δt-tinterp.t[ind]

    # Return index and elapsed time since i-th timestep
    return ind::Int, δt::TT
end

@doc raw"""
    numberofbodies(interp::TaylorInterpolant{T, U, 2}) where {T, U}

Return the number of bodies saved in a `TaylorInterpolant` produced by `PlanetaryEphemeris`.
"""
numberofbodies(L::Int) = (L - 13) ÷ 6
numberofbodies(v::Vector{T}) where {T} = numberofbodies(length(v))
numberofbodies(m::Matrix{T}) where {T} = numberofbodies(size(m, 2))
numberofbodies(interp::TaylorInterpolant) = numberofbodies(size(interp.x, 2))

# Function-like (callability) methods

@doc raw"""
    (tinterp::TaylorInterpolant{T, U, 1})(t::T) where {T, U}
    (tinterp::TaylorInterpolant{T, U, 1})(t::TT) where {T, U, TT<:TaylorInterpCallingArgs{T,U}}
    (tinterp::TaylorInterpolant{T, U, 2})(t::T) where {T, U}
    (tinterp::TaylorInterpolant{T, U, 2})(t::TT) where {T, U, TT<:TaylorInterpCallingArgs{T,U}}

Evaluate `tinterp.x` at time `t`.

See also [`getinterpindex`](@ref).
"""
function (tinterp::TaylorInterpolant{T, U, 1})(t::T) where {T, U}
    # Get index of tinterp.x that interpolates at time t
    ind::Int, δt::T = getinterpindex(tinterp, t)
    # Evaluate tinterp.x[ind] at δt
    return (tinterp.x[ind])(δt)::U
end
function (tinterp::TaylorInterpolant{T, U, 1})(t::TT) where {T, U, TT<:TaylorInterpCallingArgs{T,U}}
    # Get index of tinterp.x that interpolates at time t
    ind::Int, δt::TT = getinterpindex(tinterp, t)
    # Evaluate tinterp.x[ind] at δt
    return (tinterp.x[ind])(δt)::TT
end

function (tinterp::TaylorInterpolant{T, U, 2})(t::T) where {T, U}
    # Get index of tinterp.x that interpolates at time t
    ind::Int, δt::T = getinterpindex(tinterp, t)
    # Evaluate tinterp.x[ind] at δt
    return view(tinterp.x, ind, :)(δt)::Vector{U}
end
function (tinterp::TaylorInterpolant{T, U, 2})(t::TT) where {T, U, TT<:TaylorInterpCallingArgs{T,U}}
    # Get index of tinterp.x that interpolates at time t
    ind::Int, δt::TT = getinterpindex(tinterp, t)
    # Evaluate tinterp.x[ind] at δt
    return view(tinterp.x, ind, :)(δt)::Vector{TT}
end

function (tinterp::TaylorInterpolant{T, U, 2})(target::Int, observer::Int, t::T)  where {T, U}
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

(tinterp::TaylorInterpolant{T,U,2})(target::Int, t::TT) where {T, U, TT<:TaylorInterpCallingArgs{T,U}} = tinterp(target, 0, t)

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

@doc raw"""
    selecteph(eph::TaylorInterpolant, bodyind::Union{Int, AbstractVector{Int}}, t0::T = eph.t0,
              tf::T = eph.t0 + eph.t[end]; euler::Bool = false, ttmtdb::Bool = false) where {T <: Real}

Return a subset of `eph` with only the ephemeris of bodies `bodyind` in timerange `[t0, tf]`.
The keyword arguments allow to include lunar euler angles and/or TT-TDB.
"""
function selecteph(eph::TaylorInterpolant, bodyind::Union{Int, AbstractVector{Int}}, t0::T = eph.t0,
                   tf::T = eph.t0 + eph.t[end]; euler::Bool = false, ttmtdb::Bool = false) where {T <: Real}
    # Times
    tmin, tmax = minmax(eph.t0, eph.t0 + eph.t[end])
    @assert tmin ≤ t0 ≤ tf ≤ tmax "$tmin ≤ t0 ≤ tf ≤ $tmax"
    if issorted(eph.t)
        i_0 = searchsortedlast(eph.t, t0)
        i_f = searchsortedfirst(eph.t, tf)
    else
        i_0 = searchsortedlast(eph.t, tf, rev = true)
        i_f = searchsortedfirst(eph.t, t0, rev = true)
    end
    # Degrees of freedom
    N = numberofbodies(eph)
    @assert all(bodyind .< N) "bodyind .< $N"
    idxs = nbodyind(N, bodyind)
    if euler
        idxs = vcat(idxs, 6N+1:6N+12)
    end
    if ttmtdb
        idxs = vcat(idxs, 6N+13)
    end
    # Views
    t = view(eph.t, i_0:i_f)
    x = view(eph.x, i_0:i_f-1, idxs)
    # New TaylorInterpolant
    return TaylorInterpolant(eph.t0, t, x)
end

function join(bwd::TaylorInterpolant{T, U, 2}, fwd::TaylorInterpolant{T, U, 2}) where {T, U}
    @assert bwd.t0 == fwd.t0 "Initial time must be the same for both TaylorInterpolant"
    order_bwd = get_order(bwd.x[1, 1])
    order_fwd = get_order(fwd.x[1, 1])
    @assert order_bwd == order_fwd "Expansion order must be the same for both TaylorInterpolant"

    t0 = bwd.t0 + bwd.t[end]
    t1 = abs.(bwd.t)
    t = vcat(t1, t1[end] .+ fwd.t[2:end])
    x = vcat(reverse(bwd.x, dims = 1), fwd.x)

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
