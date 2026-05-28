# This file is part of the PlanetaryEphemeris.jl package; MIT licensed

timespan(sol::TaylorSolution) = minmax(sol.t[1], sol.t[end])

@doc raw"""
    numberofbodies(sol::TaylorSolution)

Return the number of bodies saved in a dense `TaylorSolution` produced by
`PlanetaryEphemeris`.
"""
numberofbodies(L::Int) = (L - 13) ÷ 6
numberofbodies(v::AbstractVector{T}) where {T} = numberofbodies(length(v))
numberofbodies(m::AbstractMatrix{T}) where {T} = numberofbodies(size(m, 2))
numberofbodies(sol::TaylorSolution) = numberofbodies(size(sol.p, 2))

"""
    (eph::TaylorSolution)([target [, observer] ,] t)

Evaluate `eph` at time `t`. If `target` and `observer` are given,
return the state of `target` relative to `observer`.
"""
function (eph::TaylorSolution{T,U,2})(target::Int, observer::Int, t) where {T,U}
    N = numberofbodies(eph)
    eph_t = eph(t)
    if observer == 0
        return eph_t[nbodyind(N, target)]
    else
        return eph_t[nbodyind(N, target)] - eph_t[nbodyind(N, observer)]
    end
end

(eph::TaylorSolution{T,U,2})(target::Int, t) where {T,U} = eph(target, 0, t)

"""
    selecteph(eph::TaylorSolution, bodyind [, t0, tf]; kwargs...)

Return a subset of `eph` containing only the ephemeris of the `bodyind`-th
bodies and spanning `[t0, tf]`.

# Keyword arguments

- `euler::Bool`: whether to include lunar euler angles (default: `false`).
- `ttmtdb::Bool`: whether to include TT-TDB (default: `false`).
"""
function selecteph(eph::TaylorSolution, bodyind::Union{Int, AbstractVector{Int}},
                   t0::T = first(timespan(eph)), tf::T = last(timespan(eph));
                   euler::Bool = false, ttmtdb::Bool = false) where {T <: Real}
    isnothing(eph.p) && error("selecteph requires a dense TaylorSolution")

    tmin, tmax = timespan(eph)
    @assert tmin ≤ t0 < tf ≤ tmax "$tmin ≤ t0 ≤ tf ≤ $tmax"

    N = numberofbodies(eph)
    bodyinds = bodyind isa Int ? (bodyind,) : bodyind
    @assert all(≤(N), bodyinds) "All bodyind must be smaller than $N"

    cols = nbodyind(N, bodyind)
    if euler
        cols = vcat(cols, 6N+1:6N+12)
    end
    if ttmtdb
        cols = vcat(cols, 6N+13)
    end

    if issorted(eph.t)
        j0 = searchsortedlast(eph.t, t0)
        jf = searchsortedfirst(eph.t, tf)
    else
        j0 = searchsortedlast(eph.t, tf, rev = true)
        jf = searchsortedfirst(eph.t, t0, rev = true)
    end
    j0 = max(firstindex(eph.t), j0)
    jf = min(lastindex(eph.t), jf)
    @assert j0 < jf "No TaylorSolution steps overlap the requested time span"

    t = view(eph.t, j0:jf)
    p = view(eph.p, j0:jf-1, cols)

    return TaylorSolution(collect(t), collect(p))
end

"""
    kmsec2auday(pv)
Convert a cartesian state vector from [km, km/sec] to [au, au/day].
See also [`auday2kmsec`](@ref).
"""
function kmsec2auday(x::AbstractVector)
    k = daysec / au
    y = [x[1] / au, x[2] / au, x[3] / au, x[4] * k, x[5] * k, x[6] * k]
    return y
end

"""
    auday2kmsec(pv)
Convert a cartesian state vector from [au, au/day] to [km, km/sec].
See also [`kmsec2auday`](@ref).
"""
function auday2kmsec(x::AbstractVector)
    k = au / daysec
    y = [x[1] * au, x[2] * au, x[3] * au, x[4] * k, x[5] * k, x[6] * k]
    return y
end
