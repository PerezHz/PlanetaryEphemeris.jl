"""
    DE430Params{T <: Real}

Parameters for the JPL DE430 planetary ephemeris model.

# Fields

- `N::Int`: total number of bodies.
- `N_bwd::Int`: number of bodies used to compute time-delayed
    tidal interactions.
- `jd0::T`: reference epoch [JDTDB].
- `rv::RetAlloc{Taylor1{T}}`: time-delayed tidal interactions
    buffer.
"""
struct DE430Params{T <: Real}
    N::Int
    N_bwd::Int
    jd0::T
    rv::RetAlloc{Taylor1{T}}
end

# Outer constructor
function DE430Params(jd0::T, q0::AbstractVector{T}, order::Int;
                     parse_eqs::Bool = true) where {T <: Real}
    # Total number of bodies
    N = 11 + 343
    # Number of bodies used to compute time-delayed tidal interactions
    N_bwd = 11
    # Time-delayed tidal interactions buffer
    params_bwd = (N_bwd, jd0)
    t, x, dx = init_expansions(zero(T), q0, order)
    idxs_bwd = union(nbodyind(N, 1:N_bwd),6N+1:6N+13)
    qq_bwd = [Taylor1(constant_term(x[i]), order) for i in idxs_bwd]
    dqq_bwd = [Taylor1(constant_term(dx[i]), order) for i in idxs_bwd]
    _, rv = _determine_parsing!(parse_eqs, NBP_pN_A_J23E_J23M_J2S_threads!, t, qq_bwd,
                                dqq_bwd, params_bwd)

    return DE430Params{T}(N, N_bwd, jd0, rv)
end

# Print method for DE430Params
show(io::IO, x::DE430Params) = print(io, "DE430Params")