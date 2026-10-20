"""
    DE430Constants{T <: Real}

Physical parameters of the JPL DE430 dynamical model, i.e. the mass parameters and
the extended-body, rotational and tidal parameters of the Sun, the Earth and the Moon
used by [`NBP_pN_A_J23E_J23M_J2S_threads!`](@ref) and [`DE430!`](@ref).

# Fields

- `N_ext::Int`: number of bodies in extended-body accelerations.
- `μ::Vector{T}`: mass parameters [au^3/day^2].
- `UJ_interaction::Matrix{Bool}`: whether the extended body `j` (column) interacts
    with the point mass `i` (row).
- `α_p_sun::T`, `δ_p_sun::T`: right ascension and declination of the Sun's rotation
    pole [deg].
- `RE_au::T`: Earth's radius [au].
- `J2E::T`, `J2EDOT::T`: Earth's second zonal harmonic and its rate of change [1/yr].
- `JSEM::Matrix{T}`: zonal harmonics `J_n` times the radius of the body [au] to the
    `n`-th power (row: body, column: degree).
- `n1SEM::Vector{Int}`: degree of the zonal harmonics expansion of each extended body.
- `n2M::Int`: degree of the Moon's tesseral harmonics expansion.
- `CM::Matrix{T}`, `SM::Matrix{T}`: Moon's tesseral harmonics `C_{nm}`, `S_{nm}` times
    the Moon's radius [au] to the `n`-th power.
- `R_moon::T`: Moon's radius [au].
- `k_2M::T`: Moon's potential Love number.
- `τ_M::T`: time-lag of the lunar solid-body tide [day].
- `n_moon::T`: Moon's mean motion [rad/day].
- `ITM_und::Matrix{T}`: Moon's undistorted total moment of inertia.
- `I_c::Matrix{T}`: Moon's core moment of inertia.
- `k_ν::T`: friction coefficient between the Moon's core and mantle.
- `C_c_m_A_c::T`: numerical factor of the torque on the Moon's mantle due to
    the interaction between the core and mantle.
- `k_20E::T`, `k_21E::T`, `k_22E::T`: Earth's potential Love numbers for
    long-period, diurnal and semi-diurnal deformations.
- `τ_0p::T`, `τ_1p::T`, `τ_2p::T`: Earth's orbital time-lags for long-period,
    diurnal and semi-diurnal deformations [day].
- `τ_0::T`, `τ_1::T`, `τ_2::T`: Earth's rotational time-lags for long-period,
    diurnal and semi-diurnal deformations [day].
- `tid_num_coeff::T`: overall numerical factor of the tidal acceleration of the
    Moon, `1.5 * (1 + μ[mo] / μ[ea])`.
"""
struct DE430Constants{T <: Real}
    N_ext::Int
    μ::Vector{T}
    UJ_interaction::Matrix{Bool}
    α_p_sun::T
    δ_p_sun::T
    RE_au::T
    J2E::T
    J2EDOT::T
    JSEM::Matrix{T}
    n1SEM::Vector{Int}
    n2M::Int
    CM::Matrix{T}
    SM::Matrix{T}
    R_moon::T
    k_2M::T
    τ_M::T
    n_moon::T
    ITM_und::Matrix{T}
    I_c::Matrix{T}
    k_ν::T
    C_c_m_A_c::T
    k_20E::T
    k_21E::T
    k_22E::T
    τ_0p::T
    τ_1p::T
    τ_2p::T
    τ_0::T
    τ_1::T
    τ_2::T
    tid_num_coeff::T
end

# Outer constructor (DE430 values in src/constants.jl)
DE430Constants{T}() where {T <: Real} = DE430Constants{T}(
    N_ext, T.(μ), copy(UJ_interaction), α_p_sun, δ_p_sun, RE_au, J2E, J2EDOT,
    T.(JSEM), copy(n1SEM), n2M, T.(CM), T.(SM), R_moon, k_2M, τ_M, n_moon,
    T.(ITM_und), T.(I_c), k_ν, C_c_m_A_c, k_20E, k_21E, k_22E, τ_0p, τ_1p, τ_2p,
    τ_0, τ_1, τ_2, tid_num_coeff
)

DE430Constants() = DE430Constants{Float64}()

# Print method for DE430Constants
show(io::IO, x::DE430Constants) = print(io, "DE430Constants")

"""
    DE430Params{T <: Real}

Parameters for the JPL DE430 planetary ephemeris model.

# Fields

- `N::Int`: total number of bodies.
- `N_bwd::Int`: number of bodies used to compute time-delayed
    tidal interactions.
- `jd0::T`: reference epoch [JDTDB].
- `consts::DE430Constants{T}`: physical parameters of the model.
- `rv::RetAlloc{Taylor1{T}}`: time-delayed tidal interactions
    buffer.
"""
struct DE430Params{T <: Real}
    N::Int
    N_bwd::Int
    jd0::T
    consts::DE430Constants{T}
    rv::RetAlloc{Taylor1{T}}
end

# Outer constructor
function DE430Params(jd0::T, q0::AbstractVector{T}, order::Int;
                     parse_eqs::Bool = true,
                     consts::DE430Constants{T} = DE430Constants{T}()) where {T <: Real}
    # Total number of bodies
    N = 11 + 343
    # Number of bodies used to compute time-delayed tidal interactions
    N_bwd = 11
    # Time-delayed tidal interactions buffer
    params_bwd = (N_bwd, jd0, consts)
    t, x, dx = init_expansions(zero(T), q0, order)
    idxs_bwd = union(nbodyind(N, 1:N_bwd),6N+1:6N+13)
    qq_bwd = [Taylor1(constant_term(x[i]), order) for i in idxs_bwd]
    dqq_bwd = [Taylor1(constant_term(dx[i]), order) for i in idxs_bwd]
    _, rv = _determine_parsing!(parse_eqs, NBP_pN_A_J23E_J23M_J2S_threads!, t, qq_bwd,
                                dqq_bwd, params_bwd)

    return DE430Params{T}(N, N_bwd, jd0, consts, rv)
end

# Print method for DE430Params
show(io::IO, x::DE430Params) = print(io, "DE430Params")
