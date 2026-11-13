# This file is part of the PlanetaryEphemeris.jl package; MIT licensed

# Tests of the JPL DE430/431 Earth and Moon orientation models
# Equation and table numbers refer to Folkner et al. (2014), IPN Progress Report 42-196

using PlanetaryEphemeris
using TaylorSeries
using LinearAlgebra
using Test

# Helpers

# Compare numbers or Taylor1 coefficients (absolute tolerance)
_ori_close(a::Number, b::Number; atol) = abs(a - b) ≤ atol
_ori_close(a::Taylor1, b::Taylor1; atol) = maximum(abs, a.coeffs .- b.coeffs) ≤ atol
_ori_allclose(A, B; atol) = all(_ori_close(a, b; atol) for (a, b) in zip(A, B))

# Horner evaluation of c[1] + c[2]x + c[3]x^2 + ...
function _ori_poly(x, c)
    s = zero(x) + c[end]
    for k in length(c)-1:-1:1
        s = s*x + c[k]
    end
    return s
end

# Rotation matrices, eqs. (11)-(13)
function _ori_Rx(a)
    s, c = sin(a), cos(a); o, z = one(a), zero(a)
    return [o z z; z c s; z -s c]
end
function _ori_Ry(a)
    s, c = sin(a), cos(a); o, z = one(a), zero(a)
    return [c z -s; z o z; s z c]
end
function _ori_Rz(a)
    s, c = sin(a), cos(a); o, z = one(a), zero(a)
    return [c s z; -s c z; z z o]
end

# Independent implementation of eqs. (21)-(26) with the published values in arcseconds
# (nutation amplitudes and Ω rate as in IAU 1980, see the review notes). `t` is the
# TDB time in days since J2000.0; `as2rad` converts arcseconds to radians.
function _ori_ref_angles(t, as2rad)
    T = t/36525
    Ω = _ori_poly(T, [450160.280, -6962890.539, 7.455, 0.008])*as2rad      # eq. (21)
    Δψ = -17.1996*sin(Ω)*as2rad                                            # eq. (22)
    Δϵ = 9.2025*cos(Ω)*as2rad
    ϵ̄ = _ori_poly(T, [84381.448, -46.815, -0.00059, 0.001813])*as2rad    # eq. (24)
    ζ = _ori_poly(T, [0.0, 2306.2181, 0.30188, 0.017998])*as2rad           # eq. (26)
    θ = _ori_poly(T, [0.0, 2004.3109, -0.42665, -0.041833])*as2rad
    z = _ori_poly(T, [0.0, 2306.2181, 1.09468, 0.018203])*as2rad
    Φx = (5.6754203322893470E-03 + (t/365.25)*2.7689915574483550E-04)*as2rad  # eq. (25), Table 10
    Φy = (-1.7022656914989530E-02 + (t/365.25)*(-1.2118591216559240E-03))*as2rad
    return (; Ω, Δψ, Δϵ, ϵ̄, ζ, θ, z, Φx, Φy)
end

# True pole of date, eq. (23)
function _ori_ref_pd(a)
    ϵ = a.ϵ̄ + a.Δϵ
    return [sin(a.Δψ)*sin(ϵ),
            cos(a.Δψ)*sin(ϵ)*cos(a.ϵ̄) - cos(ϵ)*sin(a.ϵ̄),
            cos(a.Δψ)*sin(ϵ)*sin(a.ϵ̄) + cos(ϵ)*cos(a.ϵ̄)]
end

# Earth pole in the inertial frame, eq. (25)
_ori_ref_pE(a) = _ori_Rz(a.ζ)*_ori_Ry(-a.θ)*_ori_Rz(a.z)*_ori_Rx(-a.Φx)*_ori_Ry(-a.Φy)*_ori_ref_pd(a)

# Inertial -> true equator of date: nutation × linear corrections × precession
function _ori_ref_c2t(a)
    N = _ori_Rx(-(a.ϵ̄ + a.Δϵ))*_ori_Rz(-a.Δψ)*_ori_Rx(a.ϵ̄)
    C = _ori_Ry(a.Φy)*_ori_Rx(a.Φx)
    P = _ori_Rz(-a.z)*_ori_Ry(a.θ)*_ori_Rz(-a.ζ)
    return N*C*P
end

# Epochs (TDB days since J2000.0) within the DE430 time span (1550-2650)
const _ORI_EPOCHS = [-160_000.0, -54_787.5, 0.0, 7_300.0, 100_000.0, 230_000.0]

@testset "Earth orientation model (DE430/431)" begin

    @testset "Rotation matrices" begin
        for a in (0.3, -1.2, 2.5)
            @test Rx(a) ≈ _ori_Rx(a)
            @test Ry(a) ≈ _ori_Ry(a)
            @test Rz(a) ≈ _ori_Rz(a)
            for R in (Rx, Ry, Rz)
                @test R(a)*R(-a) ≈ I
                @test det(R(a)) ≈ 1
            end
        end
        # Frame (not vector) rotations
        @test Rz(π/2)*[1.0, 0.0, 0.0] ≈ [0.0, -1.0, 0.0] atol = 1e-15
        @test Rx(π/2)*[0.0, 1.0, 0.0] ≈ [0.0, 0.0, -1.0] atol = 1e-15
        @test Ry(π/2)*[1.0, 0.0, 0.0] ≈ [0.0, 0.0, 1.0] atol = 1e-15
        # In-place methods vs allocating methods (all Taylor coefficients)
        α = 0.7 + 0.01*Taylor1(10)
        for (R, R!) in ((Rx, PE.Rx!), (Ry, PE.Ry!), (Rz, PE.Rz!))
            res = [zero(α) for _ in 1:3, _ in 1:3]
            auxsR = [one(α), zero(α), zero(α), zero(α)]
            R!(res, α, auxsR)
            @test _ori_allclose(res, R(α); atol = 1e-14)
        end
    end

    @testset "Constants and polynomial angles" begin
        # Taylor coefficients wrt T (Julian centuries) reproduce the published polynomials
        T = 36525*Taylor1(5)
        as2deg = 1/3600
        for (f, c) in ((PE.Ω, [125+2/60+40.280/3600, -(1934+8/60+10.539/3600), 7.455/3600, 0.008/3600]),
                       (PE.ϵ̄, [84381.448, -46.815, -0.00059, 0.001813]*as2deg),
                       (PE.Zeta, [0.0, 2306.2181, 0.30188, 0.017998]*as2deg),
                       (PE.Theta, [0.0, 2004.3109, -0.42665, -0.041833]*as2deg),
                       (PE.zeta, [0.0, 2306.2181, 1.09468, 0.018203]*as2deg))
            x = f(T)
            @test x.coeffs[1:4] ≈ deg2rad.(c) rtol = 1e-14 atol = 1e-20
            @test iszero(x.coeffs[5:end])
            # Scalar method
            @test f(0.0) ≈ deg2rad(c[1]) rtol = 1e-15
        end
        # Linear corrections: values (arcseconds) and rates (arcseconds/year) of Table 10
        rad2as(x) = rad2deg(x)*3600
        @test rad2as(PE.phi_x(0.0)) ≈ 5.6754203322893470E-03 rtol = 1e-14
        @test rad2as(PE.phi_y(0.0)) ≈ -1.7022656914989530E-02 rtol = 1e-14
        @test rad2as(PE.phi_x(yr) - PE.phi_x(0.0)) ≈ 2.7689915574483550E-04 rtol = 1e-10
        @test rad2as(PE.phi_y(yr) - PE.phi_y(0.0)) ≈ -1.2118591216559240E-03 rtol = 1e-10
        # Nutation amplitudes (IAU 1980)
        for t in _ORI_EPOCHS
            Ω = PE.Ω(t)
            abs(sin(Ω)) > 0.1 && @test rad2as(PE.Delta_psi(t))/sin(Ω) ≈ -17.1996 rtol = 1e-13
            abs(cos(Ω)) > 0.1 && @test rad2as(PE.Delta_epsilon(t))/cos(Ω) ≈ 9.2025 rtol = 1e-13
        end
    end

    @testset "Reference implementation of eqs. (21)-(26)" begin
        as2rad = big(π)/648_000
        for t in _ORI_EPOCHS
            a = _ori_ref_angles(big(t), as2rad)
            M_ref = _ori_ref_c2t(a)
            pE_ref = _ori_ref_pE(a)
            # Consistency of the reference: eq. (25) is the third row of the matrix
            @test _ori_allclose(M_ref[3, :], pE_ref; atol = 1e-60)
            # Package (Float64) vs reference (BigFloat)
            M = c2t_jpl_de430(t)
            @test _ori_allclose(M, M_ref; atol = 1e-14)
            @test _ori_allclose(PE.pole_frame(t), pE_ref; atol = 1e-14)
            @test _ori_allclose(PE.pole_date(t), _ori_ref_pd(a); atol = 1e-14)
            @test _ori_allclose(PE.nutation_iau80(t)[3, :], _ori_ref_pd(a); atol = 1e-14)
        end
        # Taylor1 epochs: all coefficients (time scaled by 1000 days)
        for t0 in _ORI_EPOCHS
            t = t0 + 1000*Taylor1(8)
            M_ref = _ori_ref_c2t(_ori_ref_angles(t, π/648_000))
            @test _ori_allclose(c2t_jpl_de430(t), M_ref; atol = 1e-14)
        end
    end

    @testset "Matrix properties" begin
        for t in _ORI_EPOCHS
            M = c2t_jpl_de430(t)
            @test M*M' ≈ I atol = 1e-14
            @test det(M) ≈ 1 atol = 1e-14
            @test t2c_jpl_de430(t) == transpose(M)
            @test PE.pole_date(t) == PE.nutation_iau80(t)[3, :]
            @test PE.pole_frame(t) ≈ M[3, :] atol = 1e-14
            @test norm(PE.pole_frame(t)) ≈ 1 atol = 1e-14
            # True pole of date, eq. (23), from the package angles
            ϵ̄, Δψ, Δϵ = PE.ϵ̄(t), PE.Delta_psi(t), PE.Delta_epsilon(t)
            ϵ = ϵ̄ + Δϵ
            p_d = [sin(Δψ)*sin(ϵ), cos(Δψ)*sin(ϵ)*cos(ϵ̄) - cos(ϵ)*sin(ϵ̄),
                   cos(Δψ)*sin(ϵ)*sin(ϵ̄) + cos(ϵ)*cos(ϵ̄)]
            @test PE.pole_date(t) ≈ p_d atol = 1e-14
        end
    end

    @testset "Pole right ascension and declination" begin
        for t in _ORI_EPOCHS
            p = PE.pole_frame(t)
            α, δ = PE.pole_radec(t)
            @test α == PE.pole_ra(t)
            @test δ == PE.pole_dec(t)
            @test -π < α ≤ π
            @test [cos(δ)*cos(α), cos(δ)*sin(α), sin(δ)] ≈ p atol = 1e-14
            # Regression: wrong quadrant whenever p_x > 0 (most epochs after 2000)
            @test PE.earth_pole_rotation(t)[3, :] ≈ c2t_jpl_de430(t)[3, :] atol = 1e-14
        end
        t = 7_300.0 + Taylor1(5)
        @test _ori_allclose(PE.earth_pole_rotation(t)[3, :], c2t_jpl_de430(t)[3, :]; atol = 1e-14)
    end

    @testset "Taylor1 vs scalar methods" begin
        s = 0.2
        for t0 in _ORI_EPOCHS
            t = t0 + 1000*Taylor1(10)
            for f in (PE.Ω, PE.ϵ̄, PE.Zeta, PE.Theta, PE.zeta, PE.Delta_psi, PE.Delta_epsilon,
                      PE.phi_x, PE.phi_y)
                x = f(t)
                @test x isa Taylor1{Float64}
                @test constant_term(x) ≈ f(t0) rtol = 1e-12 atol = 1e-18
                @test x(s) ≈ f(t0 + 1000s) rtol = 1e-12 atol = 1e-18
            end
            for f in (PE.nutation_iau80, c2t_jpl_de430, t2c_jpl_de430)
                M = f(t)
                @test constant_term.(M) ≈ f(t0) atol = 1e-14
                @test evaluate.(M, s) ≈ f(t0 + 1000s) atol = 1e-14
            end
        end
        # Low Taylor orders
        for order in 1:2
            t = 7_300.0 + Taylor1(order)
            @test constant_term.(c2t_jpl_de430(t)) ≈ c2t_jpl_de430(7_300.0) atol = 1e-14
            for f in (PE.Ω, PE.ϵ̄, PE.Zeta, PE.Theta, PE.zeta, PE.Delta_psi, PE.Delta_epsilon)
                @test constant_term(f(t)) ≈ f(7_300.0) rtol = 1e-14 atol = 1e-20
            end
        end
    end

    @testset "In-place methods" begin
        t1 = -54_787.5 + 1000*Taylor1(8)
        t2 = 7_300.0 + 1000*Taylor1(8)
        rotatBuf = PE.allocate_c2t_jpl_de430(t1)
        M = Array{Taylor1{Float64}}(undef, 3, 3, 5)
        for j in 1:3, i in 1:3
            M[i, j, ea] = zero(t1)
        end
        # Reuse the same buffers for two different epochs (no state leaks)
        for t in (t1, t2)
            PE.t2c_jpl_de430!(M, ea, t, zero(t), rotatBuf)
            @test _ori_allclose(M[:, :, ea], permutedims(c2t_jpl_de430(t)); atol = 1e-15)
            @test _ori_allclose(M[:, :, ea], t2c_jpl_de430(t, zero(t)); atol = 1e-15)
        end
        # Without buffers, the allocating method is used
        @test PE.t2c_jpl_de430!(M, ea, t2, zero(t2), nothing) == t2c_jpl_de430(t2, zero(t2))
    end

end

@testset "Moon orientation model (DE430/431)" begin

    @testset "Mantle rotation matrix" begin
        # Euler angles of the lunar mantle, Table 7
        ϕ, θ, ψ = 0.00512830031411853500, 0.38239278420173690000, 1.29416700274878300000
        A = Rz(ψ)*Rx(θ)*Rz(ϕ)
        # Form used by ITM
        @test pole_rotation(ϕ - π/2, π/2 - θ, ψ) ≈ A atol = 1e-14
        # Element-wise form used by the dynamical models
        R = [cos(ϕ)*cos(ψ)-cos(θ)*sin(ϕ)*sin(ψ)    cos(ψ)*sin(ϕ)+cos(θ)*cos(ϕ)*sin(ψ)   sin(θ)*sin(ψ);
             -cos(θ)*cos(ψ)*sin(ϕ)-cos(ϕ)*sin(ψ)   cos(θ)*cos(ϕ)*cos(ψ)-sin(ϕ)*sin(ψ)   cos(ψ)*sin(θ);
             sin(θ)*sin(ϕ)                         -cos(ϕ)*sin(θ)                      cos(θ)]
        @test R ≈ A atol = 1e-14
        # Inertial -> mantle: r_I = Rz(-ϕ)Rx(-θ)Rz(-ψ) r_PA, eq. (10)
        @test A' ≈ Rz(-ϕ)*Rx(-θ)*Rz(-ψ) atol = 1e-14
    end

    @testset "Euler angle kinematics" begin
        # moon_omega is the inverse of eq. (14)
        ϕ = Taylor1([0.1, 2e-3, 1e-5, -3e-7], 10)
        θ = Taylor1([0.4, 1e-4, -2e-6, 1e-8], 10)
        ψ = Taylor1([1.3, 0.23, 1e-4, 2e-6], 10)
        ωx, ωy, ωz = PE.moon_omega(ϕ, θ, ψ)
        dϕ = (ωx*sin(ψ) + ωy*cos(ψ))/sin(θ)
        dθ = ωx*cos(ψ) - ωy*sin(ψ)
        dψ = ωz - dϕ*cos(θ)
        @test _ori_close(dϕ, PE.ordpres_differentiate(ϕ); atol = 1e-14)
        @test _ori_close(dθ, PE.ordpres_differentiate(θ); atol = 1e-14)
        @test _ori_close(dψ, PE.ordpres_differentiate(ψ); atol = 1e-14)
    end

    @testset "Moments of inertia" begin
        # Undistorted total moment of inertia, eqs. (37)-(38) and Table 11
        A_T, B_T, C_T = PE.A_T, PE.B_T, PE.C_T
        @test (C_T - A_T)/B_T ≈ 6.3102131934887270E-04 rtol = 1e-10
        @test (B_T - A_T)/C_T ≈ 2.2773171480091860E-04 rtol = 1e-10
        @test C_T - (A_T + B_T)/2 ≈ 2.0321568464952570E-04 rtol = 1e-10
        # Second zonal harmonic of the undistorted total moment of inertia, eq. (30)
        mR2 = μ[mo]*R_moon^2
        @test (ITM_und[3, 3] - (ITM_und[1, 1] + ITM_und[2, 2])/2)/mR2 ≈ 2.0321568464952570E-04 rtol = 1e-10
        @test isdiag(ITM_und)
        # Core moment of inertia, eq. (39) and Table 11
        @test isdiag(PE.I_c)
        @test PE.I_c[1, 1] == PE.I_c[2, 2]
        @test PE.I_c[3, 3]/mR2 ≈ 0.0007*C_T rtol = 1e-14
        @test (PE.I_c[3, 3] - PE.I_c[1, 1])/PE.I_c[3, 3] ≈ 2.4623904789198150E-04 rtol = 1e-10

        # Time-dependent part of the mantle moment of inertia, eq. (41)
        x, y, z = -8.0817735147818490E-04, -1.9946299854970130E-03, -1.0872626830706890E-03
        m1 = ITM1(x, y, z)
        @test issymmetric(m1)
        @test abs(tr(m1)) ≤ 1e-14*maximum(abs, m1)
        @test m1[1, 2] ≈ -(k_2M*μ[ea]*R_moon^5)/(x^2+y^2+z^2)^2.5*x*y
        ωm = [4.5737241859914330E-05, -2.1898617456729500E-06, 2.2994486018992250E-01]
        m2 = ITM2(ωm...)
        @test issymmetric(m2)
        @test abs(tr(m2)) ≤ 1e-14*(k_2M*R_moon^5/3)*sum(abs2, ωm)
        @test m2[1, 2] ≈ (k_2M*R_moon^5/3)*ωm[1]*ωm[2]

        # Mantle moment of inertia: Moon-Earth vector rotated to the mantle frame
        q = zeros(66)
        q[3ea-2:3ea] .= [0.12, -0.92, -0.40]
        q[3mo-2:3mo] .= q[3ea-2:3ea] .+ [x, y, z]
        eulang = [0.00512830031411853500, 0.38239278420173690000, 1.29416700274878300000]
        r = Rz(eulang[3])*Rx(eulang[2])*Rz(eulang[1])*[x, y, z]
        I_m = PE.ITM(q, eulang, ωm)
        @test issymmetric(I_m)
        @test I_m ≈ (ITM_und - PE.I_c) + ITM1(r...) + ITM2(ωm...) rtol = 1e-14
        # Time-dependent part alone (a wrong frame would change it at O(1))
        @test I_m - (ITM_und - PE.I_c) ≈ ITM1(r...) + ITM2(ωm...) rtol = 1e-6
    end

    @testset "Initial conditions (Table 7)" begin
        q0 = read_initial_conditions(joinpath(pkgdir(PlanetaryEphemeris), "data", "de430ic_1969Jun28.txt"))
        N = (length(q0) - 13) ÷ 6
        mantle = [0.00512830031411853500, 0.38239278420173690000, 1.29416700274878300000,
                  0.00004573724185991433, -0.00000218986174567295, 0.22994486018992250000]
        core = [-0.00241990927040684100, 0.41101946488652730000, -0.46309468558363680000,
                -0.00661836772247824400, -0.00107295445159005100, 0.22964879652299730000]
        @test q0[6N+1:6N+6] ≈ mantle rtol = 1e-14
        @test q0[6N+7:6N+12] ≈ core rtol = 1e-14
    end

end
