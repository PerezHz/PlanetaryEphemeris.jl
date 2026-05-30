# This file is part of the PlanetaryEphemeris.jl package; MIT licensed

using PlanetaryEphemeris
using TaylorIntegration
using Dates
using Quadmath
using JLD2
using TaylorSeries
using Test

using PlanetaryEphemeris: freeparticle!, sun_posvel_pN, loadeph
using SPICE: furnsh, spkgeo
using LinearAlgebra: norm

const PKG_DATA = joinpath(pkgdir(PlanetaryEphemeris), "data")
const TEST_DATA = joinpath(pkgdir(PlanetaryEphemeris), "test", "data")

# Total number of bodies
const N = 11 + 343
# Number of years
const nyears = 30.0
# Order of Taylor expansions wrt time
const _order = 15
# Absolute tolerance
const abstol = 1E-12

# Maximum relative error
mre(x, y, z) = maximum(@. abs(x - y) / z)

@testset "Osculating elements" begin
    # Pérez-Hernández & Benet (2022) Apophis OR7 orbit
    # See Tables 2 and 3 of the Supplementary Information in
    # https://doi.org/10.1038/s43247-021-00337-x

    # Reference epoch [JDTDB]
    jd0 = 2459200.5
    # Gravitational parameters [au^3/day^2]
    μ_S, μ_A = PE.GMS, 0.0

    # Cartesian state vector [au, au/day]
    rv = [-0.17380033054708824, 0.99450381407683, -0.057081276307106465,
          -0.016259012007155557, -0.00011129354506606208, -0.0003802783732738447]
    drv = [7.12E−9, 1.94E−9, 5.41E−9, 4.79E−11, 6.72E−11, 1.38E−10]

    @test rv ≈ kmsec2auday(auday2kmsec(rv))
    @test drv ≈ kmsec2auday(auday2kmsec(drv))

    # Keplerian elements
    e, de = 0.19150886716, 1.60E-9
    q, dq = 0.74585305033, 1.54E−9                           # au
    a, da = q / (1 - e), hypot(dq/(1-e), q*de/(1-e)^2)       # au
    i, di = 3.336773201, 1.74E−7                             # deg
    Ω, dΩ = 204.04199116, 8.81E−6                            # deg
    ω, dω = 126.65396094, 9.37E−6                            # deg
    tp, dtp = 2459101.04092537, 1.17E−6                      # JDTDB
    M = rad2deg(sqrt(μ_S / a^3)) * (jd0 - tp)                # deg
    dM = hypot(-3*M*da/(2a), -M*dtp/(jd0 - tp))

    _a_ = semimajoraxis(rv..., μ_S, μ_A)
    _e_ = eccentricity(rv..., μ_S, μ_A)
    _i_ = rad2deg(inclination(rv...))
    _Ω_ = rad2deg(longascnode(rv...))
    _ω_ = rad2deg(argperi(rv..., μ_S, μ_A))
    _tp_ = timeperipass(jd0, rv..., μ_S, μ_A)
    _M_ = rad2deg(meananomaly(rv..., μ_S, μ_A))

    @test mre(_a_, a, da) < 0.05
    @test mre(_e_, e, de) < 0.05
    @test mre(_i_, i, di) < 0.05
    @test mre(_Ω_, Ω, dΩ) < 0.05
    @test mre(_ω_, ω, dω) < 0.05
    @test mre(_tp_, tp, dtp) < 0.05
    @test mre(_M_, M, dM) < 0.05

    E = eccentricanomaly(e, M)
    ν = trueanomaly(e, E)
    @test E ≈ deg2rad(eccentricanomaly(e, M, Val(true)))
    @test ν ≈ deg2rad(trueanomaly(e, E, Val(true)))

    e, M = 1 + e, 0.0
    H = deg2rad(eccentricanomaly(e, M, Val(false)))
    θ = deg2rad(trueanomaly(e, H, Val(false)))
    @test a * (1 - e^2) / (1 + e * cos(θ)) ≈ a * (1 - e * cosh(H))

end

@testset "Propagation" begin

    @testset "Initial conditions" begin
        # read_initial_conditions
        for date in (DateTime(1969, 6, 28), DateTime(2000, 1, 1, 12))
            filename = joinpath(PKG_DATA, string("de430ic_", Dates.format(date, "yyyyud"), ".txt"))
            q0_64 = read_initial_conditions(filename)
            q0_128 = Float128.(q0_64)

            @test isa(q0_64, Vector{Float64})
            @test isa(q0_128, Vector{Float128})
            @test length(q0_64) == length(q0_128) == 6*N+13
            @test q0_64 == q0_128
        end
    end

    @testset "Earth orientation model" begin
        using PlanetaryEphemeris: allocate_c2t_jpl_de430, t2c_jpl_de430!

        t = Taylor1(_order)
        xT = Taylor1(rand(_order+1))
        xTN = Taylor1([rand() * TaylorN(1) + rand() * TaylorN(2) for _ in 0:_order])
        xT1 = Taylor1([Taylor1(rand(7)) for _ in 0:_order])

        M0 = t2c_jpl_de430(t)
        MT = t2c_jpl_de430(t, zero(xT))
        MTN = t2c_jpl_de430(t, zero(xTN))
        MT1 = t2c_jpl_de430(t, zero(xT1))

        @test MT == M0
        @test constant_term.(MTN) == constant_term.(M0)
        @test constant_term.(MT1) == constant_term.(M0)

        rotatBuf = allocate_c2t_jpl_de430(t)
        MMT = Array{typeof(xT)}(undef, 3, 3, 10)
        MMTN = Array{typeof(xTN)}(undef, 3, 3, 10)
        MMT1 = Array{typeof(xT1)}(undef, 3, 3, 10)

        @. MMT[:, :, ea] = zero(xT)
        @. MMTN[:, :, ea] = zero(xTN)
        @. MMT1[:, :, ea] = zero(xT1)

        t2c_jpl_de430!(MMT, ea, t, zero(xT), rotatBuf)
        t2c_jpl_de430!(MMTN, ea, t, zero(xTN), rotatBuf)
        t2c_jpl_de430!(MMT1, ea, t, zero(xT1), rotatBuf)

        @test all(@. isapprox(constant_term(MMT[:, :, ea]), constant_term(M0)))
        @test all(@. isapprox(constant_term(MMTN[:, :, ea]), constant_term(M0)))
        @test all(@. isapprox(constant_term(MMT1[:, :, ea]), constant_term(M0)))
    end

    @testset "TaylorSolution" begin

        # Test zero TaylorSolution
        T = Float64
        U = TaylorN{T}
        @test iszero(zero(TaylorSolution{T, T, 2, Vector{T}, Matrix{T}, Matrix{Taylor1{T}}, Nothing, Nothing, Nothing}))
        @test iszero(zero(TaylorSolution{T, U, 2, Vector{T}, Matrix{U}, Matrix{Taylor1{U}}, Nothing, Nothing, Nothing}))

        # Test propagation
        tspan = (J2000, J2000 + nyears * yr)
        filename = joinpath(PKG_DATA, "de430ic_2000Jan1.txt")
        q0 = read_initial_conditions(filename)
        params = (N, J2000)
        PP = PlanetaryEphemerisProblem(freeparticle!, tspan, q0, params)
        sol = propagate(PP; order = _order, abstol)

        @test isa(string(PP), String)
        @test isa(string(sol), String)
        @test isa(sol, TaylorSolution{T, T, 2})
        @test isa(PP, PlanetaryEphemerisProblem{typeof(freeparticle!), T, Tuple{Int64, T}})
        @test sol.t[1] == 0.0
        @test sol.t[end] == nyears * yr
        @test length(sol.t) == size(sol.p, 1) + 1
        @test length(q0) == size(sol.p, 2)
        @test numberofbodies(size(sol.p, 2)) == numberofbodies(sol.p[1, :]) ==
              numberofbodies(sol.p) == numberofbodies(sol) == N

        # Evaluation
        t = Taylor1(_order)
        @test sol(sol.t[1]) == q0
        @test constant_term(sol(t)) == q0
        @test all(iszero, sol(1, 1, sol.t[1]))
        @test sol(1, sol.t[1]) == q0[nbodyind(N, 1)]

        # Convert
        sol128 = convert(Float128, sol)
        @test sol128.t == Float128.(sol.t)
        @test all(@. constant_term(sol128.p) == Float128(constant_term(sol.p)))

        # Reverse
        solrev = reverse(sol)
        @test solrev.t[1] == sol.t[end]
        @test solrev.t == reverse(sol.t)
        @test solrev(solrev.t[end]) ≈ sol(sol.t[1])

        dq = TaylorSeries.variables!("dq", order = 2, numvars = 2)
        tmid = sol.t[2] / 2
        sol1N = TaylorSolution(sol.t, sol.p .+ Taylor1(dq[1], _order))

        @test sol(tmid) isa Vector{T}
        @test sol(tmid + Taylor1(_order)) isa Vector{Taylor1{T}}
        @test sol(tmid + dq[1] + dq[1] * dq[2]) isa Vector{TaylorN{T}}
        @test sol(tmid + Taylor1([dq[1],dq[1]*dq[2]], _order)) isa Vector{Taylor1{TaylorN{T}}}
        @test sol1N(sol.t[1])() == sol(sol.t[1])
        @test sol1N(tmid)() == sol(tmid)

        # flipsign
        fsol = flipsign(sol)
        @test fsol.t[1] == sol.t[1]
        @test fsol.t == -sol.t
        @test fsol.p == sol.p(-Taylor1(_order))
        @test norm(fsol(-nyears*yr) - sol(nyears*yr)) < eps()

        # Test TaylorSolutionSerialization
        @test JLD2.writeas(typeof(sol)) == TaylorIntegration.TaylorSolutionSerialization{Float64, 2}
        jldsave("test.jld2"; sol)
        sol_file = JLD2.load("test.jld2", "sol")
        rm("test.jld2")
        @test sol_file == sol

        # Test TaylorSolutionNSerialization
        sol1N = TaylorSolution(sol.t, sol.p .* Taylor1(one(dq[1]), 25))
        @test JLD2.writeas(typeof(sol1N)) == TaylorIntegration.TaylorSolutionNSerialization{Float64, 2}
        jldsave("test.jld2"; sol1N)
        sol1N_file = JLD2.load("test.jld2", "sol1N")
        @test sol1N_file == sol1N
        rm("test.jld2")

    end

    @testset "Propagation: DE430 dynamical model" begin

        @show Threads.nthreads()

        # Planetary ephemeris problem
        tspan = (J2000, J2000 + nyears * yr)
        filename = joinpath(PKG_DATA, "de430ic_2000Jan1.txt")
        q0 = read_initial_conditions(filename)
        params = DE430Params(J2000, q0, _order)
        PP = PlanetaryEphemerisProblem(DE430!, tspan, q0, params)
        # Test propagation
        @time propagate(PP; maxsteps = 1, order = _order, abstol)
        @time sol = propagate(PP; maxsteps = 100, order = _order, abstol)

        # Solar system barycenter
        rvec_ssb, vvec_ssb, μ_star_SSB = ssb_posvel_pN(PE.μ, q0)
        rvec_sun, vvec_sun = sun_posvel_pN(PE.μ, q0)

        _rvec_ssb_ = [
            sum(PE.μ[i] * q0[j] for (i, j) in enumerate(1:3:3N)),
            sum(PE.μ[i] * q0[j] for (i, j) in enumerate(2:3:3N)),
            sum(PE.μ[i] * q0[j] for (i, j) in enumerate(3:3:3N))
        ]
        _vvec_ssb_ = [
            sum(PE.μ[i] * q0[j] for (i, j) in enumerate(3N+1:3:6N)),
            sum(PE.μ[i] * q0[j] for (i, j) in enumerate(3N+2:3:6N)),
            sum(PE.μ[i] * q0[j] for (i, j) in enumerate(3N+3:3:6N))
        ]

        @test norm(rvec_ssb) < norm(_rvec_ssb_)
        @test norm(vvec_ssb) < norm(_vvec_ssb_)
        @test μ_star_SSB ≈ sum(PE.μ)
        @test q0[sundofs] ≈ vcat(rvec_sun, vvec_sun)

        # Save results
        bodyind = 1:(11+16)
        filename = selecteph2jld2(sol, bodyind, nyears)
        recovered_sol = JLD2.load(filename, "ss16ast_eph")

        @test selecteph(sol, bodyind, euler = true, ttmtdb = true) == recovered_sol

        acceph, poteph = loadeph(recovered_sol, PE.μ)

        @test size(acceph.p, 1) == size(poteph.p, 1) == size(recovered_sol.p, 1)
        @test size(acceph.p, 2) == (size(recovered_sol.p, 2) - 13) ÷ 2
        @test size(poteph.p, 2) == (size(recovered_sol.p, 2) - 13) ÷ 6

        @test isnothing(save2jld2andcheck(filename, Dict(
           "sseph" => recovered_sol,
           "acceph" => acceph,
           "poteph" => poteph
        )))
        rm(filename)

        # Test selecteph
        t0 = sol.t[end]/3
        tf = 2*sol.t[end]/3
        idxs = vcat(nbodyind(N, [su, ea, mo]), 6N+1:6N+13)
        i_0 = searchsortedlast(sol.t, t0)
        i_f = searchsortedfirst(sol.t, tf)

        subsol = selecteph(sol, [su, ea, mo], t0, tf; euler = true, ttmtdb = true)

        @test subsol.t[1] ≤ t0
        @test subsol.t[end] ≥ tf
        @test size(subsol.p) == (i_f - i_0, length(idxs))
        @test subsol.p == sol.p[i_0:i_f-1, idxs]
        @test subsol(t0) == sol(t0)[idxs]
        @test subsol(tf) == sol(tf)[idxs]

        # Load kernels
        furnsh(
            joinpath(TEST_DATA, "naif0012.tls"),
            joinpath(TEST_DATA, "de430_2000-2002.bsp")
        )

        ttmtdb_pe = TaylorSolution(sol.t, sol.p[:, 6N+13]) # TT-TDB
        posvel_pe_su = selecteph(sol, su) # Sun
        posvel_pe_ea = selecteph(sol, ea) # Earth
        posvel_pe_mo = selecteph(sol, mo) # Moon
        posvel_pe_ma = selecteph(sol, 6) # Mars
        posvel_pe_ju = selecteph(sol, 7) # Jupiter

        ttmtdb_jpl(et) = spkgeo(1000000001, et, "J2000", 1000000000)[1][1] # TT-TDB
        posvel_jpl_su(et) = kmsec2auday(spkgeo(10, et, "J2000", 0)[1]) # Sun
        posvel_jpl_ea(et) = kmsec2auday(spkgeo(399, et, "J2000", 0)[1]) # Earth
        posvel_jpl_mo(et) = kmsec2auday(spkgeo(301, et, "J2000", 0)[1]) # Moon
        posvel_jpl_ma(et) = kmsec2auday(spkgeo(4, et, "J2000", 0)[1]) # Mars
        posvel_jpl_ju(et) = kmsec2auday(spkgeo(5, et, "J2000", 0)[1]) # Jupiter

        tv = range(sol.t[1], sol.t[end], 10)
        for t in tv
            et = t * daysec
            @show t, et
            @show abs(ttmtdb_jpl(et) - ttmtdb_pe(t))
            @show norm(posvel_jpl_su(et) - posvel_pe_su(t), Inf)
            @show norm(posvel_jpl_ea(et) - posvel_pe_ea(t), Inf)
            @show norm(posvel_jpl_mo(et) - posvel_pe_mo(t), Inf)
            @show norm(posvel_jpl_ma(et) - posvel_pe_ma(t), Inf)
            @show norm(posvel_jpl_ju(et) - posvel_pe_ju(t), Inf)

            @test abs(ttmtdb_jpl(et) - ttmtdb_pe(t)) < 1E-12
            @test norm(posvel_jpl_su(et) - posvel_pe_su(t), Inf) < 1E-14
            @test norm(posvel_jpl_ea(et) - posvel_pe_ea(t), Inf) < 1E-11
            @test norm(posvel_jpl_mo(et) - posvel_pe_mo(t), Inf) < 4E-11
            @test norm(posvel_jpl_ma(et) - posvel_pe_ma(t), Inf) < 1E-12
            @test norm(posvel_jpl_ju(et) - posvel_pe_ju(t), Inf) < 1E-13
        end

        # TO DO: test propagation with Float128

    end

end
