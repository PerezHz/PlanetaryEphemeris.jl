# This file is part of the PlanetaryEphemeris.jl package; MIT licensed

using PlanetaryEphemeris
using Dates
using Quadmath
using JLD2
using TaylorSeries
using Test
using Downloads

using PlanetaryEphemeris: freeparticle!
using SPICE: furnsh, spkgeo
using LinearAlgebra: norm

@testset "Propagation" begin

    # Total number of bodies
    local N = 11 + 343
    # Parameters
    local params = (N, J2000)
    # Number of years
    local nyears = 30.0
    # Order of Taylor expansions wrt time
    local order = 15
    # Absolute tolerance
    local abstol = 1E-12

    @testset "Initial conditions" begin
        # read_initial_conditions
        for date in (DateTime(1969, 6, 28), DateTime(2000, 1, 1, 12))
            filename = joinpath(pkgdir(PlanetaryEphemeris), "data",
                string("de430ic_", Dates.format(date, "yyyyud"), ".txt"))
            q0_64 = read_initial_conditions(filename)
            q0_128 = Float128.(q0_64)

            @test isa(q0_64, Vector{Float64})
            @test isa(q0_128, Vector{Float128})
            @test length(q0_64) == length(q0_128) == 6*N+13
            @test q0_64 == q0_128
        end
    end

    @testset "TaylorInterpolant" begin

        # Test zero TaylorInterpolant
        T = Float64
        U = TaylorN{T}
        @test iszero(zero(TaylorInterpolant{T, T, 2, Vector{T}, Matrix{Taylor1{T}}}))
        @test iszero(zero(TaylorInterpolant{T, U, 2, Vector{T}, Matrix{Taylor1{U}}}))
        @test iszero(zero(TaylorInterpolant{T, T, 2, SubArray{T, 1}, SubArray{Taylor1{T}, 2}}))
        @test iszero(zero(TaylorInterpolant{T, U, 2, SubArray{T, 1}, SubArray{Taylor1{U}, 2}}))

        # Test propagation
        tspan = (J2000, J2000 + nyears * yr)
        filename = joinpath(pkgdir(PlanetaryEphemeris), "data", "de430ic_2000Jan1.txt")
        q0 = read_initial_conditions(filename)
        PP = PlanetaryEphemerisProblem(freeparticle!, tspan, q0, params)
        sol = propagate(PP; order, abstol)

        @test sol isa TaylorInterpolant{T, T, 2}
        @test sol(sol.t0) == q0
        @test sol.t0 == 0.0
        @test sol.t[end] == nyears * yr
        @test length(sol.t) == size(sol.x, 1) + 1
        @test length(q0) == size(sol.x, 2)

        dq = TaylorSeries.set_variables("dq", order = 2, numvars = 2)
        tmid = sol.t0 + sol.t[2] / 2
        sol1N = TaylorInterpolant(sol.t0, sol.t, sol.x .+ Taylor1(dq[1], order))

        @test sol(tmid) isa Vector{T}
        @test sol(tmid + Taylor1(order)) isa Vector{Taylor1{T}}
        @test sol(tmid + dq[1] + dq[1] * dq[2]) isa Vector{TaylorN{T}}
        @test sol(tmid + Taylor1([dq[1],dq[1]*dq[2]], order)) isa Vector{Taylor1{TaylorN{T}}}
        @test sol1N(sol.t0)() == sol(sol.t0)
        @test sol1N(tmid)() == sol(tmid)

        # flipsign
        fsol = flipsign(sol)
        @test fsol.t0 == sol.t0
        @test fsol.t == -sol.t
        @test fsol.x == sol.x(-Taylor1(order))
        @test norm(fsol(-nyears*yr) - sol(nyears*yr)) < eps()

        # Test TaylorInterpolantSerialization
        @test JLD2.writeas(typeof(sol)) == PlanetaryEphemeris.TaylorInterpolantSerialization{Float64}
        jldsave("test.jld2"; sol)
        sol_file = JLD2.load("test.jld2", "sol")
        rm("test.jld2")
        @test sol_file == sol

        # Test TaylorInterpolantNSerialization
        sol1N = TaylorInterpolant(sol.t0, sol.t, sol.x .* Taylor1(one(dq[1]), 25))
        @test JLD2.writeas(typeof(sol1N)) == PlanetaryEphemeris.TaylorInterpolantNSerialization{Float64}
        jldsave("test.jld2"; sol1N)
        sol1N_file = JLD2.load("test.jld2", "sol1N")
        @test sol1N_file == sol1N
        rm("test.jld2")

    end

    @testset "Propagation: DE430 dynamical model" begin

        @show Threads.nthreads()

        # Planetary ephemeris problem
        tspan = (J2000, J2000 + nyears * yr)ß
        filename = joinpath(pkgdir(PlanetaryEphemeris), "data", "de430ic_2000Jan1.txt")
        q0 = read_initial_conditions(filename)
        PP = PlanetaryEphemerisProblem(DE430!, tspan, q0, params)
        # Test propagation
        propagate(PP; maxsteps = 1, order, abstol)
        sol = propagate(PP; maxsteps = 100, order, abstol)

        # Indices of bodies to be saved
        bodyind = 1:(11+16)
        # Save solution
        filename = selecteph2jld2(sol, bodyind, nyears)
        # Recovered solution
        recovered_sol = JLD2.load(filename, "ss16ast_eph")

        @test selecteph(sol, bodyind, euler = true, ttmtdb = true) == recovered_sol

        # Test selecteph
        t0 = sol.t0 + sol.t[end]/3
        tf = sol.t0 + 2*sol.t[end]/3
        idxs = vcat(nbodyind(N, [su, ea, mo]), 6N+1:6N+13)
        i_0 = searchsortedlast(sol.t, t0)
        i_f = searchsortedfirst(sol.t, tf)

        subsol = selecteph(sol, [su, ea, mo], t0, tf; euler = true, ttmtdb = true)

        @test subsol.t0 == sol.t0
        @test subsol.t0 + subsol.t[1] ≤ t0
        @test subsol.t0 + subsol.t[end] ≥ tf
        @test size(subsol.x) == (i_f - i_0, length(idxs))
        @test subsol.x == sol.x[i_0:i_f-1, idxs]
        @test subsol(t0) == sol(t0)[idxs]
        @test subsol(tf) == sol(tf)[idxs]

        # Kernels URLs
        LSK = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/naif0012.tls"
        TTmTDBK = "https://ssd.jpl.nasa.gov/ftp/eph/planets/bsp/TTmTDB.de430.19feb2015.bsp"
        SPK = "https://ssd.jpl.nasa.gov/ftp/eph/planets/bsp/de430_1850-2150.bsp"

        # Download kernels
        Downloads.download(LSK, "naif0012.tls")
        Downloads.download(SPK, "de430_1850-2150.bsp")
        Downloads.download(TTmTDBK, "TTmTDB.de430.19feb2015.bsp")

        # Load kernels
        furnsh("naif0012.tls", "de430_1850-2150.bsp", "TTmTDB.de430.19feb2015.bsp")

        ttmtdb_pe = TaylorInterpolant(sol.t0, sol.t, sol.x[:, 6N+13]) # TT-TDB
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

        tv = range(sol.t0, sol.t[end], 10)
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

        # Remove files
        rm.((filename, "naif0012.tls", "de430_1850-2150.bsp", "TTmTDB.de430.19feb2015.bsp"))

        # TO DO: test propagation with Float128

    end

end
