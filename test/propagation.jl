# This file is part of the PlanetaryEphemeris.jl package; MIT licensed

using PlanetaryEphemeris
using Dates
using Quadmath
using JLD2
using TaylorSeries
using Test

using PlanetaryEphemeris: initialcond, ssic_1969, ssic_2000, astic_1969, astic_2000

using Downloads
using SPICE: furnsh, spkgeo
using PlanetaryEphemeris: order, abstol
using LinearAlgebra: norm

@testset "Propagation" begin

    # Number of asteroids
    local N_ast = 343
    # Total number of bodies (solar system + asteroids)
    local N = 11 + N_ast
    # Indices of bodies to be saved
    local bodyind = 1:(11+16)
    # Starting time of integration
    local jd0 = datetime2julian(DateTime(2000,1,1,12))
    # Number of years
    local nyears = 2031.0 - year(julian2datetime(jd0))
    # Dense output
    local dense = Val(true)
    # Dynamical function
    local dynamics = DE430!

    @testset "Initial conditions" begin
        # Check initial conditions global variables
        @test size(ssic_1969) == (11, 6)
        @test size(ssic_2000) == (11, 6)

        @test isa(ssic_1969, Matrix{Float64})
        @test isa(ssic_2000, Matrix{Float64})

        @test size(astic_1969) == (N_ast, 6)
        @test size(astic_2000) == (N_ast, 6)

        @test isa(astic_1969, Matrix{Float64})
        @test isa(astic_2000, Matrix{Float64})

        # Check output of initialcond
        for jd0 in datetime2julian.( ( DateTime(1969,6,28,0,0,0), DateTime(2000,1,1,12) ) )

            q0_64 = initialcond(N, jd0)
            q0_128 = initialcond(N, convert(Float128, jd0) )

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
        # Test integration
        sol = propagate(5, jd0, nyears, dense; dynamics=PlanetaryEphemeris.freeparticle!, order, abstol)
        @test sol isa TaylorInterpolant{Float64, Float64, 2}
        q0 = initialcond(N, jd0)
        @test sol(sol.t0) == q0
        @test sol.t0 == 0.0
        @test length(sol.t) == size(sol.x, 1) + 1
        @test length(q0) == size(sol.x, 2)
        dq = TaylorSeries.set_variables("dq", order=2, numvars=2)
        tmid = sol.t0 + sol.t[2]/2
        @test sol(tmid) isa Vector{Float64}
        @test sol(tmid + Taylor1(order)) isa Vector{Taylor1{Float64}}
        @test sol(tmid + dq[1] + dq[1]*dq[2]) isa Vector{TaylorN{Float64}}
        @test sol(tmid + Taylor1([dq[1],dq[1]*dq[2]], order)) isa Vector{Taylor1{TaylorN{Float64}}}
        sol1N = TaylorInterpolant(sol.t0, sol.t, sol.x .+ Taylor1(dq[1], 25))
        @test sol1N(sol.t0)() == sol(sol.t0)
        @test sol1N(tmid)() == sol(tmid)
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

        # Float64

        # Test integration
        sol64 = propagate(100, jd0, nyears, dense; dynamics, order, abstol)
        # Save solution
        filename = selecteph2jld2(sol64, bodyind, nyears)
        # Recovered solution
        recovered_sol64 = JLD2.load(filename, "ss16ast_eph")

        @test selecteph(sol64, bodyind, euler = true, ttmtdb = true) == recovered_sol64

        # Test selecteph
        t0 = sol64.t0 + sol64.t[end]/3
        tf = sol64.t0 + 2*sol64.t[end]/3
        idxs = vcat(nbodyind(N, [su, ea, mo]), 6N+1:6N+13)
        i_0 = searchsortedlast(sol64.t, t0)
        i_f = searchsortedfirst(sol64.t, tf)

        subsol = selecteph(sol64, [su, ea, mo], t0, tf; euler = true, ttmtdb = true)

        @test subsol.t0 == sol64.t0
        @test subsol.t0 + subsol.t[1] ≤ t0
        @test subsol.t0 + subsol.t[end] ≥ tf
        @test size(subsol.x) == (i_f - i_0, length(idxs))
        @test subsol.x == sol64.x[i_0:i_f-1, idxs]
        @test subsol(t0) == sol64(t0)[idxs]
        @test subsol(tf) == sol64(tf)[idxs]

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

        ttmtdb_pe = TaylorInterpolant(sol64.t0, sol64.t, sol64.x[:, 6N+13]) # TT-TDB
        posvel_pe_su = selecteph(sol64,su) # Sun
        posvel_pe_ea = selecteph(sol64,ea) # Earth
        posvel_pe_mo = selecteph(sol64,mo) # Moon
        posvel_pe_ma = selecteph(sol64,6) # Mars
        posvel_pe_ju = selecteph(sol64,7) # Jupiter

        ttmtdb_jpl(et) = spkgeo(1000000001, et, "J2000", 1000000000)[1][1] # TT-TDB
        posvel_jpl_su(et) = kmsec2auday(spkgeo(10, et, "J2000", 0)[1]) # Sun
        posvel_jpl_ea(et) = kmsec2auday(spkgeo(399, et, "J2000", 0)[1]) # Earth
        posvel_jpl_mo(et) = kmsec2auday(spkgeo(301, et, "J2000", 0)[1]) # Moon
        posvel_jpl_ma(et) = kmsec2auday(spkgeo(4, et, "J2000", 0)[1]) # Mars
        posvel_jpl_ju(et) = kmsec2auday(spkgeo(5, et, "J2000", 0)[1]) # Jupiter

        tv = range(sol64.t0, sol64.t[end], 10)
        for t in tv
            et = t * daysec
            @show t, et
            @show abs(ttmtdb_jpl(et) - ttmtdb_pe(t))
            @show norm(posvel_jpl_su(et) - posvel_pe_su(t), Inf)
            @show norm(posvel_jpl_ea(et) - posvel_pe_ea(t), Inf)
            @show norm(posvel_jpl_mo(et) - posvel_pe_mo(t), Inf)
            @show norm(posvel_jpl_ma(et) - posvel_pe_ma(t), Inf)
            @show norm(posvel_jpl_ju(et) - posvel_pe_ju(t), Inf)

            @test abs(ttmtdb_jpl(et) - ttmtdb_pe(t)) < 1e-12
            @test norm(posvel_jpl_su(et) - posvel_pe_su(t), Inf) < 1e-14
            @test norm(posvel_jpl_ea(et) - posvel_pe_ea(t), Inf) < 1e-11
            @test norm(posvel_jpl_mo(et) - posvel_pe_mo(t), Inf) < 1e-11
            @test norm(posvel_jpl_ma(et) - posvel_pe_ma(t), Inf) < 1e-12
            @test norm(posvel_jpl_ju(et) - posvel_pe_ju(t), Inf) < 1e-13
        end

        # Remove files
        rm.((filename, "naif0012.tls", "de430_1850-2150.bsp", "TTmTDB.de430.19feb2015.bsp"))

        # Float 128
        #=
        # Test integration
        sol128 = propagate(1, Float128(jd0), nyears, dense; dynamics = dynamics, order = order, abstol = abstol)
        # Save solution
        filename = selecteph2jld2(sol128, bodyind, nyears)
        # Recovered solution
        recovered_sol128 = JLD2.load(filename, "ss16ast_eph")
        # Remove file
        rm(filename)

        @test selecteph(sol128, bodyind, euler = true, ttmtdb = true) == recovered_sol128
        =#
    end

end
