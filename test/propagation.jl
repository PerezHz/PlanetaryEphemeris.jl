# This file is part of the PlanetaryEphemeris.jl package; MIT licensed

using PlanetaryEphemeris
using Dates
using Quadmath
using TaylorIntegration
using JLD2
using Test

using PlanetaryEphemeris: initialcond, ssic_1969, ssic_2000, astic_1969, astic_2000

@testset "Initial conditions" begin

    # Number of asteroids
    local N_ast = 343
    # Total number of bodies (solar system + asteroids)
    local N = 11 + N_ast
    # Check initial conditions global variables
    @test size(ssic_1969) == (11, 6)
    @test size(ssic_2000) == (11, 6)

    @test isa(ssic_1969, Matrix{Float64})
    @test isa(ssic_2000, Matrix{Float64})

    @test size(astic_1969) == (N_ast, 6)
    @test size(astic_2000) == (N_ast, 6)

    @test isa(astic_1969, Matrix{Float64})
    @test isa(astic_2000, Matrix{Float64})

    # Number of bodies
    local N = 11 + 343
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

using Downloads
using SPICE: furnsh, spkgeo
using PlanetaryEphemeris: order, abstol
using LinearAlgebra: norm

@testset "Propagation" begin

    @show Threads.nthreads()

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

    # Float64

    # Test integration
    sol64 = propagate(10, jd0, nyears, dense; dynamics = dynamics, order = order, abstol = abstol)
    # Save solution
    filename = selecteph2jld2(sol64, bodyind, nyears)
    # Recovered solution
    recovered_sol64 = JLD2.load(filename, "ss16ast_eph")

    @test selecteph(sol64, bodyind, euler = true, ttmtdb = true) == recovered_sol64

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

    etv = range(sol64.t0, sol64.t[end], 5)
    for et in eachindex(etv)
        @test abs(ttmtdb_pe(et) - ttmtdb_jpl(et)) < 1e-18
        @test norm(posvel_jpl_su(et) - posvel_pe_su(et), Inf) < 1e-17
        @test norm(posvel_jpl_ea(et) - posvel_pe_ea(et), Inf) < 1e-14
        @test norm(posvel_jpl_mo(et) - posvel_pe_mo(et), Inf) < 1e-14
        @test norm(posvel_jpl_ma(et) - posvel_pe_ma(et), Inf) < 1e-14
        @test norm(posvel_jpl_ju(et) - posvel_pe_ju(et), Inf) < 1e-14
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
