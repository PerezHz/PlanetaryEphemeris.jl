# This file is part of the PlanetaryEphemeris.jl package; MIT licensed

using PlanetaryEphemeris
using Dates
using Quadmath
using TaylorIntegration
using TaylorSeries
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

using PlanetaryEphemeris: jd0, nyears, dense, dynamics, order, abstol

@testset "Propagation" begin

    @show Threads.nthreads()

    # Number of asteroids
    local N_ast = 343
    # Total number of bodies (solar system + asteroids)
    local N = 11 + N_ast
    # Indices of bodies to be saved
    local bodyind = 1:(11+16)
    # Years to be integrated
    local nyears = 31

    # Float64

    # Test integration
    sol64 = propagate(1, jd0, nyears, dense; dynamics = dynamics, order = order, abstol = abstol)
    # Save solution
    filename = selecteph2jld2(sol64, bodyind, nyears)
    # Recovered solution
    recovered_sol64 = JLD2.load(filename, "ss16ast_eph")
    # Remove file
    rm(filename)

    @test selecteph(sol64, bodyind, euler = true, ttmtdb = true) == recovered_sol64

    # Float 128

    # Test integration
    sol128 = propagate(1, Float128(jd0), nyears, dense; dynamics = dynamics, order = order, abstol = abstol)
    # Save solution
    filename = selecteph2jld2(sol128, bodyind, nyears)
    # Recovered solution
    recovered_sol128 = JLD2.load(filename, "ss16ast_eph")
    # Remove file
    rm(filename)

    @test selecteph(sol128, bodyind, euler = true, ttmtdb = true) == recovered_sol128
end
