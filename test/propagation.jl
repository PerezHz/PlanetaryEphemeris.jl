# This file is part of the PlanetaryEphemeris.jl package; MIT licensed

using PlanetaryEphemeris
using Dates 
using Quadmath
using TaylorIntegration
using TaylorSeries 
using JLD2
using Test

using PlanetaryEphemeris: initialcond

@testset "Initial conditions" begin  

    # Special method of julian2datetime for Float128
    local date = DateTime( rand(1950:2050), rand(1:12), rand(1:28) )
    local jd0 = datetime2julian(date)
    local jd1 = convert(Float128, jd0)

    @test jd0 == jd1 
    @test julian2datetime(jd0) == julian2datetime(jd1)

    # initialcond
    local N = 11 + 343

    for jd0 in datetime2julian.( ( DateTime(1969,6,28,0,0,0), DateTime(2000,1,1,12) ) )
        
        q0_64 = initialcond(N, jd0)
        q0_128 = initialcond(N, convert(Float128, jd0) )

        @test isa(q0_64, Vector{Float64})
        @test isa(q0_128, Vector{Float128})
        @test q0_64 == q0_128

    end 

end 

@testset "Interpolation" begin
    
    # Kepler problem
    local order = 28

    @taylorize function kepler1!(dq, q, p, t)
        local μ = -1.0
        r = sqrt(q[1]^2+q[2]^2)
        r_p3d2 = r^3

        dq[1] = q[3]
        dq[2] = q[4]
        dq[3] = μ * q[1] / r_p3d2
        dq[4] = μ * q[2] / r_p3d2

        return nothing
    end

    # Float64 
    local abstol_64 = 1.0e-20
    local t0_64 = 0.0
    local tf_64 = 2π*100.0

    q0_64 = [0.2, 0.0, 0.0, 3.0]

    tv_f64, xv_f64, polynV_f64 = taylorinteg(kepler1!, q0_64, t0_64, tf_64, order, abstol_64, Val(true), maxsteps = 500000, 
                                             parse_eqs = false)
    interp_f64 = TaylorInterpolant(tv_f64[1], tv_f64[:], polynV_f64[2:end, :])

    @test isa(interp_f64, TaylorInterpolant{Float64, Float64, 2})

    tv_t64, xv_t64, polynV_t64 = taylorinteg(kepler1!, q0_64, t0_64, tf_64, order, abstol_64, Val(true), maxsteps = 500000, 
                                             parse_eqs = true)
    interp_t64 = TaylorInterpolant(tv_t64[1], tv_t64[:], polynV_t64[2:end, :])

    @test isa(interp_t64, TaylorInterpolant{Float64, Float64, 2})

    @test interp_f64 == interp_t64

    # Float128 
    local abstol_128 = Float128(1.0e-20)
    local t0_128 = Float128(0.0)
    local tf_128 = Float128(2π*100.0)

    q0_128 = Float128.([0.2, 0.0, 0.0, 3.0])

    tv_f128, xv_f128, polynV_f128 = taylorinteg(kepler1!, q0_128, t0_128, tf_128, order, abstol_128, Val(true), maxsteps = 500000, 
                                                parse_eqs = false)
    interp_f128 = TaylorInterpolant(tv_f128[1], tv_f128[:], polynV_f128[2:end, :])

    @test isa(interp_f128, TaylorInterpolant{Float128, Float128, 2})

    tv_t128, xv_t128, polynV_t128 = taylorinteg(kepler1!, q0_128, t0_128, tf_128, order, abstol_128, Val(true), maxsteps = 500000, 
                                                parse_eqs = true)
    interp_t128 = TaylorInterpolant(tv_t128[1], tv_t128[:], polynV_t128[2:end, :])

    @test isa(interp_t128, TaylorInterpolant{Float128, Float128, 2})

    @test interp_f128 == interp_t128

    # Crossed tests
    @test typeof(convert(Float128, interp_t64)) == typeof(interp_t128)
    @test typeof(convert(Float128, interp_f64)) == typeof(interp_f128)

    @test typeof(convert(Float64, interp_t128)) == typeof(interp_t64)
    @test typeof(convert(Float64, interp_f128)) == typeof(interp_f64)

end 

@testset "Read/write .jld2 files" begin
    
    # Create a matrix of random Taylor series
    M = Matrix{Taylor1{Float64}}(undef, 100, 100)
    for i in eachindex(M)
        M[i] = Taylor1(rand(25), 25)
    end 

    JLD2.save("test.jld2", "M", M)
    recovered_M = JLD2.load("test.jld2", "M")
    rm("test.jld2")

    @test M == recovered_M

end 
