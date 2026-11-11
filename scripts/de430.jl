using ArgParse, PlanetaryEphemeris, Dates, JLD2, Printf

# NAIF ids of the Sun, Mercury, Venus, the Earth, the Moon, and the barycenters of Mars,
# Jupiter, Saturn, Uranus, Neptune and Pluto; i.e. the first 11 bodies of the ephemeris
const PLANETS_SPKIDS = [10, 199, 299, 399, 301, 4, 5, 6, 7, 8, 9]
# NAIF ids of the bodies in the JPL DE430 model, in the order of the integration
const DE430_SPKIDS = vcat(PLANETS_SPKIDS, 2_000_000 .+ de430_343ast_ids)

function parse_commandline()

    s = ArgParseSettings()

    # Program name (for usage & help screen)
    s.prog = "de430.jl"
    # Desciption (for help screen)
    s.description = "Integrate JPL DE430 ephemeris"

    @add_arg_table! s begin
        "--d0"
            help = "Initial date"
            arg_type = DateTime
            default = DateTime(2000, 1, 1, 12)
        "--q0"
            help = "Initial conditions file"
            arg_type = String
            default = joinpath(pkgdir(PlanetaryEphemeris), "data", "de430ic_2000Jan1.txt")
        "--nyears"
            help = "Number of years"
            arg_type = Float64
            default = 10.0
        "--maxsteps"
            help = "Maximum number of steps"
            arg_type = Int
            default = 1_000_000
        "--order"
            help = "Order of Taylor expansions with respect to time"
            arg_type = Int
            default = 25
        "--abstol"
            help = "Absolute tolerance"
            arg_type = Float64
            default = 1E-20
        "--parse_eqs"
            help = "Whether to use the specialized method of `jetcoeffs` or not"
            arg_type = Bool
            default = true
        "--bodyind"
            help = "Body indices in output"
            arg_type = UnitRange{Int}
            default = 1:(11+16)
        "--output", "-o"
            help = "Output (.jld2) ephemeris file"
            arg_type = String
    end

    s.epilog = """
        examples:\n
        \n
        # Multi-threaded\n
        julia -t 10 --project de430.jl --d0 "2000-01-01T12:00:00" --nyears 100.0 -o de430.jld2\n
        \n
        # Single-threaded\n
        julia --project de430.jl --d0 "2000-01-01T12:00:00" --nyears 100.0 -o de430.jld2\n
        \n
    """

    return parse_args(s)
end

computationtime(x::DateTime, y::DateTime) = @sprintf("%.2f", (y - x).value / 60_000)

printitle(s::AbstractString, d::AbstractString) = println(d ^ length(s),
    '\n', s, '\n', d ^ length(s))

function main()

    printitle("Integrate JPL DE430 ephemeris", "=")

    # Parse arguments from commandline
    parsed_args = parse_commandline()

    # Number of threads
    println("• Number of threads: ", Threads.nthreads())

    # Initial date
    d0::DateTime = parsed_args["d0"]
    println("• Initial date: ", d0)

    # Initial conditions file
    filename::String = parsed_args["q0"]
    println("• Initial conditions file: ", filename)

    # Number of years
    nyears::Float64 = parsed_args["nyears"]
    println("• Number of years: ", nyears)

    # Maximum number of steps
    maxsteps::Int = parsed_args["maxsteps"]
    println("• Maximum number of steps: ", maxsteps)

    # Order of Taylor polynomials
    order::Int = parsed_args["order"]
    println("• Order of Taylor polynomials: ", order)

    # Absolute tolerance
    abstol::Float64 = parsed_args["abstol"]
    println("• Absolute tolerance: ", abstol)

    # Wheter to use @taylorize or not
    parse_eqs::Bool = parsed_args["parse_eqs"]
    println("• Use @taylorize: ", parse_eqs)

    # Body indices in output
    bodyind::UnitRange{Int} = parsed_args["bodyind"]
    println("• Body indices in output: ", bodyind)

    # Output file (checked before integrating, so a long run is not lost)
    isnothing(parsed_args["output"]) && error("Missing output file (--output)")
    output::String = parsed_args["output"]
    println("• Output file: ", output)

    # Planetary ephemeris problem
    jd0 = datetime2julian(d0)
    tspan = (jd0, jd0 + nyears * yr)
    q0 = read_initial_conditions(filename)
    params = DE430Params(jd0, q0, order)
    PE = PlanetaryEphemerisProblem(DE430!, tspan, q0, params)

    printitle("Integrator warmup", "-")
    @time propagate(PE; maxsteps = 1, order, abstol, parse_eqs)

    printitle("Full integration", "-")
    println( "Initial time of integration: ", string(julian2datetime(PE.tspan[1])) )
    println( "Final time of integration: ", string(julian2datetime(PE.tspan[2])))
    @time sol = propagate(PE; maxsteps, order, abstol, parse_eqs)

    # Bodies in output, lunar mantle and core Euler angles, and TT-TDB
    sseph = selecteph(sol, bodyind, euler = true, ttmtdb = true)
    # NAIF ids and gravitational parameters [au^3/day^2] of the bodies in output
    spkids = DE430_SPKIDS[bodyind]
    μ = PlanetaryEphemeris.μ[bodyind]
    # Ephemeris model and model parameters used by NEOs.jl, with the names and units of
    # the JPL DE430 technical comments file (same as scripts/spice_to_pe.jl): Earth's
    # equatorial radius [km] and second and third zonal harmonics
    model = "DE430"
    params = Dict{String, Float64}("RE" => RE, "J2E" => J2E,
        "J3E" => PlanetaryEphemeris.J3E)

    # Save results
    jldsave(output; sseph, spkids, μ, model, params)
    println("• Saved output to: ", output)

    return nothing
end

main()
