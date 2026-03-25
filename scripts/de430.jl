using ArgParse, PlanetaryEphemeris, Dates, Printf

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
    end

    s.epilog = """
        examples:\n
        \n
        # Multi-threaded\n
        julia -t 10 --project de430.jl --d0 "2000-01-01T12:00:00" --nyears 100.0\n
        \n
        # Single-threaded\n
        julia --project de430.jl --d0 "2000-01-01T12:00:00" --nyears 100.0\n
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
    println("• Body indices in output: ", parse_eqs)

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

    # Save results
    selecteph2jld2(sol, bodyind, nyears)

    return nothing
end

main()
