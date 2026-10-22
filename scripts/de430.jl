using ArgParse, PlanetaryEphemeris, Dates, Printf

# Parse a vector of body indices given as a comma-separated list of indices and
# ranges, e.g. "1:27", "1,2,3", "1:11,55,66" or "[1:2:11, 55]"
function ArgParse.parse_item(::Type{AbstractVector{Int}}, x::AbstractString)
    v = Int[]
    for item in split(strip(x, ['[', ']', ' ']), r"[,\s]+"; keepempty = false)
        r = parse.(Int, split(item, ':'))
        if length(r) == 1
            push!(v, r[1])
        elseif length(r) == 2
            append!(v, r[1]:r[2])
        elseif length(r) == 3
            append!(v, r[1]:r[2]:r[3])
        else
            throw(ArgumentError("Invalid body indices: $x"))
        end
    end
    return v
end

function parse_commandline()

    s = ArgParseSettings()

    # Program name (for usage & help screen)
    s.prog = "de430.jl"
    # Desciption (for help screen)
    s.description = "Integrate JPL DE430 ephemeris"

    @add_arg_table! s begin
        "--params", "-p"
            help = "Parameters file"
            arg_type = String
            default = joinpath(pkgdir(PlanetaryEphemeris), "data", "de430params.txt")
        "--initcond", "-i"
            help = "Initial conditions file"
            arg_type = String
            default = joinpath(pkgdir(PlanetaryEphemeris), "data", "de430ic.txt")
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
            help = "Body indices in output, e.g. 1:27 or 1:11,55,66 (default: Sun, \
                planets, Moon, Pluto and the 16 most massive asteroids)"
            arg_type = AbstractVector{Int}
            default = 1:(11+16)
    end

    s.epilog = """
        examples:\n
        \n
        # Multi-threaded\n
        julia -t 10 --project de430.jl -p ../data/de430params.txt -i ../data/de430ic.txt --nyears 100.0\n
        \n
        # Single-threaded\n
        julia --project de430.jl -p ../data/de430params.txt -i ../data/de430ic.txt --nyears 100.0\n
        \n
        # Save the Sun, planets, Moon, Pluto and Ceres only\n
        julia -t 10 --project de430.jl --nyears 100.0 --bodyind 1:12\n
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

    # Parameters file
    paramsfile::String = parsed_args["params"]
    println("• Parameters file: ", paramsfile)

    # Initial conditions file
    icfile::String = parsed_args["initcond"]
    println("• Initial conditions file: ", icfile)

    # Initial conditions and their epoch [JD TDB]
    jd0, q0 = read_initial_conditions(icfile)
    println("• Initial date: ", julian2datetime(jd0))

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
    bodyind::AbstractVector{Int} = parsed_args["bodyind"]
    println("• Body indices in output: ", bodyind)

    # Planetary ephemeris problem
    tspan = (jd0, jd0 + nyears * yr)
    consts = DE430Constants{Float64}(paramsfile)
    params = DE430Params(jd0, q0, order; parse_eqs, consts)
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
