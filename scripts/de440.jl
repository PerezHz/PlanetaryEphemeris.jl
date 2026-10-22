using ArgParse, PlanetaryEphemeris, Dates, JLD2, Printf

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

# Indices of the Sun, planets, Moon and Pluto, followed by the indices of the 16 most
# massive asteroids in data/de440ic.txt (bodies 12 to 420 are sorted by decreasing mass
# parameter, but include KBOs and the Kuiper belt masses):
# Ceres, Vesta, Pallas, Hygiea, Davida, Interamnia, Europa, Sylvia, Eunomia, Juno,
# Psyche, Camilla, Thisbe, Iris, Euphrosyne and Cybele
const DE440_BODYIND = vcat(1:11, [55, 66, 69, 80, 81, 82, 83, 84, 85, 86, 88, 89, 90, 91,
    92, 93])

function parse_commandline()

    s = ArgParseSettings()

    # Program name (for usage & help screen)
    s.prog = "de440.jl"
    # Desciption (for help screen)
    s.description = "Integrate JPL DE440 ephemeris"

    @add_arg_table! s begin
        "--params", "-p"
            help = "Parameters file"
            arg_type = String
            default = joinpath(pkgdir(PlanetaryEphemeris), "data", "de440params.txt")
        "--initcond", "-i"
            help = "Initial conditions file"
            arg_type = String
            default = joinpath(pkgdir(PlanetaryEphemeris), "data", "de440ic.txt")
        "--output", "-o"
            help = "Output (.jld2) ephemeris file"
            arg_type = String
            required = true
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
            default = DE440_BODYIND
    end

    s.epilog = """
        examples:\n
        \n
        # Multi-threaded\n
        julia -t 10 --project de440.jl -p ../data/de440params.txt -i ../data/de440ic.txt --nyears 100.0 -o de440.jld2\n
        \n
        # Single-threaded\n
        julia --project de440.jl -p ../data/de440params.txt -i ../data/de440ic.txt --nyears 100.0 -o de440.jld2\n
        \n
        # Save the Sun, planets, Moon, Pluto and Ceres only\n
        julia -t 10 --project de440.jl --nyears 100.0 -o de440.jld2 --bodyind 1:11,55\n
        \n
    """

    return parse_args(s)
end

computationtime(x::DateTime, y::DateTime) = @sprintf("%.2f", (y - x).value / 60_000)

printitle(s::AbstractString, d::AbstractString) = println(d ^ length(s),
    '\n', s, '\n', d ^ length(s))

function main()

    printitle("Integrate JPL DE440 ephemeris", "=")

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

    # Output ephemeris file
    output::String = parsed_args["output"]
    println("• Output ephemeris file: ", output)

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
    consts = DE440Constants{Float64}(paramsfile)
    params = DE440Params(jd0, q0, order; parse_eqs, consts)
    PE = PlanetaryEphemerisProblem(DE440!, tspan, q0, params)

    printitle("Integrator warmup", "-")
    @time propagate(PE; maxsteps = 1, order, abstol, parse_eqs)

    printitle("Full integration", "-")
    println( "Initial time of integration: ", string(julian2datetime(PE.tspan[1])) )
    println( "Final time of integration: ", string(julian2datetime(PE.tspan[2])))
    @time sol = propagate(PE; maxsteps, order, abstol, parse_eqs)

    # Save results
    # Select the `bodyind`-th bodies and the lunar mantle and core Euler angles and
    # angular velocities (DE440 does not integrate TT-TDB)
    sseph = selecteph(sol, bodyind; euler = true, ttmtdb = false)
    println("Saving ephemeris to file: ", output)
    jldsave(output; sseph)

    return nothing
end

main()
