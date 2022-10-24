# Multi-threaded:
# JULIA_NUM_THREADS=<number-of-threads> julia --project=@. integrate_ephemeris.jl --help
# Single-threaded:
# julia --project=@. integrate_ephemeris.jl --help  

using ArgParse, PlanetaryEphemeris, Dates

function parse_commandline()
    s = ArgParseSettings()

    # Program name (for usage & help screen)
    s.prog = "integrate_ephemeris.jl"  
    # Desciption (for help screen)
    s.description = "Integrates JPL DE430 Ephemeris" 

    @add_arg_table! s begin
        "--maxsteps"
            help = "Maximum number of steps during integration"
            arg_type = Int
            default = 1_000_000 
        "--jd0"
            help = "Starting time of integration"
            arg_type = Float64
            # default = datetime2julian(DateTime(1969,6,28,0,0,0))
            default = datetime2julian(DateTime(2000,1,1,12))
        "--dense"
            help = "Whether to output the Taylor polynomial solutions obtained at each time step or not"
            arg_type = Bool
            default = true
        "--dynamics"
            help = "Dynamical model function"
            arg_type = Function
            default = DE430!
        "--nast"
            help = "Number of asteroid perturbers"
            arg_type = Int
            default = 343 # 16
        "--quadmath"
            help = "Whether to use quadruple precision or not"
            arg_type = Bool
            default = false
        "--bodyind"
            help = "Body indices in output"
            arg_type = UnitRange{Int}
            default = 1:(11+16) # 1:(11+nast)
        "--order"
            help = "Order of Taylor polynomials expansions during integration"
            arg_type = Int
            default = 25
        "--abstol"
            help = "Absolute tolerance"
            arg_type = Float64
            default = 1.0E-20
    end

    return parse_args(s)
end

function main()
    parsed_args = parse_commandline()

    println("Number of threads: ", Threads.nthreads())

    nyears = 2031.0 - year(julian2datetime(parsed_args["jd0"]))
    println("jd0: ", parsed_args["jd0"])
    println("J2000: ", J2000)
    println("jd0 - J2000: ", parsed_args["jd0"]-J2000)
    println("dynamics: ", parsed_args["dynamics"])

    # bodyind = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 17, 18, 19, 20, 21, 22, 23, 25, 27, 28, 30, 40, 41, 46, 55, 62, 73, 113, 115, 277, 322] # SS + 25 ast perturbers

    # Integrator warmup
    PlanetaryEphemeris.propagate(1, parsed_args["jd0"], nyears, output = false, 
                                 dense = parsed_args["dense"], dynamics = parsed_args["dynamics"],
                                 nast = parsed_args["nast"], quadmath = parsed_args["quadmath"],
                                 bodyind = parsed_args["bodyind"], order = parsed_args["order"], 
                                 abstol = parsed_args["abstol"])
    println("*** Finished warmup")

    # Perform full integration
    PlanetaryEphemeris.propagate(parsed_args["maxsteps"], parsed_args["jd0"], nyears, 
                                 dense = parsed_args["dense"], dynamics = parsed_args["dynamics"],
                                 nast = parsed_args["nast"], quadmath = parsed_args["quadmath"],
                                 bodyind = parsed_args["bodyind"], order = parsed_args["order"],
                                 abstol = parsed_args["abstol"])
    println("*** Finished full integration")

end

main()
