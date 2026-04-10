using ArgParse, PlanetaryEphemeris, Dates, TaylorSeries, LinearAlgebra
using SPICE, Printf, JLD2
using PlanetaryEphemeris: loadeph

const DensePropagation2{T, U} = TaylorInterpolant{T, U, 2, Vector{T}, Matrix{Taylor1{U}}}

function parse_commandline()

    s = ArgParseSettings()

    # Program name (for usage & help screen)
    s.prog = "spk2jld2.jl"
    # Desciption (for help screen)
    s.description = "Convert a JPL .bsp planetary ephemeris kernel \
        into a .jld2 file compatible with PlanetaryEphemeris.jl"

    @add_arg_table! s begin
        "--start", "-s"
            help = "start of interval [TDB]"
            arg_type = DateTime
        "--end", "-e"
            help = "end of interval [TDB]"
            arg_type = DateTime
        "--input", "-i"
            help = "input .bsp kernel"
            arg_type = String
            default = "de430.bsp"
        "--output", "-o"
            help = "output .jld2 file"
            arg_type = String
            default = "de430.jld2"
    end

    s.epilog = """
        examples:\n
        \n
        julia --project spk2jld2.jl -i de430.bsp -o de430.jld2\n
        \n
    """

    return parse_args(s)
end

computationtime(x::DateTime, y::DateTime) = @sprintf("%.2f", (y - x).value / 60_000)

printitle(s::AbstractString, d::AbstractString) = println(d ^ length(s),
    '\n', s, '\n', d ^ length(s))

function increase_order(x::Taylor1, order::Int)
    y = Taylor1(TS.numtype(x), order)
    TS.zero!(y)
    y.coeffs[1:min(get_order(x), order)+1] .= x.coeffs
    return y
end

function main()

    # Parse arguments from commandline
    parsed_args = parse_commandline()

    printitle("Convert a JPL kernel into a file compatible with PlanetaryEphemeris.jl", "=")

    # Number of workers and threads
    println("• Detected 1 worker with ", Threads.nthreads(), " thread(s)")

    # Start of interval
    d0::DateTime = parsed_args["start"]
    println("• Start of interval: ", d0, " TDB")

    # End of interval
    df::DateTime = parsed_args["end"]
    println("• End of interval: ", df, " TDB")

    # Input .bsp kernel
    input::String = parsed_args["input"]
    println("• Input .bsp kernel: ", input)

    # Output .jld2 file
    output::String = parsed_args["output"]
    println("• Output .jld2 file: ", output)

    # Global initial time
    global_initial_time = now()
    println("• Run started at ", global_initial_time)

    # Generate Chebyshev polynomials
    order = 20
    t = Taylor1(order)
    chebyshev = Vector{Taylor1{Int}}(undef, order+1)
    chebyshev[1] = one(t)
    chebyshev[2] = t
    for j in 3:order+1
        chebyshev[j] = 2 * chebyshev[2] * chebyshev[j-1] - chebyshev[j-2]
    end

    # Pre-allocate memory
    dict = Dict{Tuple{Int, Int}, DensePropagation2{Float64, Float64}}()

    # Open the DAF file for reading
    handle = dafopr(input)

    # Begin a forward search through the segments
    dafbfs(handle)

    # Find the next (forward) array in the current DAF
    while daffna()
        # Get the summary of the current segment
        sum = dafgs()
        # Unpack an array summary into its double precision and integer components
        # da[1], da[2] = initial/final epoch [TDB seconds since J2000]
        # ia[1] = target ID, ia[2] = center ID,
        # ia[3] = frame ID, ia[4] = SPK type,
        # ia[5], da[6] = starting/ending location of the data in the DAF
        da, ia = dafus(sum, 2, 6)
        target, center = ia[1], ia[2]

        # Get raw data between the start and end addresses
        # These are the actual coefficients and the time metadata
        start_addr, end_addr = ia[5:6]
        raw_data = dafgda(handle, start_addr, end_addr)

        # For SPK Type 2 (Fixed-length Chebyshev polynomials):
        # A four-number `directory' at the end of the segment contains the
        # information needed to determine the location of the record corresponding
        # to a particular epoch.
        # 1. INIT is the initial epoch of the first record, given in ephemeris seconds past J2000.
        # 2. INTLEN is the length of the interval covered by each record, in seconds.
        # 3. RSIZE is the total size of (number of array elements in) each record.
        # 4. N is the number of records contained in the segment.
        INIT = raw_data[end-3]
        INTLEN = raw_data[end-2]
        RSIZE = Int(raw_data[end-1])
        N = Int(raw_data[end])
        order = Int( ( RSIZE - 2 ) / 3 - 1 )

        # Initial time
        t0 = INIT / daysec

        # Vector of times
        times = Vector{Float64}(undef, N+1)
        times[1] = zero(Float64)
        for i in 2:N+1
            times[i] = times[i-1] + INTLEN
        end
        @. times = times / daysec

        # Matrix of coefficients
        Cs = chebyshev[1:order+1](2 * (Taylor1(order) * daysec / INTLEN) - 1)
        M = Matrix{Taylor1{Float64}}(undef, target == 1000000001 ? 1 : 6, N)
        for (j, x) in enumerate(Iterators.partition(raw_data, RSIZE))
            length(x) < RSIZE && continue
            # MID, RADIUS = x[1], x[2]
            if target == 1000000001
                M[1, j] = dot(x[3+0*order:3+1*order], Cs)
            else
                M[1, j] = dot(x[3+0*order:3+1*order], Cs) / au
                M[2, j] = dot(x[4+1*order:4+2*order], Cs) / au
                M[3, j] = dot(x[5+2*order:5+3*order], Cs) / au
                @. M[4:6, j] = PE.ordpres_differentiate(M[1:3, j])
            end
        end
        dict[(target, center)] = TaylorInterpolant(t0, times, collect(transpose(M)))
    end

    # Close the DAF file
    dafcls(handle)

    # Make all the interpolants of the same order
    order = maximum(get_order, values(dict))
    for TI in values(dict)
        @. TI.x = increase_order(TI.x, order)
    end

    # Initial and final times [TDB seconds since J2000]
    t0 = datetime2julian(d0) - J2000
    tf = datetime2julian(df) - J2000
    # Vector of times
    times = mapreduce(x -> x.t0 .+ x.t, vcat, values(dict))
    unique!(times)
    sort!(times)
    clamp!(times, t0, tf)
    unique!(times)
    sort!(times)
    @. times = times - t0
    # Expand every interpolant at every t in times
    for (i, TI) in dict
        M = Matrix{Taylor1{Float64}}(undef, size(TI.x, 2), length(times)-1)
        for j in axes(M, 2)
            M[:, j] .= TI(t0 + times[j] + Taylor1(order))
        end
        dict[i] = TaylorInterpolant(t0, times, collect(transpose(M)))
    end

    # Bodies to be included in the output
    spkids = [
        10, 199, 299, 399, 301, 4, 5, 6, 7, 8, 9, 2000001, 2000004, 2000002,
        2000010, 2000031, 2000704, 2000511, 2000015, 2000003, 2000016,
        2000065, 2000088, 2000048, 2000052, 2000451, 2000087
    ]
    # Positions and velocities
    N = length(spkids)
    M = Matrix{Taylor1{Float64}}(undef, length(times) - 1, 6N+13)
    for (i, id) in enumerate(spkids)
        if haskey(dict, (id, 0))
            @. M[:, nbodyind(N, i)] = dict[(id, 0)].x
        elseif haskey(dict, (id, 10))
            @. M[:, nbodyind(N, i)] = dict[(id, 10)].x + dict[(10, 0)].x
        else
            jd = digits(id)[end]
            @. M[:, nbodyind(N, i)] = dict[(id, jd)].x + dict[(jd, 0)].x
        end
    end
    # Lunar Euler angles
    zero_q = zero(M[1])
    for j in 6N+1:6N+12
        for i in axes(M, 1)
            M[i, j] = zero(zero_q)
        end
    end
    # TT-TDB
    @. M[:, end] = dict[((1000000001, 1000000000))].x
    # Assemble the global TaylorInterpolant
    sseph = TaylorInterpolant(t0, times, M)

    # Save output
    jldsave(output; ss16ast_eph = sseph)
    println("• Saved output to: ", output)

    # Final time
    global_final_time = now()
    println("• Run started ", global_initial_time, " and finished ", global_final_time)
    global_computation_time = computationtime(global_initial_time, global_final_time)
    println("• Total computation time was: ", global_computation_time, " min")

    return nothing
end

main()