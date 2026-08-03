using ArgParse, Ephemerides, TaylorIntegration
using PlanetaryEphemeris, Dates, JLD2, Printf
using Ephemerides: spk_links, element_id, initial_time, final_time,
                   get_daf, file_id, get_segment, segment_list, header

# Conversion factors
const AU = 1.495978707E8
const DAYSEC = 86_400

# Included bodies IDs
const PLANETS_SPKIDS = [10, 199, 299, 399, 301, 4, 5, 6, 7, 8, 9]
const TTMTDB_SPKID = 1000000001

# Regular expressions
const FLOAT_REGEX = r"(?<gm>[0-9]\.[0-9]{16,18}[eED][-+][0-9]{2})"
const MA_REGEX = r"MA(?<id>[0-9]{4})\s+" * FLOAT_REGEX

#= DE430/DE431
const PERTURBERS_SPKIDS = [
    10, 199, 299, 399, 301, 4, 5, 6, 7, 8, 9,
    2000001, 2000004, 2000002, 2000010,
    2000031, 2000704, 2000511, 2000015,
    2000003, 2000016, 2000065, 2000088,
    2000048, 2000052, 2000451, 2000087,
    1000000001
]
=#
#= DE440/DE441
const PERTURBERS_SPKIDS = [
    10, 199, 299, 399, 301, 4, 5, 6, 7, 8, 9,
    2000001, 2000004, 2000002, 2000010,
    2000511, 2000704, 2000052, 2000087,
    2000015, 2000003, 2000016, 2000107,
    2000088, 2000007, 2000031, 2000065,
    1000000001
]
=#

function parse_commandline()

    s = ArgParseSettings()

    # Program name (for usage & help screen)
    s.prog = "de4xx2jld2.jl"
    # Desciption (for help screen)
    s.description = "Convert a JPL .bsp planetary ephemeris kernel \
        into a .jld2 file compatible with PlanetaryEphemeris.jl"

    @add_arg_table! s begin
        "--tech", "-t"
            help = "tech comments file"
            arg_type = String
        "--bsp", "-b"
            help = "planetary ephemeris .bsp kernel"
            arg_type = String
        "--output", "-o"
            help = "output .jld2 file"
            arg_type = String
        "--start", "-s"
            help = "start of interval [TDB]"
            arg_type = DateTime
            default = typemin(DateTime)
        "--end", "-e"
            help = "end of interval [TDB]"
            arg_type = DateTime
            default = typemax(DateTime)
        "--order", "-d"
            help = "order of Taylor expansions with respect to time"
            arg_type = Int
            default = 25
        "--abstol", "-a"
            help = "absolute tolerance"
            arg_type = Float64
            default = 1E-20
    end

    s.epilog = """
        examples:\n
        \n
        julia --project de4xx2jld2.jl -t de440_tech-comments.txt -b de440.bsp -o de440.jld2\n
        \n
    """

    return parse_args(s)
end

computationtime(x::DateTime, y::DateTime) = @sprintf("%.2f", (y - x).value / 60_000)

printitle(s::AbstractString, d::AbstractString) = println(d ^ length(s),
    '\n', s, '\n', d ^ length(s))

function spkid2gm(id::Int)
    if id == 10
        return "GMS"
    elseif id == 199 || id == 299 || id == 399
        i = digits(id)[end]
        return "GM$i"
    elseif id == 301
        return "GMM"
    elseif 4 ≤ id ≤ 9
        return "GM$id"
    else
        throw(ArgumentError("Invalid id ($id)"))
    end
end

function globalproperties(eph::EphemerisProvider, spkids::AbstractVector{Int})
    links = spk_links(eph)
    order = 0
    dt = typemax(Float64)
    et0, etf = typemax(Float64), typemin(Float64)
    fromtos = Vector{Vector{Int}}(undef, length(spkids))
    for (i, to) in enumerate(spkids)
        @assert haskey(links, to)
        if haskey(links[to], 1000000000)
            sublinks = links[to][1000000000]
            fromtos[i] = [1000000000, to]
        elseif haskey(links[to], 0)
            sublinks = links[to][0]
            fromtos[i] = [0, to]
        else
            mid = first(keys(links[to]))
            @assert haskey(links[mid], 0)
            sublinks = vcat(links[to][mid], links[mid][0])
            fromtos[i] = [0, mid, to]
        end
        for link in sublinks
            eid = element_id(link)
            et0 = min(et0, initial_time(link))
            etf = max(etf, final_time(link))
            daf = get_daf(eph, file_id(link))
            seg = get_segment(segment_list(daf), 1, eid)
            head = header(seg)
            dt = min(dt, head.tlen)
            order = max(order, head.order)
        end
    end
    return order, dt / DAYSEC, et0 / DAYSEC, etf / DAYSEC, fromtos
end

function evaleph(eph::EphemerisProvider, fromtos::AbstractVector, t::Number)
    et = t * DAYSEC
    N = length(fromtos) - 1
    y = Vector{typeof(t)}(undef, 6N + 1)
    for (i, fromto) in enumerate(fromtos)
        if fromto[end] == 1000000001
            y[end] = ephem_vector6(eph, fromto[1], fromto[2], et)[1]
        else
            idxs = nbodyind(N, i)
            y[idxs] = kmsec2auday(ephem_vector6(eph, fromto[1], fromto[2], et))
            for j in 3:length(fromto)
                y[idxs] += kmsec2auday(ephem_vector6(eph, fromto[j-1], fromto[j], et))
            end
        end
    end
    return y
end

function evaleph!(dq, q, params, t)
    local eph, fromtos = params
    local F = evaleph(eph, fromtos, t)
    for i in eachindex(dq)
        dq[i] = PE.ordpres_differentiate(F[i])
    end
    return nothing
end

function TaylorIntegration._allocate_jetcoeffs!(
        ::Val{evaleph!}, t::Taylor1{_T}, q::AbstractArray{Taylor1{_S}, _N},
        dq::AbstractArray{Taylor1{_S}, _N}, params
    ) where {_T <: Real, _S <: Number, _N}
    return TaylorIntegration.RetAlloc{Taylor1{_S}}(
        Taylor1{_S}[], [Array{Taylor1{_S}, 1}(undef, 0)], [Array{Taylor1{_S}, 2}(undef, 0, 0)],
        [Array{Taylor1{_S}, 3}(undef, 0, 0, 0)], [Array{Taylor1{_S}, 4}(undef, 0, 0, 0, 0)]
    )
end

function TaylorIntegration.jetcoeffs!(
        ::Val{evaleph!}, t::Taylor1{_T}, q::AbstractArray{Taylor1{_S}, _N},
        dq::AbstractArray{Taylor1{_S}, _N}, params,
        __ralloc::TaylorIntegration.RetAlloc{Taylor1{_S}}
    ) where {_T <: Real, _S <: Number, _N}
    order = TS.order(t)
    local eph, fromtos = params
    local F = evaleph(eph, fromtos, t)
    for ord = 0:order - 1
        ordnext = ord + 1
        for i = eachindex(dq)
            TS.differentiate!(dq[i], F[i], ord)
        end
        for __idx = eachindex(q)
            TaylorIntegration.solcoeff!(q[__idx], dq[__idx], ordnext)
        end
    end
    return nothing
end

function main()
    # Parse arguments from commandline
    parsed_args = parse_commandline()

    printitle("Create a JLD2 file for JPL DE4xx Planetary Ephemerides", "=")

    # Global initial time
    global_initial_time = now()
    println("• Run started at ", global_initial_time)

    # Number of workers and threads
    println("• Detected 1 worker with ", Threads.nthreads(), " thread(s)")

    # Tech comments file
    techfile::String = parsed_args["tech"]
    println("• Tech comments file: ", techfile)

    # Planetary ephemeris .bsp kernel
    kernelfile::String = parsed_args["bsp"]
    println("• Planetary ephemeris .bsp kernel: ", kernelfile)

    # Output .jld2 file
    output::String = parsed_args["output"]
    println("• Output .jld2 file: ", output)

    # Parse planets mass parameters
    text = read(techfile, String)
    spkids = Vector{Int}(undef, length(PLANETS_SPKIDS))
    μ = Vector{Float64}(undef, length(PLANETS_SPKIDS))
    for (i, spkid) in enumerate(PLANETS_SPKIDS)
        name = spkid2gm(spkid)
        re = Regex("$name\\s+") * FLOAT_REGEX
        m = match(re, text)
        spkids[i] = spkid
        μ[i] = parse(Float64, replace(m["gm"], r"[eED]" => 'E'))
    end

    # Parse perturbers mass parameters
    ms = Set{Tuple{Int, Float64}}()
    for m in eachmatch(MA_REGEX, text)
        spkid = 2000000 + parse(Int, m["id"])
        gm = parse(Float64, replace(m["gm"], r"[eED]" => 'E'))
        push!(ms, (spkid, gm))
    end
    filter!(x -> x[1] < 2008000, ms)
    ps = partialsort(collect(ms), 1:16, by = last, rev = true)
    append!(spkids, first.(ps), [TTMTDB_SPKID])
    append!(μ, last.(ps))
    println("• Included bodies IDs: \n", spkids)
    println("• Included bodies mass parameters: \n", μ)

    # Load planetary ephemerides
    eph = EphemerisProvider(kernelfile)
    # Global properties
    gps = globalproperties(eph, spkids)
    fromtos = gps[end]

    order::Int = parsed_args["order"]
    order = max(order, gps[1])
    println("• Order of Taylor expansions with respect to time: ", order)

    abstol::Float64 = parsed_args["abstol"]
    println("• Abstolute tolerance: ", abstol)

    maxstepsize::Float64 = gps[2]
    println("• Maximum allowed timestep: ", maxstepsize)

    d0::DateTime = parsed_args["start"]
    df::DateTime = parsed_args["end"]
    t0 = max(gps[3], datetime2julian(d0) - PE.J2000)
    tf = min(gps[4], datetime2julian(df) - PE.J2000)
    d0 = julian2datetime(t0 + PE.J2000)
    df = julian2datetime(tf + PE.J2000)
    println("• Ephemerides time interval [TDB]: ", d0, " - ", df)

    # Solar system ephemerides
    q0 = evaleph(eph, fromtos, t0)
    params = (eph, fromtos)
    maxsteps = ceil(Int, (tf - t0) / maxstepsize)
    sseph = taylorinteg(evaleph!, q0, t0, tf, order, abstol, params; maxstepsize, maxsteps)

    # Acceleration and Newtonian potential ephemerides
    acceph, poteph = PlanetaryEphemeris.loadeph(sseph, μ)

    # Save output
    jldsave(output; spkids, μ, sseph, acceph, poteph)

    # Final time
    global_final_time = now()
    println("• Run started ", global_initial_time, " and finished ", global_final_time)
    global_computation_time = computationtime(global_initial_time, global_final_time)
    println("• Total computation time was: ", global_computation_time, " min")

    return nothing
end

main()