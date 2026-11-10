using ArgParse, Dates, JLD2, Printf, TaylorSeries
using PlanetaryEphemeris
using PlanetaryEphemeris: au, daysec, nbodyind

#=
Convert JPL planetary ephemeris SPICE kernels (.bsp, and optionally the lunar
orientation binary PCK .bpc) into a .jld2 file compatible with PlanetaryEphemeris.jl
and NEOs.jl.

Method
------
SPICE kernels store each body's ephemeris as fixed-length Chebyshev records, with
time grids and degrees that differ from body to body. This script:

1. Reads the raw Chebyshev records directly from the DAF files (SPK types 2 and 3,
   binary PCK type 2); no SPICE toolkit or Ephemerides.jl is needed.
2. Builds a common time grid given by the union of the record boundaries (and segment
   boundaries) of every segment used, so no Taylor polynomial ever straddles a
   Chebyshev record.
3. On each grid interval, converts every Chebyshev record into a Taylor polynomial
   about the left endpoint of the interval. This is an exact change of polynomial
   basis (up to floating point rounding); positions are summed along the SPK chain
   of each body (e.g. Earth = SSB->EMB + EMB->Earth), velocities are the analytical
   derivatives (or the velocity records for type 3 segments).
4. Truncates every polynomial at the highest Chebyshev degree found in the kernels
   (optionally zero-padded to a larger order with --order).

Output layout (same as the ephemerides integrated by PlanetaryEphemeris.jl)
--------------------------------------------------------------------------------
Columns of `ss16ast_eph.p` (N bodies, see `spkids`):
    1:3N        barycentric positions [au]
    3N+1:6N     barycentric velocities [au/day]
    6N+1:6N+3   lunar mantle Euler angles φ, θ, ψ (3-1-3, ICRF -> PA) [rad]
    6N+4:6N+6   lunar mantle angular velocity in the mantle frame [rad/day];
                not a polynomial in time, so it is approximated on each
                interval by Chebyshev interpolation
    6N+7:6N+12  lunar core Euler angles and angular velocity (not available in
                SPICE kernels; filled with zeros)
    6N+13       TT-TDB [s]
Time `ss16ast_eph.t` is in TDB days since J2000. Columns 6N+1:6N+6 are zero when no
lunar PCK is given.

Additional keys in the output file:
    spkids::Vector{Int}      NAIF ids of the N bodies, in ephemeris order
    μ::Vector{Float64}       gravitational parameters [au^3/day^2], same order
    model::String            name of the ephemeris model (e.g. "DE440")
    params::Dict{String, Float64}  model parameters used by NEOs.jl, with the names and
                             units of the technical comments file (see MODEL_PARAMS)
    kernels::Vector{String}  basenames of the input kernels
=#

# Bodies: Sun, Mercury, Venus, Earth, Moon, Mars, Jupiter, Saturn, Uranus, Neptune, Pluto
const PLANETS_SPKIDS = [10, 199, 299, 399, 301, 4, 5, 6, 7, 8, 9]
const SSB_SPKID = 0
const TTMTDB_SPKID = 1000000001
const TTMTDB_CENTER = 1000000000
const J2000_FRAME = 1
# DateTime of the J2000 epoch (TDB)
const J2000_DATETIME = DateTime(2000, 1, 1, 12)
# Breakpoints closer than this [s] are merged
const BREAKPOINT_TOL = 1E-3
# Model parameters written to the output file, with the names and units of the technical
# comments file: Earth's equatorial radius [km] and second and third zonal harmonics
const MODEL_PARAMS = ["RE", "J2E", "J3E"]

isasteroid(id::Integer) = (2_000_000 < id < 3_000_000) || (20_000_000 < id < 30_000_000)

function gmname(id::Int)
    id == 10 && return "GMS"
    id == 301 && return "GMM"
    id in (199, 299, 399) && return "GM$(id ÷ 100)"
    1 ≤ id ≤ 9 && return "GM$id"
    isasteroid(id) && return "MA" * lpad(id < 3_000_000 ? id - 2_000_000 :
        id - 20_000_000, 4, '0')
    throw(ArgumentError("No gravitational parameter name for body $id"))
end

function parse_commandline()
    s = ArgParseSettings()
    s.prog = "spice_to_pe.jl"
    s.description = "Convert JPL planetary ephemeris SPICE kernels into a .jld2 file \
        compatible with PlanetaryEphemeris.jl and NEOs.jl"
    @add_arg_table! s begin
        "--bsp", "-b"
            help = "SPK kernel(s); later kernels have higher priority (as in SPICE)"
            arg_type = String
            nargs = '+'
            required = true
        "--tech", "-t"
            help = "technical comments file (source of the gravitational parameters)"
            arg_type = String
            required = true
        "--pck", "-p"
            help = "lunar orientation binary PCK (e.g. moon_pa_de440_200625.bpc)"
            arg_type = String
        "--output", "-o"
            help = "output .jld2 file"
            arg_type = String
            required = true
        "--start", "-s"
            help = "start of interval [TDB]; default: start of common coverage"
            arg_type = DateTime
        "--end", "-e"
            help = "end of interval [TDB]; default: end of common coverage"
            arg_type = DateTime
        "--asteroids", "-a"
            help = "comma-separated NAIF ids of the asteroids to include, in output \
                order; default: every asteroid in the kernels, by decreasing GM"
            arg_type = String
        "--order", "-d"
            help = "order of the Taylor polynomials; default: highest Chebyshev \
                degree in the kernels (larger values are zero-padded)"
            arg_type = Int
            default = -1
        "--model", "-m"
            help = "name of the ephemeris model (e.g. DE440); default: inferred from \
                the technical comments file"
            arg_type = String
        "--check", "-c"
            help = "number of epochs used to validate the output against \
                direct Chebyshev evaluation of the kernels (0 to skip)"
            arg_type = Int
            default = 1000
    end
    s.epilog = """
        examples:\n
        \n
        julia --project spice_to_pe.jl -b de440.bsp -t de440_tech-comments.txt -o de440.jld2\n
        \n
        julia -t 4 --project spice_to_pe.jl -b de440t.bsp sb441-n16.bsp -t de440_tech-comments.txt -p moon_pa_de440_200625.bpc -s 1900-01-01 -e 2100-01-01 -o de440.jld2\n
    """
    return parse_args(s)
end

printitle(s::AbstractString, d::AbstractString) = println(d ^ length(s), '\n', s,
    '\n', d ^ length(s))

computationtime(x::DateTime, y::DateTime) = @sprintf("%.2f", (y - x).value / 60_000)

# TDB seconds since J2000 <-> DateTime (exact up to the DateTime millisecond resolution)
datetime2et(d::DateTime) = (d - J2000_DATETIME).value / 1_000
et2datetime(et::Real) = J2000_DATETIME + Millisecond(round(Int, 1_000 * et))

## DAF reader

# A segment of Chebyshev records of a DAF file (SPK types 2/3 or binary PCK type 2)
struct ChebyshevSegment
    file::String         # source file
    priority::Int        # higher priority segments win (SPICE convention)
    target::Int          # SPK: target body; PCK: frame class id
    center::Int          # SPK: center body; PCK: reference frame
    frame::Int           # reference frame
    type::Int            # SPK/PCK data type
    et0::Float64         # coverage start [TDB s since J2000]
    etf::Float64         # coverage end [TDB s since J2000]
    init::Float64        # start of the first record [TDB s since J2000]
    intlen::Float64      # length of each record [s]
    nrec::Int            # number of records
    ncomp::Int           # number of interpolated components
    ncoeff::Int          # number of Chebyshev coefficients per component
    records::Matrix{Float64} # (2 + ncomp*ncoeff) × nrec
end

degree(s::ChebyshevSegment) = s.ncoeff - 1
edge(s::ChebyshevSegment) = (s.target, s.center)

function readdaf(path::AbstractString, priority0::Int)
    segments = ChebyshevSegment[]
    skipped = Tuple{Int, Int, Int}[]
    kind = :unknown
    open(path, "r") do io
        rec = read(io, 1024)
        length(rec) == 1024 || error("$path is too short to be a DAF file")
        idword = String(rec[1:8])
        startswith(idword, "DAF/") || error("$path is not a DAF file (ID word: $idword)")
        kind = strip(idword[5:8]) == "SPK" ? :spk : strip(idword[5:8]) == "PCK" ?
            :pck : error("$path is neither an SPK nor a binary PCK file ($idword)")
        fmt = String(rec[89:96])
        swap = if fmt == "LTL-IEEE"
            ENDIAN_BOM != 0x04030201
        elseif fmt == "BIG-IEEE"
            ENDIAN_BOM == 0x04030201
        else
            error("$path has an unsupported binary format ($fmt)")
        end
        fix(x) = swap ? bswap(x) : x
        geti(b, o) = Int(fix(reinterpret(Int32, b[o+1:o+4])[1]))
        getd(b, o) = fix(reinterpret(Float64, b[o+1:o+8])[1])
        nd, ni = geti(rec, 8), geti(rec, 12)
        (nd == 2 && ni == (kind == :spk ? 6 : 5)) || error("$path has unexpected \
            summary format (ND = $nd, NI = $ni)")
        fward = geti(rec, 76)
        ss = nd + cld(ni, 2)
        recno = fward
        while recno != 0
            seek(io, 1024 * (recno - 1))
            srec = read(io, 1024)
            next, nsum = round(Int, getd(srec, 0)), round(Int, getd(srec, 16))
            for k in 0:nsum-1
                o = 24 + 8 * ss * k
                et0, etf = getd(srec, o), getd(srec, o + 8)
                ints = [geti(srec, o + 16 + 4j) for j in 0:ni-1]
                if kind == :spk
                    target, center, frame, type, a0, a1 = ints
                else
                    target, frame, type, a0, a1 = ints
                    center = frame
                end
                ncomp = (type == 2) ? 3 : (type == 3 && kind == :spk) ? 6 : 0
                if ncomp == 0
                    push!(skipped, (target, center, type))
                    continue
                end
                # Segment data (addresses are 1-based double precision words)
                seek(io, 8 * (a0 - 1))
                data = fix.(reinterpret(Float64, read(io, 8 * (a1 - a0 + 1))))
                init, intlen, rsize, nrec = data[end-3], data[end-2],
                    round(Int, data[end-1]), round(Int, data[end])
                (rsize - 2) % ncomp == 0 || error("Inconsistent record size \
                    ($rsize) in segment ($target, $center) of $path")
                ncoeff = (rsize - 2) ÷ ncomp
                records = reshape(data[1:rsize*nrec], rsize, nrec)
                push!(segments, ChebyshevSegment(path, priority0 + length(segments) + 1,
                    target, center, frame, type, et0, etf, init, intlen, nrec, ncomp,
                    ncoeff, records))
            end
            recno = next
        end
    end
    return kind, segments, skipped
end

# Piecewise selection of the highest priority segment of a (target, center) pair:
# segs[i] is used in [times[i], times[i+1]] (nothing if no segment covers it)
struct SegmentSchedule
    times::Vector{Float64}
    segs::Vector{Union{Nothing, ChebyshevSegment}}
end

function SegmentSchedule(segs::AbstractVector{ChebyshevSegment})
    times = sort!(unique!(vcat([s.et0 for s in segs], [s.etf for s in segs])))
    # Segments sorted by decreasing priority
    bypriority = sort(segs, by = s -> -s.priority)
    sel = Vector{Union{Nothing, ChebyshevSegment}}(undef, length(times) - 1)
    for i in eachindex(sel)
        a, b = times[i], times[i+1]
        j = findfirst(s -> s.et0 ≤ a && b ≤ s.etf, bypriority)
        sel[i] = isnothing(j) ? nothing : bypriority[j]
    end
    return SegmentSchedule(times, sel)
end

# Highest priority segment covering [a, b]
function findsegment(sch::SegmentSchedule, a::Real, b::Real = a)
    i = clamp(searchsortedlast(sch.times, a), 1, length(sch.segs))
    s = sch.segs[i]
    (isnothing(s) || !(s.et0 ≤ a && b ≤ s.etf)) && error("No segment covers \
        [$a, $b] s since J2000")
    return s
end

# Index of the record of `s` containing `et`
function recordindex(s::ChebyshevSegment, et::Real)
    k = floor(Int, (et - s.init) / s.intlen) + 1
    return clamp(k, 1, s.nrec)
end

# Chebyshev coefficients of component `comp` of record `k`
coeffs(s::ChebyshevSegment, k::Int, comp::Int) =
    view(s.records, 2 + (comp-1)*s.ncoeff + 1 : 2 + comp*s.ncoeff, k)

## Coverage

# Union of the coverage intervals of `segs`, as a sorted list of disjoint intervals
function coverage(segs::AbstractVector{ChebyshevSegment})
    ivs = sort!([(s.et0, s.etf) for s in segs])
    merged = [ivs[1]]
    for (a, b) in Iterators.drop(ivs, 1)
        if a ≤ merged[end][2]
            merged[end] = (merged[end][1], max(merged[end][2], b))
        else
            push!(merged, (a, b))
        end
    end
    return merged
end

function intersectcoverage(A, B)
    C = Tuple{Float64, Float64}[]
    for (a1, b1) in A, (a2, b2) in B
        a, b = max(a1, a2), min(b1, b2)
        a < b && push!(C, (a, b))
    end
    return sort!(C)
end

## SPK chains

# Sequence of (target, center) edges from `body` to the solar system barycenter
function spkchain(edges::Dict{Int, Vector{Int}}, body::Int)
    chain = Tuple{Int, Int}[]
    cur = body
    visited = Set{Int}()
    while cur != SSB_SPKID
        cur in visited && error("Cyclic SPK chain for body $body")
        push!(visited, cur)
        haskey(edges, cur) || error("The kernels have no segment for body $cur, \
            needed to compute the barycentric state of body $body")
        cs = unique(edges[cur])
        c = SSB_SPKID in cs ? SSB_SPKID : length(cs) == 1 ? cs[1] :
            error("Ambiguous center for body $cur: $cs")
        push!(chain, (cur, c))
        cur = c
    end
    return chain
end

## Gravitational parameters

const FLOAT_PATTERN = "([+-]?(?:\\d+\\.\\d*|\\.\\d+)(?:[eEdD][+-]?\\d+)?)"

function techconstant(text::AbstractString, name::AbstractString)
    re = Regex("(?<![A-Za-z0-9_])" * name * "(?![A-Za-z0-9_])\\s*=?\\s*" * FLOAT_PATTERN)
    values = [parse(Float64, replace(m[1], r"[dD]" => 'E')) for m in eachmatch(re, text)]
    isempty(values) && error("Constant $name not found in the technical comments file")
    allequal(values) || error("Constant $name has inconsistent values in the \
        technical comments file: $values")
    return first(values)
end

# Name of the ephemeris model, inferred from the name of the technical comments file
# `techfile` (e.g. de440_tech-comments.txt) or, failing that, from its contents `text`
function modelname(techfile::AbstractString, text::AbstractString)
    m = match(r"de(\d{3,4})"i, basename(techfile))
    isnothing(m) || return "DE" * m[1]
    ms = unique(m[1] for m in eachmatch(r"\bDE[- ]?(\d{3,4})\b", text))
    length(ms) == 1 && return "DE" * ms[1]
    error("Could not infer the ephemeris model from $techfile; use --model")
end

## Chebyshev -> Taylor

# Taylor coefficients (in τ) of Σ c[k+1] T_k(α + βτ); `p` must have length ≥ length(c)
function cheb2taylor!(p::Vector{Float64}, T0::Vector{Float64}, T1::Vector{Float64},
                      T2::Vector{Float64}, c::AbstractVector{Float64}, α::Float64,
                      β::Float64)
    L, n = length(p), length(c)
    @assert L ≥ n
    fill!(p, 0.0)
    p[1] = c[1]
    n == 1 && return p
    fill!(T0, 0.0); T0[1] = 1.0
    fill!(T1, 0.0); T1[1] = α; T1[2] = β
    fill!(T2, 0.0)
    @inbounds for j in 1:2
        p[j] += c[2] * T1[j]
    end
    @inbounds for k in 2:n-1
        # T_k = 2(α + βτ) T_{k-1} - T_{k-2} (degree k)
        T2[1] = 2α * T1[1] - T0[1]
        for j in 2:k+1
            T2[j] = 2 * (α * T1[j] + β * T1[j-1]) - T0[j]
        end
        for j in 1:k+1
            p[j] += c[k+1] * T2[j]
        end
        T0, T1, T2 = T1, T2, T0
    end
    return p
end

# Derivative of the polynomial with coefficients `p`, with the same length
function polyderivative(p::Vector{Float64})
    L = length(p)
    dp = zeros(L)
    for j in 1:L-1
        dp[j] = j * p[j+1]
    end
    return dp
end

struct Workspace
    T0::Vector{Float64}
    T1::Vector{Float64}
    T2::Vector{Float64}
end
Workspace(L::Int) = Workspace(zeros(L), zeros(L), zeros(L))

# Taylor expansions, about `a` and in TDB days, of the first `ncomp` components of
# `segs` over the interval [a, b] [TDB s since J2000]; also returns their derivatives
# (from the velocity records for SPK type 3 segments)
function expandedge(segs::SegmentSchedule, a::Float64, b::Float64, L::Int,
                    ws::Workspace; ncomp::Int = 3)
    s = findsegment(segs, a, b)
    k = recordindex(s, (a + b) / 2)
    MID, RADIUS = s.records[1, k], s.records[2, k]
    @assert MID - RADIUS - BREAKPOINT_TOL ≤ a && b ≤ MID + RADIUS + BREAKPOINT_TOL
    α, β = (a - MID) / RADIUS, daysec / RADIUS
    x = [cheb2taylor!(zeros(L), ws.T0, ws.T1, ws.T2, coeffs(s, k, i), α, β)
         for i in 1:ncomp]
    v = if s.ncomp == 6
        # SPK type 3 velocity records [km/s] -> [km/day]
        [cheb2taylor!(zeros(L), ws.T0, ws.T1, ws.T2, coeffs(s, k, 3 + i), α, β) .* daysec
         for i in 1:ncomp]
    else
        [polyderivative(x[i]) for i in 1:ncomp]
    end
    return x, v
end

## Direct evaluation (used for validation)

# Value and derivative (wrt x) of Σ c[k+1] T_k(x)
function chebvalder(c::AbstractVector{Float64}, x::Float64)
    T0, T1 = 1.0, x
    dT0, dT1 = 0.0, 1.0
    f, df = c[1], 0.0
    if length(c) > 1
        f += c[2] * T1
        df += c[2] * dT1
    end
    for k in 3:length(c)
        T2 = 2x * T1 - T0
        dT2 = 2T1 + 2x * dT1 - dT0
        f += c[k] * T2
        df += c[k] * dT2
        T0, T1, dT0, dT1 = T1, T2, dT1, dT2
    end
    return f, df
end

# Components and their time derivatives [per day] of `segs` at `et`
function evaledge(segs::SegmentSchedule, et::Float64; ncomp::Int = 3)
    s = findsegment(segs, et)
    k = recordindex(s, et)
    MID, RADIUS = s.records[1, k], s.records[2, k]
    x = (et - MID) / RADIUS
    vals = zeros(ncomp)
    ders = zeros(ncomp)
    for i in 1:ncomp
        f, df = chebvalder(coeffs(s, k, i), x)
        vals[i] = f
        ders[i] = s.ncomp == 6 ? chebvalder(coeffs(s, k, 3 + i), x)[1] * daysec :
            df * daysec / RADIUS
    end
    return vals, ders
end

# Lunar mantle angular velocity in the mantle frame [rad/day], from the Euler angles
# (φ, θ, ψ) and their time derivatives [rad/day]; inverse of eqs. (14)-(16) of
# Folkner et al. (2014)
function mantleangvel(ang, dang)
    φ, θ, ψ = ang
    dφ, dθ, dψ = dang
    return [dφ * sin(θ) * sin(ψ) + dθ * cos(ψ),
            dφ * sin(θ) * cos(ψ) - dθ * sin(ψ),
            dψ + dφ * cos(θ)]
end

# Lunar mantle angular velocity over [a, b] [TDB s since J2000], as Taylor polynomials
# of degree L-1 about `a` and in TDB days. The angular velocity is not a polynomial in
# time, so it is approximated by Chebyshev interpolation at L nodes (near-minimax),
# and then converted to the Taylor basis
function expandangvel(lunarsegs::SegmentSchedule, a::Float64, b::Float64, L::Int,
                      ws::Workspace)
    nodes = [cos(π * (k + 0.5) / L) for k in 0:L-1]
    W = [mantleangvel(evaledge(lunarsegs, (a + b) / 2 + x * (b - a) / 2)...)
         for x in nodes]
    hdays = (b - a) / daysec
    ω = Vector{Vector{Float64}}(undef, 3)
    for i in 1:3
        c = [2 / L * sum(W[k+1][i] * cos(π * j * (k + 0.5) / L) for k in 0:L-1)
             for j in 0:L-1]
        c[1] /= 2
        ω[i] = cheb2taylor!(zeros(L), ws.T0, ws.T1, ws.T2, c, -1.0, 2 / hdays)
    end
    return ω
end

## Main

function main()
    parsed_args = parse_commandline()

    printitle("Convert JPL SPICE kernels into a PlanetaryEphemeris.jl .jld2 file", "=")
    global_initial_time = now()
    println("• Run started at ", global_initial_time)
    println("• Detected 1 worker with ", Threads.nthreads(), " thread(s)")

    # Read kernels
    bspfiles::Vector{String} = parsed_args["bsp"]
    techfile::String = parsed_args["tech"]
    pckfile::Union{Nothing, String} = parsed_args["pck"]
    output::String = parsed_args["output"]

    spksegs = Dict{Tuple{Int, Int}, Vector{ChebyshevSegment}}()
    skipped = Tuple{Int, Int, Int}[]
    priority = 0
    for file in bspfiles
        kind, segs, skip = readdaf(file, priority)
        kind == :spk || error("$file is not an SPK file")
        priority += length(segs) + length(skip)
        for s in segs
            push!(get!(spksegs, edge(s), ChebyshevSegment[]), s)
        end
        append!(skipped, skip)
        println("• SPK kernel: ", file, " (", length(segs), " segments)")
    end
    for v in values(spksegs)
        sort!(v, by = s -> s.priority)
    end
    edges = Dict{Int, Vector{Int}}()
    for (t, c) in keys(spksegs)
        push!(get!(edges, t, Int[]), c)
    end
    for (t, c, type) in unique(skipped)
        @warn "Skipped unsupported segment ($t wrt $c, SPK type $type)"
    end

    lunarsegs = ChebyshevSegment[]
    if isnothing(pckfile)
        @warn "No lunar PCK given; lunar Euler angles will be filled with zeros"
    else
        kind, segs, skip = readdaf(pckfile, 0)
        kind == :pck || error("$pckfile is not a binary PCK file")
        ids = unique(s.target for s in segs if 31000 ≤ s.target < 32000)
        length(ids) == 1 || error("Expected exactly one lunar frame in $pckfile, \
            found $ids")
        lunarsegs = sort!(filter(s -> s.target == ids[1], segs), by = s -> s.priority)
        all(s -> s.center == J2000_FRAME, lunarsegs) || error("Lunar PCK angles must be \
            referred to the J2000 (ICRF) frame")
        println("• Lunar PCK: ", pckfile, " (frame ", ids[1], ", ",
            length(lunarsegs), " segments)")
    end

    # Technical comments
    text = read(techfile, String)
    println("• Technical comments file: ", techfile)

    # Bodies
    asteroids = if isnothing(parsed_args["asteroids"])
        ids = sort!(unique(t for (t, _) in keys(spksegs) if isasteroid(t)))
        sort!(ids, by = id -> (-techconstant(text, gmname(id)), id))
    else
        [parse(Int, strip(x)) for x in split(parsed_args["asteroids"], ',')]
    end
    allunique(asteroids) || error("Repeated asteroid ids: $asteroids")
    all(isasteroid, asteroids) || error("Invalid asteroid ids: $asteroids")
    spkids = vcat(PLANETS_SPKIDS, asteroids)
    N = length(spkids)
    μ = [techconstant(text, gmname(id)) for id in spkids]

    # Ephemeris model and model parameters
    model::String = isnothing(parsed_args["model"]) ? modelname(techfile, text) :
        parsed_args["model"]
    params = Dict{String, Float64}(name => techconstant(text, name)
        for name in MODEL_PARAMS)
    println("• Ephemeris model: ", model)
    println("• Model parameters: ", join(["$name = $(params[name])"
        for name in MODEL_PARAMS], ", "))

    # Chains of SPK edges
    chains = [spkchain(edges, id) for id in spkids]
    haskey(spksegs, (TTMTDB_SPKID, TTMTDB_CENTER)) || error("The kernels have no \
        TT-TDB segment ($TTMTDB_SPKID wrt $TTMTDB_CENTER); include e.g. de440t.bsp")
    usededges = unique(vcat(chains..., [(TTMTDB_SPKID, TTMTDB_CENTER)]))
    for e in usededges, s in spksegs[e]
        (e[1] == TTMTDB_SPKID || s.frame == J2000_FRAME) || error("Segment $e of \
            $(s.file) is not referred to the J2000 (ICRF) frame (frame = $(s.frame))")
    end

    println("• Included bodies (", N, "):")
    @printf("    %-10s %-24s %-8s %-22s %s\n", "SPK id", "μ [au³/day²]", "GM name",
        "chain", "segments (type, record length [days], degree)")
    for (i, id) in enumerate(spkids)
        segsinfo = join([@sprintf("(%d, %g, %d)", s.type, s.intlen / daysec, degree(s))
            for e in chains[i] for s in spksegs[e]], " ")
        chainstr = join([string(e[1]) for e in chains[i]], "→") * "→0"
        @printf("    %-10d %-24.16e %-8s %-22s %s\n", id, μ[i], gmname(id), chainstr,
            segsinfo)
    end

    # Common coverage = intersection of the coverage of every edge used
    cov = [(-Inf, Inf)]
    for e in usededges
        cov = intersectcoverage(cov, coverage(spksegs[e]))
    end
    isempty(lunarsegs) || (cov = intersectcoverage(cov, coverage(lunarsegs)))
    isempty(cov) && error("The kernels have no common coverage")
    println("• Common coverage [TDB]: ", join(["$(et2datetime(a)) - \
        $(et2datetime(b))" for (a, b) in cov], ", "))

    d0::Union{Nothing, DateTime} = parsed_args["start"]
    df::Union{Nothing, DateTime} = parsed_args["end"]
    if isnothing(d0) && isnothing(df)
        length(cov) == 1 || error("The common coverage has gaps; use --start and \
            --end to select an interval without gaps")
        w0, wf = cov[1]
    else
        w0 = isnothing(d0) ? -Inf : datetime2et(d0)
        wf = isnothing(df) ? Inf : datetime2et(df)
        if isinf(w0) || isinf(wf)
            iv = findfirst(x -> (isinf(w0) || x[1] ≤ w0 ≤ x[2]) &&
                (isinf(wf) || x[1] ≤ wf ≤ x[2]), cov)
            isnothing(iv) && error("The requested interval is outside the common \
                coverage")
            w0 = isinf(w0) ? cov[iv][1] : w0
            wf = isinf(wf) ? cov[iv][2] : wf
        end
        w0 < wf || error("--start must be earlier than --end")
        any(x -> x[1] ≤ w0 && wf ≤ x[2], cov) || error("The requested interval \
            [$(et2datetime(w0)), $(et2datetime(wf))] is not contained in the common \
            coverage")
    end
    println("• Ephemeris interval [TDB]: ", et2datetime(w0), " - ", et2datetime(wf))

    # Time grid: union of the record and segment boundaries of every segment used
    bps = [w0, wf]
    for segs in vcat([spksegs[e] for e in usededges], [lunarsegs])
        for s in segs
            a, b = max(s.et0, w0), min(s.etf, wf)
            a < b || continue
            push!(bps, a, b)
            for k in floor(Int, (a - s.init) / s.intlen):ceil(Int, (b - s.init) / s.intlen)
                t = s.init + k * s.intlen
                a < t < b && push!(bps, t)
            end
        end
    end
    sort!(bps)
    grid = [bps[1]]
    for t in Iterators.drop(bps, 1)
        t - grid[end] > BREAKPOINT_TOL && push!(grid, t)
    end
    grid[end] = wf
    nint = length(grid) - 1
    steps = diff(grid) ./ daysec
    @printf("• Time grid: %d intervals (min/max length: %.6g/%.6g days)\n", nint,
        minimum(steps), maximum(steps))

    # Order of the Taylor polynomials
    usedsegs = vcat([spksegs[e] for e in usededges]..., lunarsegs)
    maxdeg = maximum(degree, usedsegs)
    order::Int = parsed_args["order"]
    order = order < 0 ? maxdeg : order
    order ≥ maxdeg || error("--order ($order) must be ≥ the highest Chebyshev degree \
        in the kernels ($maxdeg)")
    println("• Highest Chebyshev degree in the kernels: ", maxdeg)
    println("• Order of the Taylor polynomials: ", order)
    L = order + 1
    ncols = 6N + 13
    @printf("• Estimated size of the coefficients: %.3f GB\n", nint * ncols * L * 8 / 1E9)

    # Taylor polynomials
    schedules = Dict(e => SegmentSchedule(spksegs[e]) for e in usededges)
    lunarsched = isempty(lunarsegs) ? nothing : SegmentSchedule(lunarsegs)
    P = Matrix{Taylor1{Float64}}(undef, nint, ncols)
    zeropoly() = Taylor1(zeros(L), order)
    Threads.@threads for j in 1:nint
        a, b = grid[j], grid[j+1]
        ws = Workspace(L)
        cache = Dict{Tuple{Int, Int}, NTuple{2, Vector{Vector{Float64}}}}()
        for (i, chain) in enumerate(chains)
            x = [zeros(L) for _ in 1:3]
            v = [zeros(L) for _ in 1:3]
            for e in chain
                xe, ve = get!(() -> expandedge(schedules[e], a, b, L, ws), cache, e)
                for c in 1:3
                    x[c] .+= xe[c]
                    v[c] .+= ve[c]
                end
            end
            idxs = nbodyind(N, i)
            for c in 1:3
                P[j, idxs[c]] = Taylor1(x[c] ./ au, order)
                P[j, idxs[3+c]] = Taylor1(v[c] ./ au, order)
            end
        end
        # Lunar mantle Euler angles and angular velocity (in the mantle frame)
        if isempty(lunarsegs)
            for c in 6N+1:6N+6
                P[j, c] = zeropoly()
            end
        else
            ang, _ = expandedge(lunarsched, a, b, L, ws)
            ω = expandangvel(lunarsched, a, b, L, ws)
            for c in 1:3
                P[j, 6N+c] = Taylor1(ang[c], order)
                P[j, 6N+3+c] = Taylor1(ω[c], order)
            end
        end
        # Lunar core (not available in SPICE kernels)
        for c in 6N+7:6N+12
            P[j, c] = zeropoly()
        end
        # TT-TDB [s]
        tt, _ = expandedge(schedules[(TTMTDB_SPKID, TTMTDB_CENTER)], a, b, L, ws;
            ncomp = 1)
        P[j, 6N+13] = Taylor1(tt[1], order)
    end
    sseph = TaylorSolution(grid ./ daysec, P)

    # Validation against direct evaluation of the Chebyshev records
    ncheck::Int = parsed_args["check"]
    if ncheck > 0
        errx, errxrel, errv, errtt, errang, errω = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
        for k in 1:ncheck
            # Deterministic, well spread epochs; dyadic fractions of each interval keep
            # the epochs exactly representable both in days and in seconds
            j = clamp(round(Int, k / ncheck * nint + 0.5), 1, nint)
            u = clamp(round(Int, 1024 * mod(k * 0.6180339887498949, 1)), 1, 1023) / 1024
            tdays = sseph.t[j] + u * (sseph.t[j+1] - sseph.t[j])
            et = tdays * daysec
            y = sseph(tdays)
            for (i, chain) in enumerate(chains)
                x, v = zeros(3), zeros(3)
                for e in chain
                    xe, ve = evaledge(schedules[e], et)
                    x .+= xe
                    v .+= ve
                end
                idxs = nbodyind(N, i)
                Δx = maximum(abs, y[idxs[1:3]] .- x ./ au)
                errx = max(errx, Δx)
                errxrel = max(errxrel, Δx / (maximum(abs, x) / au))
                errv = max(errv, maximum(abs, y[idxs[4:6]] .- v ./ au))
            end
            tt, _ = evaledge(schedules[(TTMTDB_SPKID, TTMTDB_CENTER)], et;
                ncomp = 1)
            errtt = max(errtt, abs(y[6N+13] - tt[1]))
            if !isempty(lunarsegs)
                ang, dang = evaledge(lunarsched, et)
                errang = max(errang, maximum(abs, y[6N+1:6N+3] .- ang))
                errω = max(errω, maximum(abs, y[6N+4:6N+6] .- mantleangvel(ang, dang)))
            end
        end
        println("• Maximum differences w.r.t. direct Chebyshev evaluation at ", ncheck,
            " epochs:")
        @printf("    positions: %.3e au (%.3e km); relative: %.3e\n", errx, errx * au,
            errxrel)
        @printf("    velocities: %.3e au/day (%.3e km/s)\n", errv, errv * au / daysec)
        @printf("    TT-TDB: %.3e s\n", errtt)
        if !isempty(lunarsegs)
            @printf("    lunar Euler angles: %.3e rad\n", errang)
            @printf("    lunar angular velocity: %.3e rad/day\n", errω)
        end
    end

    # Save output
    kernels = basename.(vcat(bspfiles, isnothing(pckfile) ? String[] : [pckfile]))
    jldsave(output; ss16ast_eph = sseph, spkids, μ, model, params, kernels)
    println("• Saved output to: ", output)

    global_final_time = now()
    println("• Run started ", global_initial_time, " and finished ", global_final_time)
    println("• Total computation time was: ",
        computationtime(global_initial_time, global_final_time), " min")

    return nothing
end

main()
