# Headless GR backend (save plots to files without opening a window)
ENV["GKSwstype"] = "100"

using ArgParse, Dates, JLD2, LinearAlgebra, Plots, Printf, SPICE, TaylorIntegration,
      TaylorSeries
using PlanetaryEphemeris
using PlanetaryEphemeris: J2000, auday2kmsec, daysec, nbodyind, numberofbodies, yr

# NAIF ids of the Solar System barycenter, the Sun, the Earth, the Moon and TT-TDB
const SSB_ID = 0
const SUN_ID = 10
const EARTH_ID = 399
const MOON_ID = 301
const TTMTDB_ID = 1000000001
const TTMTDB_CENTER = 1000000000
# NAIF ids of the planets, the Moon and Pluto in PlanetaryEphemeris.jl ephemerides
const PLANETS_IDS = [10, 199, 299, 399, 301, 4, 5, 6, 7, 8, 9]
# Name of the lunar orientation frame defined for the comparison
const PA_FRAME = "MOON_PA_PE_COMPARISON"

function parse_commandline()

    s = ArgParseSettings()

    # Program name (for usage & help screen)
    s.prog = "spice_vs_pe.jl"
    # Desciption (for help screen)
    s.description = "Compare a PlanetaryEphemeris.jl ephemeris file against SPICE kernels"

    @add_arg_table! s begin
        "--ephemeris", "-j"
            help = "PlanetaryEphemeris.jl ephemeris (.jld2) file"
            arg_type = String
            required = true
        "--kernel", "-k"
            help = "SPICE SPK kernel(s) (.bsp); later kernels have higher priority"
            arg_type = String
            nargs = '+'
            required = true
        "--pck", "-p"
            help = "lunar orientation binary PCK (.bpc), e.g. moon_pa_de440_200625.bpc; \
                if omitted, the lunar mantle Euler angles are not compared"
            arg_type = String
        "--ids"
            help = "comma-separated NAIF ids of the bodies in the ephemeris file; only \
                needed if the file does not store them (keys `ids` or `spkids`)"
            arg_type = String
        "--outdir", "-o"
            help = "output directory for the plots"
            arg_type = String
            required = true
        "--start", "-s"
            help = "start of the comparison [TDB]; default: start of the common time span"
            arg_type = DateTime
        "--end", "-e"
            help = "end of the comparison [TDB]; default: end of the common time span"
            arg_type = DateTime
        "--step"
            help = "time step between comparison points [day]"
            arg_type = Float64
            default = 1.0
        "--maxpoints"
            help = "maximum number of points per curve in the plots; longer curves are \
                reduced to their min/max envelope (0 to plot every point)"
            arg_type = Int
            default = 5_000
        "--label", "-l"
            help = "label of the SPICE ephemeris in the plots (e.g. \"JPL DE440\"); \
                default: names of the kernels"
            arg_type = String
    end

    s.epilog = """
        examples:\n
        \n
        julia --project spice_vs_pe.jl -j de440.jld2 -k de440.bsp -o plots_de440\n
        \n
        julia --project spice_vs_pe.jl -j de430.jld2 -k de430t.bsp sb431-n16s.bsp -p moon_pa_de430_1550-2650.bpc -l "JPL DE430" -o plots_de430\n
        \n
    """

    return parse_args(s)
end

printitle(s::AbstractString, d::AbstractString) = println(d ^ length(s),
    '\n', s, '\n', d ^ length(s))

# TDB days since J2000 <-> DateTime
datetime2days(d::DateTime) = datetime2julian(d) - J2000
days2datetime(t::Real) = julian2datetime(J2000 + t)

# Load the ephemeris and the NAIF ids of its bodies from a .jld2 file. Supported
# layouts: `sseph` + `ids` (PlanetaryEphemeris.jl integrations) and `ss16ast_eph` +
# `spkids` (scripts/spice_to_pe.jl). The ids can also be given explicitly (`idsarg`)
function loadpe(filename::AbstractString, idsarg::Union{Nothing, AbstractString})
    file = JLD2.load(filename)
    key = findfirst(k -> haskey(file, k), ["sseph", "ss16ast_eph"])
    isnothing(key) && error("$filename has no ephemeris (keys `sseph` or `ss16ast_eph`)")
    sseph = file[["sseph", "ss16ast_eph"][key]]
    N = numberofbodies(sseph)
    ids = if !isnothing(idsarg)
        split(idsarg, ',')
    elseif haskey(file, "ids")
        file["ids"]
    elseif haskey(file, "spkids")
        file["spkids"]
    elseif isdefined(PlanetaryEphemeris, :de430_343ast_ids)
        @warn "$filename does not store the NAIF ids of its bodies; assuming the \
            planets, the Moon, Pluto and the first $(N - length(PLANETS_IDS)) \
            asteroids of PlanetaryEphemeris.de430_343ast_ids"
        vcat(PLANETS_IDS, 2_000_000 .+
            PlanetaryEphemeris.de430_343ast_ids[1:N-length(PLANETS_IDS)])
    else
        error("$filename does not store the NAIF ids of its bodies; use --ids")
    end
    ids = [id isa Integer ? Int(id) : parse(Int, strip(id)) for id in ids]
    @assert length(ids) == N "Number of bodies in $filename ($N) does not match its \
        number of NAIF ids ($(length(ids)))"
    return sseph, ids
end

# Whether all the polynomials in columns `cols` of `sseph` are identically zero
istrivial(sseph::TaylorSolution, cols) = all(iszero, view(sseph.p, :, cols))

# Name of a body according to SPICE, or its NAIF id if SPICE does not know it
function bodyname(id::Int)
    name = try
        bodc2n(id)
    catch
        nothing
    end
    return isnothing(name) ? string(id) : name
end

# Time span [days since J2000 TDB] covered by `kernels` for all the NAIF ids in `ids`
function kernelspan(kernels::AbstractVector{String}, ids::AbstractVector{Int})
    tmin, tmax = -Inf, Inf
    for id in ids
        cover = SpiceDoubleCell(2_000)
        for kernel in kernels
            spkcov!(cover, kernel, id)
        end
        window = collect(cover)
        tmin = max(tmin, first(window) / daysec)
        tmax = min(tmax, last(window) / daysec)
    end
    return tmin, tmax
end

# Lunar mantle Euler angle rates [rad/day] from the Euler angles `(ϕ, θ, ψ)` [rad] and the
# angular velocity `ω` [rad/day] in the mantle frame
# See equation (14) in page 9 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract
function eulerrates(ϕ, θ, ψ, ω)
    dϕ = (ω[1] * sin(ψ) + ω[2] * cos(ψ)) / sin(θ)
    dθ = ω[1] * cos(ψ) - ω[2] * sin(ψ)
    dψ = ω[3] - dϕ * cos(θ)
    return dϕ, dθ, dψ
end

# Min/max envelope of the columns of `y` (sampled at `x`) for plotting: if there are more
# than `maxpoints` rows, split them into `maxpoints ÷ 2` bins and keep, in each bin, the
# points with the minimum and maximum value of each column (so no extremum is lost).
# Returns matrices `X`, `Y` with one column per column of `y`
function envelope(x::AbstractVector, y::AbstractVecOrMat, maxpoints::Int)
    y = reshape(y, size(y, 1), :)
    n, m = size(y)
    (maxpoints ≤ 0 || n ≤ maxpoints) && return repeat(x, 1, m), Matrix(y)
    nbins = max(1, maxpoints ÷ 2)
    edges = round.(Int, range(0, n, length = nbins + 1))
    X = Matrix{eltype(x)}(undef, 2nbins, m)
    Y = Matrix{eltype(y)}(undef, 2nbins, m)
    for c in 1:m, b in 1:nbins
        idxs = edges[b]+1:edges[b+1]
        imin = idxs[argmin(view(y, idxs, c))]
        imax = idxs[argmax(view(y, idxs, c))]
        i1, i2 = minmax(imin, imax)
        X[2b-1, c], X[2b, c] = x[i1], x[i2]
        Y[2b-1, c], Y[2b, c] = y[i1, c], y[i2, c]
    end
    return X, Y
end

# Save `fig` to `outdir/name` and report it
function savefigure(fig, outdir::AbstractString, name::AbstractString)
    filename = joinpath(outdir, name)
    savefig(fig, filename)
    println("Saved ", filename)
    return nothing
end

# Load the PlanetaryEphemeris.jl ephemeris and the SPICE kernels, and compute the
# differences PlanetaryEphemeris.jl - SPICE every `step` days. Keeping this in its own
# function lets the (large) ephemeris be garbage collected before plotting
function compare(ephfile::String, kernels::Vector{String}, pck::Union{Nothing, String},
                 idsarg::Union{Nothing, String}, step::Float64,
                 d0::Union{Nothing, DateTime}, df::Union{Nothing, DateTime})

    # PlanetaryEphemeris.jl ephemeris
    sseph, ids = loadpe(ephfile, idsarg)
    N = length(ids)
    println("• Number of bodies: ", N)
    println("• NAIF ids: ", join(ids, ", "))

    # Indices of the Sun, the Earth and the Moon in the ephemeris
    isu, iea, imo = (findfirst(==(id), ids) for id in (SUN_ID, EARTH_ID, MOON_ID))
    @assert !any(isnothing, (isu, iea, imo)) "The ephemeris must contain the Sun, \
        the Earth and the Moon"
    # Center of the comparison of each body (index in the ephemeris, 0 for the Solar
    # System barycenter): barycentric for the Sun, geocentric for the Moon and
    # heliocentric for the other bodies
    centers = [i == isu ? 0 : i == imo ? iea : isu for i in 1:N]
    center_ids = [c == 0 ? SSB_ID : ids[c] for c in centers]
    frames = [c == 0 ? "barycentric" : c == iea ? "geocentric" : "heliocentric"
        for c in centers]

    # SPICE kernels
    foreach(furnsh, kernels)
    kernel_ids = SpiceIntCell(10_000)
    foreach(k -> spkobj!(kernel_ids, k), kernels)
    kernel_ids = collect(kernel_ids)
    missing_ids = setdiff(ids, kernel_ids)
    @assert isempty(missing_ids) "Bodies not in the SPICE kernels: \
        $(join(missing_ids, ", "))"

    # TT-TDB is compared only if both the kernels and the ephemeris have it
    compare_ttmtdb = if TTMTDB_ID ∉ kernel_ids
        println("• TT-TDB: not compared (not in the SPICE kernels)")
        false
    elseif istrivial(sseph, 6N+13)
        println("• TT-TDB: not compared (identically zero in the ephemeris)")
        false
    else
        println("• TT-TDB: compared")
        true
    end

    # The lunar mantle Euler angles are compared only if a lunar binary PCK is given
    # and the ephemeris has non trivial lunar mantle Euler angles
    compare_eulang = if isnothing(pck)
        println("• Lunar mantle Euler angles: not compared (no lunar binary PCK)")
        false
    elseif istrivial(sseph, 6N+1:6N+6)
        println("• Lunar mantle Euler angles: not compared (identically zero in the \
            ephemeris)")
        false
    else
        println("• Lunar mantle Euler angles: compared")
        true
    end

    # Lunar orientation binary PCK: define a PCK-based frame (centered at the Moon) for
    # the frame class ID in the file, since no frame kernel is loaded
    pa_id = 0
    if compare_eulang
        furnsh(pck)
        pck_ids = SpiceIntCell(10)
        pckfrm!(pck_ids, pck)
        @assert length(pck_ids) == 1 "$pck must contain exactly one frame class ID"
        pa_id = pck_ids[1]
        pipool("FRAME_$PA_FRAME", [pa_id])
        pcpool("FRAME_$(pa_id)_NAME", [PA_FRAME])
        pipool("FRAME_$(pa_id)_CLASS", [2])
        pipool("FRAME_$(pa_id)_CLASS_ID", [pa_id])
        pipool("FRAME_$(pa_id)_CENTER", [MOON_ID])
        println("• Lunar orientation frame class ID: ", pa_id)
    end

    # Time span covered by both PlanetaryEphemeris.jl and SPICE
    pe_tmin, pe_tmax = extrema(sseph.t)
    jpl_tmin, jpl_tmax = kernelspan(kernels,
        compare_ttmtdb ? vcat(ids, TTMTDB_ID) : ids)
    if compare_eulang
        pck_cover = SpiceDoubleCell(2_000)
        pckcov!(pck_cover, pck, pa_id)
        pck_window = collect(pck_cover)
        jpl_tmin = max(jpl_tmin, first(pck_window) / daysec)
        jpl_tmax = min(jpl_tmax, last(pck_window) / daysec)
    end
    tmin, tmax = max(pe_tmin, jpl_tmin), min(pe_tmax, jpl_tmax)
    isnothing(d0) || (tmin = max(tmin, datetime2days(d0)))
    isnothing(df) || (tmax = min(tmax, datetime2days(df)))
    @assert tmin < tmax "PlanetaryEphemeris.jl and SPICE time spans do not overlap"
    println("• PlanetaryEphemeris.jl time span: ", days2datetime(pe_tmin), " - ",
        days2datetime(pe_tmax))
    println("• SPICE time span: ", days2datetime(jpl_tmin), " - ",
        days2datetime(jpl_tmax))
    println("• Comparison time span: ", days2datetime(tmin), " - ", days2datetime(tmax))

    # J2000 state [km, km/s] of the `i`-th body relative to the `c`-th body (`c = 0`
    # for the Solar System barycenter) in the full PlanetaryEphemeris.jl state `x`
    function pe_ephemeris(i::Int, c::Int, x)
        xi = x[nbodyind(N, i)]
        xc = c == 0 ? zero(xi) : x[nbodyind(N, c)]
        return auday2kmsec(xi .- xc)
    end
    # J2000 state [km, km/s] of the body with NAIF id `id` relative to the body with
    # NAIF id `c` in the SPICE kernels at time `t` [days since J2000 TDB]
    jpl_ephemeris(id::Int, c::Int, t) = spkgeo(id, t * daysec, "J2000", c)[1]
    # TT-TDB [s] in the SPICE kernels at time `t` [days since J2000 TDB]
    jpl_ttmtdb(t) = spkgeo(TTMTDB_ID, t * daysec, "J2000", TTMTDB_CENTER)[1][1]
    # Lunar mantle Euler angles [rad] and their rates [rad/day] in the full
    # PlanetaryEphemeris.jl state `x`
    function pe_eulang(x)
        ϕ, θ, ψ = x[6N+1], x[6N+2], x[6N+3]
        return [ϕ, θ, ψ, eulerrates(ϕ, θ, ψ, x[6N+4:6N+6])...]
    end
    # Lunar mantle Euler angles [rad] and their rates [rad/day] in the SPICE binary PCK
    # at time `t` [days since J2000 TDB]. The rotation from J2000 to the mantle frame is
    # R = R_z(ψ) R_x(θ) R_z(ϕ) and dR/dt = -[ω×] R, with ω the angular velocity in the
    # mantle frame
    function jpl_eulang(t)
        S = sxform("J2000", PA_FRAME, t * daysec)
        R, dR = S[1:3, 1:3], S[4:6, 1:3]
        ψ, θ, ϕ = m2eul(R, 3, 1, 3)
        W = -(dR * transpose(R)) .* daysec
        ω = [W[3, 2], W[1, 3], W[2, 1]]
        return [ϕ, θ, ψ, eulerrates(ϕ, θ, ψ, ω)...]
    end

    # Differences PlanetaryEphemeris.jl - SPICE
    printitle("Differences PlanetaryEphemeris.jl - SPICE", "-")
    ts = collect(tmin:step:tmax)
    last(ts) < tmax && push!(ts, tmax)
    # States [km, km/s]
    Δ = Array{Float64}(undef, length(ts), N, 6)
    # Moon's geocentric range [km]
    Δρ = Vector{Float64}(undef, length(ts))
    # TT-TDB [s]
    Δttmtdb = Vector{Float64}(undef, compare_ttmtdb ? length(ts) : 0)
    # Lunar mantle Euler angles [rad] and their rates [rad/day]
    Δeulang = Array{Float64}(undef, compare_eulang ? length(ts) : 0, 6)
    @time for (j, t) in enumerate(ts)
        x = sseph(t)
        for i in 1:N
            Δ[j, i, :] .= pe_ephemeris(i, centers[i], x) .-
                jpl_ephemeris(ids[i], center_ids[i], t)
        end
        Δρ[j] = norm(pe_ephemeris(imo, iea, x)[1:3]) -
            norm(jpl_ephemeris(MOON_ID, EARTH_ID, t)[1:3])
        compare_ttmtdb && (Δttmtdb[j] = x[6N+13] - jpl_ttmtdb(t))
        compare_eulang && (Δeulang[j, :] .= pe_eulang(x) .- jpl_eulang(t))
    end
    # Euler angles are only defined modulo 2π
    Δeulang[:, 1:3] .= rem2pi.(Δeulang[:, 1:3], RoundNearest)
    # Positions and ranges in m, velocities in mm/s, TT-TDB in ns and Euler angles and
    # their rates in mas and mas/day
    Δ[:, :, 1:3] .*= 1E3
    Δ[:, :, 4:6] .*= 1E6
    Δρ .*= 1E3
    Δttmtdb .*= 1E9
    Δeulang .*= rad2deg(1) * 3600 * 1E3

    return (; ids, frames, ts, Δ, Δρ, Δttmtdb, Δeulang, compare_ttmtdb, compare_eulang)
end

function main()

    printitle("Compare a PlanetaryEphemeris.jl ephemeris against SPICE kernels", "=")

    # Parse arguments from commandline
    parsed_args = parse_commandline()

    # PlanetaryEphemeris.jl ephemeris file
    ephfile::String = parsed_args["ephemeris"]
    println("• PlanetaryEphemeris.jl ephemeris file: ", ephfile)

    # SPICE kernels
    kernels::Vector{String} = parsed_args["kernel"]
    println("• SPICE kernel(s): ", join(kernels, ", "))

    # Lunar orientation binary PCK
    pck::Union{Nothing, String} = parsed_args["pck"]
    println("• Lunar orientation binary PCK: ", something(pck, "none"))

    # Output directory
    outdir::String = parsed_args["outdir"]
    println("• Output directory: ", outdir)

    # Time step
    step::Float64 = parsed_args["step"]
    println("• Time step [day]: ", step)

    # Maximum number of points per curve in the plots
    maxpoints::Int = parsed_args["maxpoints"]

    # Label of the SPICE ephemeris
    label::String = something(parsed_args["label"], join(basename.(kernels), " + "))
    println("• Label of the SPICE ephemeris: ", label)

    # Differences PlanetaryEphemeris.jl - SPICE
    (; ids, frames, ts, Δ, Δρ, Δttmtdb, Δeulang, compare_ttmtdb, compare_eulang) =
        compare(ephfile, kernels, pck, parsed_args["ids"], step, parsed_args["start"],
            parsed_args["end"])
    N = length(ids)
    # Free the memory used by the ephemeris
    GC.gc()

    # Table of maximum absolute differences
    printitle("Maximum absolute differences PlanetaryEphemeris.jl - SPICE", "-")
    @printf("%5s  %-10s %-20s %-12s %12s %12s %12s %12s %12s %12s\n", "i", "NAIF id",
        "body", "frame", "x [m]", "y [m]", "z [m]", "vx [mm/s]", "vy [mm/s]",
        "vz [mm/s]")
    for i in 1:N
        @printf("%5d  %-10d %-20s %-12s", i, ids[i], bodyname(ids[i]), frames[i])
        for k in 1:6
            @printf(" %12.4e", maximum(abs, view(Δ, :, i, k)))
        end
        println()
    end
    println()
    @printf("Moon geocentric range [m]: %12.4e\n", maximum(abs, Δρ))
    compare_ttmtdb && @printf("TT-TDB [ns]:               %12.4e\n",
        maximum(abs, Δttmtdb))
    if compare_eulang
        println()
        @printf("%-20s %12s %12s %12s %12s %12s %12s\n", "Lunar mantle", "ϕ [mas]",
            "θ [mas]", "ψ [mas]", "dϕ [mas/d]", "dθ [mas/d]", "dψ [mas/d]")
        @printf("%-20s", "Euler angles")
        for k in 1:6
            @printf(" %12.4e", maximum(abs, view(Δeulang, :, k)))
        end
        println()
    end

    # Plots
    printitle("Plots", "-")
    mkpath(outdir)
    years = 2000 .+ ts ./ yr
    title_suffix = "PlanetaryEphemeris.jl - $label"
    # Common plot attributes
    single = (size = (900, 400), left_margin = 5Plots.mm, bottom_margin = 5Plots.mm,
        dpi = 300, label = false, xlabel = "Year (TDB)")
    double = (layout = (2, 1), size = (900, 700), left_margin = 5Plots.mm, dpi = 300)
    for i in 1:N
        name = bodyname(ids[i])
        p1 = plot(envelope(years, Δ[:, i, 1:3], maxpoints)...; label = ["x" "y" "z"],
            ylabel = "Δ position [m]",
            title = "$name ($(ids[i])), $(frames[i]): $title_suffix")
        p2 = plot(envelope(years, Δ[:, i, 4:6], maxpoints)...;
            label = ["vx" "vy" "vz"], ylabel = "Δ velocity [mm/s]", xlabel = "Year (TDB)")
        fig = plot(p1, p2; double...)
        savefigure(fig, outdir, @sprintf("%03d_%d.png", i, ids[i]))
    end
    # Moon's geocentric range
    fig = plot(envelope(years, Δρ, maxpoints)...; ylabel = "Δ range [m]",
        title = "Moon geocentric range: $title_suffix", single...)
    savefigure(fig, outdir, "moon_range.png")
    # TT-TDB
    if compare_ttmtdb
        fig = plot(envelope(years, Δttmtdb, maxpoints)...; ylabel = "Δ(TT-TDB) [ns]",
            title = "TT-TDB: $title_suffix", single...)
        savefigure(fig, outdir, "ttmtdb.png")
    end
    # Lunar mantle Euler angles and their rates
    if compare_eulang
        p1 = plot(envelope(years, Δeulang[:, 1:3], maxpoints)...;
            label = ["ϕ" "θ" "ψ"], ylabel = "Δ angle [mas]",
            title = "Lunar mantle Euler angles: $title_suffix")
        p2 = plot(envelope(years, Δeulang[:, 4:6], maxpoints)...;
            label = ["dϕ/dt" "dθ/dt" "dψ/dt"], ylabel = "Δ rate [mas/day]",
            xlabel = "Year (TDB)")
        fig = plot(p1, p2; double...)
        savefigure(fig, outdir, "moon_eulang.png")
    end

    return nothing
end

main()
