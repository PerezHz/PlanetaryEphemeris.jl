# Headless GR backend (save plots to files without opening a window)
ENV["GKSwstype"] = "100"

using ArgParse, PlanetaryEphemeris, Dates, JLD2, Plots, Printf, SPICE

function parse_commandline()

    s = ArgParseSettings()

    # Program name (for usage & help screen)
    s.prog = "jpl_vs_pe_de440.jl"
    # Desciption (for help screen)
    s.description = "Compare a PlanetaryEphemeris.jl integration of the JPL DE440 model \
        against a JPL DE440 SPICE kernel"

    @add_arg_table! s begin
        "--bwd", "-b"
            help = "Backward integration (.jld2) ephemeris file"
            arg_type = String
            required = true
        "--fwd", "-f"
            help = "Forward integration (.jld2) ephemeris file"
            arg_type = String
            required = true
        "--kernel", "-k"
            help = "JPL DE440 SPICE kernel (.bsp) file"
            arg_type = String
            required = true
        "--outdir", "-o"
            help = "Output directory for the plots"
            arg_type = String
            required = true
        "--step"
            help = "Time step between comparison points [day]"
            arg_type = Float64
            default = 1.0
    end

    s.epilog = """
        examples:\n
        \n
        julia --project jpl_vs_pe_de440.jl -b de440_bwd.jld2 -f de440_fwd.jld2 -k de440.bsp -o plots_de440\n
        \n
    """

    return parse_args(s)
end

printitle(s::AbstractString, d::AbstractString) = println(d ^ length(s),
    '\n', s, '\n', d ^ length(s))

# Load the ephemeris `sseph` and the NAIF ids `ids` of its bodies from a .jld2 file
function loadpe(filename::AbstractString)
    sseph, ids = JLD2.load(filename, "sseph", "ids")
    @assert numberofbodies(sseph) == length(ids) "Number of bodies in $filename \
        does not match its number of NAIF ids"
    return sseph, parse.(Int, ids)
end

# Name of a body according to SPICE, or its NAIF id if SPICE does not know it
function bodyname(id::Int)
    name = try
        bodc2n(id)
    catch
        nothing
    end
    return isnothing(name) ? string(id) : name
end

function main()

    printitle("Compare JPL DE440 and PlanetaryEphemeris.jl", "=")

    # Parse arguments from commandline
    parsed_args = parse_commandline()

    # Backward integration file
    bwdfile::String = parsed_args["bwd"]
    println("• Backward integration file: ", bwdfile)

    # Forward integration file
    fwdfile::String = parsed_args["fwd"]
    println("• Forward integration file: ", fwdfile)

    # SPICE kernel
    kernel::String = parsed_args["kernel"]
    println("• SPICE kernel: ", kernel)

    # Output directory
    outdir::String = parsed_args["outdir"]
    println("• Output directory: ", outdir)

    # Time step
    step::Float64 = parsed_args["step"]
    println("• Time step [day]: ", step)

    # PlanetaryEphemeris.jl ephemerides
    sseph_bwd, ids = loadpe(bwdfile)
    sseph_fwd, ids_fwd = loadpe(fwdfile)
    @assert ids == ids_fwd "Backward and forward integrations contain different bodies"
    N = length(ids)
    # Epoch [days since J2000 TDB] and time span
    t0 = sseph_fwd.t[1]
    @assert sseph_bwd.t[1] == t0 "Backward and forward integrations start at \
        different epochs"
    tmin, tmax = minimum(sseph_bwd.t), maximum(sseph_fwd.t)
    println("• Epoch: ", julian2datetime(J2000 + t0))
    println("• Time span: ", julian2datetime(J2000 + tmin), " - ",
        julian2datetime(J2000 + tmax))
    println("• Number of bodies: ", N)

    # SPICE kernel
    furnsh(kernel)
    kernel_ids = SpiceIntCell(10_000)
    spkobj!(kernel_ids, kernel)
    missing_ids = setdiff(ids, collect(kernel_ids))
    @assert isempty(missing_ids) "Bodies not in $kernel: $(join(missing_ids, ", "))"

    # Full PlanetaryEphemeris.jl state [au, au/day] at time `t` [days since J2000 TDB]
    pe_state(t) = t ≥ t0 ? sseph_fwd(t) : sseph_bwd(t)
    # Barycentric J2000 state [km, km/s] of the `i`-th body in the
    # PlanetaryEphemeris.jl ephemeris at time `t` [days since J2000 TDB]; the full
    # state `x` at `t` can be passed to avoid evaluating it once per body
    pe_ephemeris(i::Int, t, x = pe_state(t)) = auday2kmsec(x[nbodyind(N, i)])
    # Barycentric J2000 state [km, km/s] of the body with NAIF id `id` in the
    # SPICE kernel at time `t` [days since J2000 TDB]
    jpl_ephemeris(id::Int, t) = spkgeo(id, t * daysec, "J2000", 0)[1]

    # Differences PlanetaryEphemeris.jl - JPL [km, km/s]
    printitle("Differences PlanetaryEphemeris.jl - JPL", "-")
    ts = collect(tmin:step:tmax)
    Δ = Array{Float64}(undef, length(ts), N, 6)
    @time for (j, t) in enumerate(ts)
        x = pe_state(t)
        for i in 1:N
            Δ[j, i, :] .= pe_ephemeris(i, t, x) .- jpl_ephemeris(ids[i], t)
        end
    end
    # Positions in m and velocities in mm/s
    Δ[:, :, 1:3] .*= 1E3
    Δ[:, :, 4:6] .*= 1E6

    # Table of maximum absolute differences
    printitle("Maximum absolute differences PlanetaryEphemeris.jl - JPL", "-")
    @printf("%5s  %-20s %12s %12s %12s %12s %12s %12s\n", "i", "body", "x [m]",
        "y [m]", "z [m]", "vx [mm/s]", "vy [mm/s]", "vz [mm/s]")
    for i in 1:N
        @printf("%5d  %-20s", i, bodyname(ids[i]))
        for k in 1:6
            @printf(" %12.4e", maximum(abs, view(Δ, :, i, k)))
        end
        println()
    end

    # Plots
    printitle("Plots", "-")
    mkpath(outdir)
    years = 2000 .+ ts ./ yr
    for i in 1:N
        name = bodyname(ids[i])
        p1 = plot(years, Δ[:, i, 1:3], label = ["x" "y" "z"], ylabel = "Δ position [m]",
            title = "$name ($(ids[i])): PlanetaryEphemeris.jl - JPL DE440")
        p2 = plot(years, Δ[:, i, 4:6], label = ["vx" "vy" "vz"],
            ylabel = "Δ velocity [mm/s]", xlabel = "Year (TDB)")
        fig = plot(p1, p2, layout = (2, 1), size = (900, 700), left_margin = 5Plots.mm)
        filename = joinpath(outdir, @sprintf("%03d_%d.png", i, ids[i]))
        savefig(fig, filename)
        println("Saved ", filename)
    end

    return nothing
end

main()
