using ArgParse, PlanetaryEphemeris, Dates
using Downloads, SPICE, Printf

const TECH_COMMENTS_URL = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de430_tech-comments.txt"
const LEAP_SECONDS_KERNEL_URL = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/naif0012.tls"
const PLANETS_KERNEL_URL = "https://ssd.jpl.nasa.gov/ftp/eph/planets/bsp/de430_1850-2150.bsp"
const ASTEROIDS_KERNEL_URL = "https://ssd.jpl.nasa.gov/ftp/eph/small_bodies/asteroids_de430/ast343de430.bsp"
const TTMTDB_KERNEL_URL = "https://ssd.jpl.nasa.gov/ftp/eph/planets/bsp/TTmTDB.de430.19feb2015.bsp"

function parse_commandline()

    s = ArgParseSettings()

    # Program name (for usage & help screen)
    s.prog = "de430ic.jl"
    # Desciption (for help screen)
    s.description = "Generate the initial conditions file for the DE430 ephemeris"

    @add_arg_table! s begin
        "--output", "-o"
            help = "Output file"
            arg_type = String
            default = "de430ic.txt"
    end

    s.epilog = """
        examples:\n
        \n
        julia --project de430ic.jl -o de430ic.txt\n
        \n
    """

    return parse_args(s)
end

computationtime(x::DateTime, y::DateTime) = @sprintf("%.2f", (y - x).value / 60_000)

printitle(s::AbstractString, d::AbstractString) = println(d ^ length(s),
    '\n', s, '\n', d ^ length(s))

# Coordinates and units of the positions and velocities of each body
const POSVEL_COORDS = ("x", "y", "z", "vx", "vy", "vz")
const POSVEL_UNITS = ("au", "au", "au", "au/day", "au/day", "au/day")

# Coordinates and units of the lunar mantle and core Euler angles and angular velocities
const LUNAR_COORDS = ("phi", "theta", "psi", "omega_x", "omega_y", "omega_z")
const LUNAR_UNITS = ("rad", "rad", "rad", "rad/day", "rad/day", "rad/day")

# Line of an initial conditions file
icline(index, id, coord, value, unit) = string(lpad(index, 7), "  ", rpad(id, 13),
    rpad(coord, 12), rpad(@sprintf("%+.16E", value), 26), unit)

"""
    write_initial_conditions(filename, model, jd0, sources, ids, rvs, extras)

Write the initial conditions of a planetary ephemeris `model` at epoch `jd0` [JD TDB],
with one coordinate per line in the order used by PlanetaryEphemeris: the positions
of all bodies, the velocities of all bodies and the additional coordinates. `ids`
are the bodies' ids, `rvs` their positions [au] and velocities [au/day] (one row per
body) and `extras` a vector of `(id, coordinate, value, unit)` tuples.
"""
function write_initial_conditions(filename::AbstractString, model::AbstractString,
        jd0::Real, sources::AbstractVector{<:AbstractString},
        ids::AbstractVector{<:AbstractString}, rvs::AbstractMatrix{<:Real},
        extras::AbstractVector)
    N = length(ids)
    ncoords = 6N + length(extras)
    open(filename, "w") do io
        println(io, "# format   = PlanetaryEphemeris.jl initial conditions v1")
        println(io, "# model    = ", model)
        println(io, "# epoch    = ", jd0, "  # JD TDB (", julian2datetime(jd0), ")")
        println(io, "# frame    = ICRF")
        println(io, "# nbodies  = ", N)
        println(io, "# ncoords  = ", ncoords)
        println(io, "# sources  = ", join(sources, ", "))
        println(io, "# created  = ", round(now(), Second), " by scripts/",
            basename(@__FILE__))
        println(io, "#")
        println(io, "#", lpad("index", 6), "  ", rpad("id", 13), rpad("coordinate", 12),
            rpad("value", 26), "unit")
        index = 0
        # Positions (k = 1:3) and velocities (k = 4:6) of all bodies
        for ks in (1:3, 4:6), i in 1:N, k in ks
            index += 1
            println(io, icline(index, ids[i], POSVEL_COORDS[k], rvs[i, k], POSVEL_UNITS[k]))
        end
        # Additional coordinates
        for (id, coord, value, unit) in extras
            index += 1
            println(io, icline(index, id, coord, value, unit))
        end
    end
    return nothing
end

float_regex(c::Char = 'E') =  Regex("([\\s|\\+|\\-]\\d\\.\\d+$c[\\+|\\-]\\d{2})")

function parse_constant(name::AbstractString, text::AbstractString; c::Char = 'E')
    re = Regex("$name\\s+") * float_regex(c)
    m = match(re, text)
    x = first(m)
    if c != 'E'
        x = replace(x, c => 'E')
    end
    return parse(Float64, x)
end

function planet_id2gm(id::Int)
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

function main()

    # Parse arguments from commandline
    parsed_args = parse_commandline()
    # Output file
    output::String = parsed_args["output"]

    printitle("Generate the initial conditions file for the DE430 ephemeris", "=")

    # Global initial time
    global_initial_time = now()
    println("• Run started at ", global_initial_time)

    # Download technical comments file
    print("• Downloading technical comments file... ")
    Downloads.download(TECH_COMMENTS_URL, "de430_tech-comments.txt")
    tech_comments = read("de430_tech-comments.txt", String)
    rm("de430_tech-comments.txt")
    println("Done")

    # Download leap seconds kernel
    print("• Downloading leap seconds kernel... ")
    Downloads.download(LEAP_SECONDS_KERNEL_URL, "naif0012.tls")
    furnsh("naif0012.tls")
    rm("naif0012.tls")
    println("Done")

    # Download planets kernel
    print("• Downloading planets kernel... ")
    Downloads.download(PLANETS_KERNEL_URL, "de430_1850-2150.bsp")
    furnsh("de430_1850-2150.bsp")
    rm("de430_1850-2150.bsp")
    println("Done")

    # Download asteroids kernel
    print("• Downloading asteroids kernel... ")
    Downloads.download(ASTEROIDS_KERNEL_URL, "ast343de430.bsp")
    asteroids_ids = SpiceIntCell(1000)
    spkobj!(asteroids_ids, "ast343de430.bsp")
    furnsh("ast343de430.bsp")
    rm("ast343de430.bsp")
    println("Done")

    # Download TT-TDB kernel
    print("• Downloading TT-TDB kernel... ")
    Downloads.download(TTMTDB_KERNEL_URL, "TTmTDB.de430.19feb2015.bsp")
    furnsh("TTmTDB.de430.19feb2015.bsp")
    rm("TTmTDB.de430.19feb2015.bsp")
    println("Done")

    # Number of bodies (11 planets + 343 asteroids)
    N = 11 + 343
    # Pre-allocate memory
    ids = Vector{String}(undef, N)
    gms = Vector{Float64}(undef, N)
    rvs = Matrix{Float64}(undef, N, 6)
    # Initial epoch
    jd0 = parse_constant("JDEPOC", tech_comments)
    et = (jd0 - J2000) * daysec

    # Planets
    planets_ids = [10, 199, 299, 399, 301, 4, 5, 6, 7, 8, 9]
    @. ids[1:11] = string(planets_ids)
    for (i, id) in enumerate(planets_ids)
        name = planet_id2gm(id)
        gms[i] = parse_constant(name, tech_comments; c = 'D')
        rvs[i, :] .= kmsec2auday(spkgeo(id, et, "J2000", 0)[1])
    end

    # Asteroids
    @. ids[12:N] = string(asteroids_ids)
    for (i, id) in enumerate(asteroids_ids)
        name = "MA" * lpad(id - 2000000, 4, '0')
        gms[11+i] = parse_constant(name, tech_comments; c = 'E')
        rvs[11+i, :] .= kmsec2auday(spkgeo(id, et, "J2000", 0)[1])
    end
    # Sort asteroids by decreasing mass parameter
    perm = 11 .+ sortperm(view(gms, 12:N), rev = true)
    ids[12:N] .= ids[perm]
    rvs[12:N, :] .= rvs[perm, :]

    # Additional coordinates: lunar mantle and core Euler angles and angular velocities
    extras = Tuple{String, String, Float64, String}[]
    mantle_names = ["PHI", "THT", "PSI", "OMEGAX", "OMEGAY", "OMEGAZ"]
    for (name, coord, unit) in zip(mantle_names, LUNAR_COORDS, LUNAR_UNITS)
        push!(extras, ("LunarMantle", coord, parse_constant(name, tech_comments), unit))
    end
    core_names = ["PHIC", "THTC",  "PSIC", "OMGCX", "OMGCY", "OMGCZ"]
    for (name, coord, unit) in zip(core_names, LUNAR_COORDS, LUNAR_UNITS)
        push!(extras, ("LunarCore", coord, parse_constant(name, tech_comments), unit))
    end
    # and TT-TDB
    ttmtdb = spkgeo(1000000001, et, "J2000", 1000000000)[1][1]
    push!(extras, ("1000000001", "TT-TDB", ttmtdb, "s"))

    # Save initial conditions
    print("• Saving initial conditions to $output... ")
    sources = basename.([PLANETS_KERNEL_URL, ASTEROIDS_KERNEL_URL, TTMTDB_KERNEL_URL,
        TECH_COMMENTS_URL])
    write_initial_conditions(output, "DE430", jd0, sources, ids, rvs, extras)
    println("Done")

    # Final time
    global_final_time = now()
    println("• Run started ", global_initial_time, " and finished ", global_final_time)
    global_computation_time = computationtime(global_initial_time, global_final_time)
    println("• Total computation time was: ", global_computation_time, " min")

    return nothing
end

main()