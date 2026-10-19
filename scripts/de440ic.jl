using ArgParse, PlanetaryEphemeris, Dates
using Downloads, SPICE, Printf

const TECH_COMMENTS_URL = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440_tech-comments.txt"
const LEAP_SECONDS_KERNEL_URL = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/naif0012.tls"
const PLANETS_KERNEL_URL = "https://ssd.jpl.nasa.gov/ftp/eph/planets/bsp/de440t.bsp"
const ASTEROIDS_KERNEL_URL = "https://ssd.jpl.nasa.gov/ftp/eph/small_bodies/asteroids_de441/sb441-n373s.bsp"

# Mean obliquity of the ecliptic in degrees
# See equation (5-153) in page 5-61 of https://doi.org/10.1002/0471728470
const ϵ = deg2rad(84_381.448 / 3_600)
# Rotation matrix from the ecliptic plane to the equatorial plane
const M_ecl2eq = hcat(vcat(Rx(-ϵ), zeros(3, 3)), vcat(zeros(3, 3), Rx(-ϵ)))

# See https://ssd.jpl.nasa.gov/ftp/eph/planets/bsp/README.txt
const KBOs_SPKID_TO_ID = Dict{Int, Int}(
    2136199 => 8001,
    2136108 => 8002,
    2090377 => 8003,
    2136472 => 8004,
    2050000 => 8005,
    2084522 => 8006,
    2090482 => 8007,
    2020000 => 8008,
    2055637 => 8009,
    2028978 => 8010,
    2307261 => 8011,
    2174567 => 8012,
    3361580 => 8013,
    3308265 => 8014,
    2055565 => 8015,
    2145452 => 8016,
    2090568 => 8017,
    2208996 => 8018,
    2225088 => 8019,
    2019521 => 8020,
    2120347 => 8021,
    2278361 => 8022,
    3525142 => 8023,
    2230965 => 8024,
    2042301 => 8025,
    2455502 => 8026,
    3545742 => 8027,
    2523639 => 8028,
    2528381 => 8029,
    3515022 => 8030
)

function parse_commandline()

    s = ArgParseSettings()

    # Program name (for usage & help screen)
    s.prog = "de440ic.jl"
    # Desciption (for help screen)
    s.description = "Generate the initial conditions file for the DE440 ephemeris"

    @add_arg_table! s begin
        "--output", "-o"
            help = "Output file"
            arg_type = String
            default = "de440ic.txt"
    end

    s.epilog = """
        examples:\n
        \n
        julia --project de440ic.jl -o de440ic.txt\n
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

function spkid2ma(id::Int32)
    # Kuiper belt
    if id < 0
        return string("MA82", lpad(abs(id), 2, '0'))
    # Asteroid
    elseif 0 < id < 2010000
        return string("MA", lpad(id - 2000000, 4, '0'))
    # KBO
    elseif 2010000 < id && haskey(KBOs_SPKID_TO_ID, id)
        return string("MA", KBOs_SPKID_TO_ID[id])
    else
        throw(ArgumentError("Cannot convert SPKID to MA"))
    end
end

function main()

    # Parse arguments from commandline
    parsed_args = parse_commandline()
    # Output file
    output::String = parsed_args["output"]

    printitle("Generate the initial conditions file for the DE440 ephemeris", "=")

    # Global initial time
    global_initial_time = now()
    println("• Run started at ", global_initial_time)

    # Download technical comments file
    print("• Downloading technical comments file... ")
    Downloads.download(TECH_COMMENTS_URL, "de440_tech-comments.txt")
    tech_comments = read("de440_tech-comments.txt", String)
    rm("de440_tech-comments.txt")
    println("Done")

    # Download leap seconds kernel
    print("• Downloading leap seconds kernel... ")
    Downloads.download(LEAP_SECONDS_KERNEL_URL, "naif0012.tls")
    furnsh("naif0012.tls")
    rm("naif0012.tls")
    println("Done")

    # Download planets kernel
    print("• Downloading planets kernel... ")
    Downloads.download(PLANETS_KERNEL_URL, "de440t.bsp")
    furnsh("de440t.bsp")
    rm("de440t.bsp")
    println("Done")

    # Download asteroids kernel
    print("• Downloading asteroids kernel... ")
    Downloads.download(ASTEROIDS_KERNEL_URL, "sb441-n373s.bsp")
    asteroids_ids = SpiceIntCell(1000)
    spkobj!(asteroids_ids, "sb441-n373s.bsp")
    append!(asteroids_ids, -1:-1:-36)
    furnsh("sb441-n373s.bsp")
    rm("sb441-n373s.bsp")
    println("Done")

    # Number of bodies (11 planets + 343 asteroids + 30 KBOs + 36 masses
    # representing the Kuiper belt)
    N = 11 + 343 + 30 + 36
    # Pre-allocate memory
    ids = Vector{String}(undef, N)
    gms = Vector{Float64}(undef, N)
    rvs = Matrix{Float64}(undef, N, 6)
    # Initial epoch
    jd0 = parse_constant("JDEPOC", tech_comments; c = 'D')
    et = (jd0 - J2000) * daysec

    # Planets
    planets_ids = [10, 199, 299, 399, 301, 4, 5, 6, 7, 8, 9]
    @. ids[1:11] = string(planets_ids)
    for (i, id) in enumerate(planets_ids)
        name = planet_id2gm(id)
        gms[i] = parse_constant(name, tech_comments; c = 'e')
        rvs[i, :] .= kmsec2auday(spkgeo(id, et, "J2000", 0)[1])
    end

    # Asteroids, KBOs and Kuiper belt
    for (i, id) in enumerate(asteroids_ids)
        ids[11+i] = id < 0 ? string("Kuiper", abs(id)) : string(id)
        name = spkid2ma(id)
        gms[11+i] = parse_constant(name, tech_comments; c = 'e')
        # Kuiper belt
        if id < 0
            a, μ = 44.0, sum(view(gms, 11+1:11+373))
            ν = deg2rad(10 * (abs(id)-1))
            q = [a * cos(ν), a * sin(ν), 0.0, -sqrt(μ/a) * sin(ν), sqrt(μ/a) * cos(ν), 0.0]
            rvs[11+i, :] .= M_ecl2eq * q
        # Asteroid or KBO
        else
            rvs[11+i, :] .= kmsec2auday(spkgeo(id, et, "J2000", 0)[1])
        end
    end
    # Sort asteroids, KBOs and Kuiper belt masses by decreasing mass parameter
    perm = 11 .+ sortperm(view(gms, 12:N), rev = true)
    ids[12:N] .= ids[perm]
    rvs[12:N, :] .= rvs[perm, :]

    # Additional coordinates: lunar mantle and core Euler angles and angular velocities
    extras = Tuple{String, String, Float64, String}[]
    mantle_names = ["PHI", "THT", "PSI", "OMEGAX", "OMEGAY", "OMEGAZ"]
    for (name, coord, unit) in zip(mantle_names, LUNAR_COORDS, LUNAR_UNITS)
        value = parse_constant(name, tech_comments; c = 'D')
        push!(extras, ("LunarMantle", coord, value, unit))
    end
    core_names = ["PHIC", "THTC",  "PSIC", "OMGCX", "OMGCY", "OMGCZ"]
    for (name, coord, unit) in zip(core_names, LUNAR_COORDS, LUNAR_UNITS)
        value = parse_constant(name, tech_comments; c = 'D')
        push!(extras, ("LunarCore", coord, value, unit))
    end

    # Save initial conditions
    print("• Saving initial conditions to $output... ")
    sources = basename.([PLANETS_KERNEL_URL, ASTEROIDS_KERNEL_URL, TECH_COMMENTS_URL])
    write_initial_conditions(output, "DE440", jd0, sources, ids, rvs, extras)
    println("Done")

    # Final time
    global_final_time = now()
    println("• Run started ", global_initial_time, " and finished ", global_final_time)
    global_computation_time = computationtime(global_initial_time, global_final_time)
    println("• Total computation time was: ", global_computation_time, " min")

    return nothing
end

main()