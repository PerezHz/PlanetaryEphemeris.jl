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

    # Number of bodies
    N = 11 + 343 + 2 + 1 # planets + asteroids + lunar core/mantle + TT-TDB
    # Pre-allocate memory
    ids = Vector{String}(undef, N)
    qs = Matrix{Float64}(undef, N, 7)
    # Initial epoch
    jd0 = parse_constant("JDEPOC", tech_comments)
    et = (jd0 - J2000) * daysec

    # Planets
    planets_ids = [10, 199, 299, 399, 301, 4, 5, 6, 7, 8, 9]
    @. ids[1:11] = string(planets_ids)
    for (i, id) in enumerate(planets_ids)
        name = planet_id2gm(id)
        qs[i, 1] = parse_constant(name, tech_comments; c = 'D')
        qs[i, 2:end] .= kmsec2auday(spkgeo(id, et, "J2000", 0)[1])
    end

    # Asteroids
    @. ids[12:354] = string(asteroids_ids)
    for (i, id) in enumerate(asteroids_ids)
        name = "MA" * lpad(id - 2000000, 4, '0')
        qs[11+i, 1] = parse_constant(name, tech_comments; c = 'E')
        qs[11+i, 2:end] .= kmsec2auday(spkgeo(id, et, "J2000", 0)[1])
    end
    perm = sortperm(view(qs, 11+1:11+343, 1), rev = true)
    permute!(view(ids, 11+1:11+343), perm)
    for i in axes(qs, 2)
        permute!(view(qs, 11+1:11+343, i), perm)
    end

    # Lunar mantle euler angles
    mantle_names = ["PHI", "THT", "PSI", "OMEGAX", "OMEGAY", "OMEGAZ"]
    ids[end-2, 1] = "LunarMantle"
    qs[end-2, 1] = NaN
    @. qs[end-2, 2:end] = parse_constant(mantle_names, tech_comments)

    # Lunar core euler angles
    core_names = ["PHIC", "THTC",  "PSIC", "OMGCX", "OMGCY", "OMGCZ"]
    ids[end-1, 1] = "LunarCore"
    qs[end-1, 1] = NaN
    @. qs[end-1, 2:end] = parse_constant(core_names, tech_comments)

    # TT-TDB
    ids[end] = "1000000001"
    qs[end, :] .= NaN, spkgeo(1000000001, et, "J2000", 1000000000)[1][1], NaN, NaN,
        NaN, NaN, NaN

    # Save initial conditions vector
    print("• Saving initial conditions to $output... ")
    colnames = ["Spice ID", "GM [au^3/day^2]", "x [au]", "y [au]", "z [au]",
        "vx [au/day]", "vy [au/day]", "vz [au/day]"]
    open(output, "w") do file
        write(file, join(rpad(c, 28) for c in colnames), '\n')
        lines = Vector{String}(undef, N)
        for i in eachindex(ids)
            x = [rpad(@sprintf("%+.16E", y), 28) for y in qs[i, :]]
            n = i == N ? "" : "\n"
            lines[i] = string(rpad(ids[i], 28), join(x), n)
        end
        write(file, join(lines))
    end
    println("Done")

    # Final time
    global_final_time = now()
    println("• Run started ", global_initial_time, " and finished ", global_final_time)
    global_computation_time = computationtime(global_initial_time, global_final_time)
    println("• Total computation time was: ", global_computation_time, " min")

    return nothing
end

main()