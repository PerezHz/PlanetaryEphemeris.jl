# This file is part of the PlanetaryEphemeris.jl package; MIT licensed

"""
    read_initial_conditions(filename)

Read a vector of initial conditions from a file and return its epoch [JD TDB]
and the vector, i.e. `jd0, q0 = read_initial_conditions(filename)`.

# Extended help

The file must follow the `PlanetaryEphemeris.jl initial conditions v1` format,
written by `scripts/de430ic.jl` and `scripts/de440ic.jl`:
- A header of `# key = value` lines with, at least, keys `format`, `epoch` (the
  epoch of the initial conditions [JD TDB]) and `ncoords` (the length of the
  initial conditions vector). Other lines starting with `#` and blank lines are
  ignored.
- One `index id coordinate value unit` line per coordinate, where `index` is
  the position of `value` in the initial conditions vector.

For a planetary ephemeris model with `N` bodies, the vector of initial conditions
consists of `3N` positions [au], `3N` velocities [au/day] and the additional
coordinates of the model, e.g. 3 lunar mantle Euler angles [rad], 3 lunar mantle
angular velocities [rad/day], 3 lunar core Euler angles [rad], 3 lunar core
angular velocities [rad/day] and, for DE430, the initial value of TT-TDB [sec].
"""
function read_initial_conditions(filename::AbstractString)
    header = Dict{String, String}()
    idxs, vals = Int[], Float64[]
    for line in eachline(filename)
        line = strip(line)
        isempty(line) && continue
        # Header
        if startswith(line, '#')
            m = match(r"^#\s*(\w+)\s*=\s*([^#]*?)\s*(#.*)?$", line)
            isnothing(m) || (header[m[1]] = m[2])
            continue
        end
        # Coordinate
        fields = split(line)
        push!(idxs, parse(Int, fields[1]))
        push!(vals, parse(Float64, fields[4]))
    end
    # Check header
    fmt = get(header, "format", "")
    fmt == "PlanetaryEphemeris.jl initial conditions v1" ||
        throw(ArgumentError("Unknown initial conditions format: \"$fmt\""))
    haskey(header, "epoch") || throw(ArgumentError("Missing epoch in header"))
    jd0 = parse(Float64, header["epoch"])
    haskey(header, "ncoords") || throw(ArgumentError("Missing ncoords in header"))
    ncoords = parse(Int, header["ncoords"])
    # Check indices
    sort(idxs) == 1:ncoords || throw(ArgumentError("Indices of $filename are \
        not a permutation of 1:$ncoords"))
    # Fill the initial conditions vector
    q = Vector{Float64}(undef, ncoords)
    q[idxs] .= vals

    return jd0, q
end
