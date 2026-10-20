# This file is part of the PlanetaryEphemeris.jl package; MIT licensed

"""
    read_parameters(filename)

Read the parameters of a planetary ephemeris model from a file and return them
as a `Dict{String, Union{Float64, Vector{Float64}}}`.

# Extended help

The file must follow the `PlanetaryEphemeris.jl parameters v1` format, written by
`scripts/de430ic.jl`:
- A header of `# key = value` lines with, at least, keys `format` and `nparams`
  (the number of parameter lines). Other lines starting with `#` and blank lines
  are ignored.
- One `name value unit source` line per parameter, where `name` is the name of
  the parameter in the technical comments file of the ephemeris. Array elements
  are written as `name[i]`, and the indices of each array must be `1:n`.
"""
function read_parameters(filename::AbstractString)
    header = Dict{String, String}()
    scalars = Dict{String, Float64}()
    arrays = Dict{String, Dict{Int, Float64}}()
    nparams = 0
    for line in eachline(filename)
        line = strip(line)
        isempty(line) && continue
        # Header
        if startswith(line, '#')
            m = match(r"^#\s*(\w+)\s*=\s*([^#]*?)\s*(#.*)?$", line)
            isnothing(m) || (header[m[1]] = m[2])
            continue
        end
        # Parameter
        nparams += 1
        fields = split(line)
        m = match(r"^(\w+)(?:\[(\d+)\])?$", fields[1])
        isnothing(m) && throw(ArgumentError("Invalid parameter name: $(fields[1])"))
        name, value = String(m[1]), parse(Float64, fields[2])
        if isnothing(m[2])
            haskey(scalars, name) && throw(ArgumentError("Repeated parameter: $name"))
            scalars[name] = value
        else
            i = parse(Int, m[2])
            array = get!(arrays, name, Dict{Int, Float64}())
            haskey(array, i) && throw(ArgumentError("Repeated parameter: $name[$i]"))
            array[i] = value
        end
    end
    # Check header
    fmt = get(header, "format", "")
    fmt == "PlanetaryEphemeris.jl parameters v1" ||
        throw(ArgumentError("Unknown parameters format: \"$fmt\""))
    haskey(header, "nparams") || throw(ArgumentError("Missing nparams in header"))
    nparams == parse(Int, header["nparams"]) || throw(ArgumentError("Number of \
        parameters in $filename ($nparams) does not match its header"))
    # Assemble arrays
    params = Dict{String, Union{Float64, Vector{Float64}}}(scalars)
    for (name, array) in arrays
        haskey(params, name) && throw(ArgumentError("Parameter $name is both a \
            scalar and an array"))
        n = length(array)
        sort!(collect(keys(array))) == 1:n || throw(ArgumentError("Indices of \
            $name are not 1:$n"))
        params[name] = [array[i] for i in 1:n]
    end

    return params
end
