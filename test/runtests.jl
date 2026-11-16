# This file is part of the PlanetaryEphemeris.jl package; MIT licensed

testfiles = (
    "orientation.jl",
    "propagation.jl",
    )

for file in testfiles
    include(file)
end
