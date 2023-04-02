# This file is part of the PlanetaryEphemeris.jl package; MIT licensed

testfiles = (
    "propagation.jl",
    )

for file in testfiles
    include(file)
end
