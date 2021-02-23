function initialcond(N)
    # output from JPL Horizons

    q0 = Array{Float64}(undef, 6N) #initial condition array

    # fill initial conditions for Sun, Moon, planets and Pluto
    # order is the following:
    # Mercury, Venus, Earth-Moon barycenter, Moon (wrt geocenter), Mars, Jupiter, Saturn, Uranus, Neptune, Pluto
    smppic = readdlm( joinpath(dirname(pathof(PlanetaryEphemeris)), "ss11ic.txt") )
    for i in 1:11
        q0[3i-2:3i]         .= smppic[i, 1:3]
        q0[3(N+i)-2:3(N+i)] .= smppic[i, 4:6]
    end
    # fill initial conditions for 343 asteroids used in integration of JPL DE430 ephemeris
    ast343ic = readdlm( joinpath(dirname(pathof(PlanetaryEphemeris)), "ast343ic.txt") )
    for i in 12:N
        q0[3i-2:3i]         .= ast343ic[i-11, 1:3]
        q0[3(N+i)-2:3(N+i)] .= ast343ic[i-11, 4:6]
    end
    return q0
end
