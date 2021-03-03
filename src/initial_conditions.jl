function initialcond(N)
    # output from JPL Horizons

    q0 = Array{Float64}(undef, 6N+7) #initial condition array

    # fill initial conditions for Sun, Moon, planets and Pluto
    # order is the following:
    # Mercury, Venus, Earth-Moon barycenter, Moon (wrt geocenter), Mars, Jupiter, Saturn, Uranus, Neptune, Pluto
    smppic = readdlm( joinpath(dirname(pathof(PlanetaryEphemeris)), "ss11ic.txt") )
    for i in 1:11
        ith_body_ind = nbodyind(N, i)
        q0[ith_body_ind[1:3]] .= smppic[i, 1:3]
        q0[ith_body_ind[4:6]] .= smppic[i, 4:6]
    end
    # fill initial conditions for 343 asteroids used in integration of JPL DE430 ephemeris
    ast343ic = readdlm( joinpath(dirname(pathof(PlanetaryEphemeris)), "ast343ic.txt") )
    for i in 12:N
        ith_body_ind = nbodyind(N, i)
        q0[ith_body_ind[1:3]] .= ast343ic[i-11, 1:3]
        q0[ith_body_ind[4:6]] .= ast343ic[i-11, 4:6]
    end
    # # 1969-Jun-28.0 (TDB)
    # q0[6N+1:6N+3] .= [0.00512830031411853500, 0.38239278420173690000, 1.29416700274878300000] # Euler angles
    # q0[6N+4:6N+6] .= [0.00004573724185991433, -0.00000218986174567295, 0.22994486018992250000] # lunar mantle angular velocity vector (mantle frame)
    # q0[6N+7] = -0.00016266592104301078 # DE430 TT-TDB at initial epoch
    # J2000
    q0[6N+1:6N+3] .= [-0.05414766382529318, 0.42485573826608863, 2564.2582726674223] # Euler angles
    q0[6N+4:6N+6] .= [2.3013404932266894e-6, -6.600217715260397e-5, 0.22999341817950175] # lunar mantle angular velocity vector (mantle frame)
    q0[6N+7] = 9.930292723454279e-5 # DE430 TT-TDB at initial epoch
    return q0
end
