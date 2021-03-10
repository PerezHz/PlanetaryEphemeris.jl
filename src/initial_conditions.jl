function initialcond(N)
    # output from JPL Horizons

    q0 = Array{Float64}(undef, 6N+13) #initial condition array

    # fill initial conditions for Sun, Moon, planets and Pluto
    # order is the following:
    # Mercury, Venus, Earth-Moon barycenter, Moon (wrt geocenter), Mars, Jupiter, Saturn, Uranus, Neptune, Pluto
    smppic = readdlm( joinpath(dirname(pathof(PlanetaryEphemeris)), "ss11ic_1969Jun28.txt") ) # "ss11ic.txt"
    for i in 1:11
        ith_body_ind = nbodyind(N, i)
        q0[ith_body_ind[1:3]] .= smppic[i, 1:3]
        q0[ith_body_ind[4:6]] .= smppic[i, 4:6]
    end
    # fill initial conditions for 343 asteroids used in integration of JPL DE430 ephemeris
    ast343ic = readdlm( joinpath(dirname(pathof(PlanetaryEphemeris)), "ast343ic_1969Jun28.txt") ) # "ast343ic.txt"
    for i in 12:N
        ith_body_ind = nbodyind(N, i)
        q0[ith_body_ind[1:3]] .= ast343ic[i-11, 1:3]
        q0[ith_body_ind[4:6]] .= ast343ic[i-11, 4:6]
    end
    # 1969-Jun-28.0 (TDB)
    q0[6N+1:6N+3] .= [0.00512830031411853500, 0.38239278420173690000, 1.29416700274878300000] # lunar mantle Euler angles
    q0[6N+4:6N+6] .= [0.00004573724185991433, -0.00000218986174567295, 0.22994486018992250000] # lunar mantle angular velocity vector (mantle frame)
    q0[6N+7:6N+9] .= [-0.00241990927040684100, 0.41101946488652730000, -0.46309468558363680000] # lunar core Euler angles
    q0[6N+10:6N+12] .= [-0.00661836772247824400, -0.00107295445159005100, 0.22964879652299730000] # lunar core angular velocity vector (mantle frame)
    q0[6N+13] = -0.00016266592104301078 # DE430 TT-TDB at initial epoch

    # # J2000
    # q0[6N+1:6N+3] .= [-0.05414766382529318, 0.42485573826608863, 2564.2582726674223] # Euler angles
    # q0[6N+4:6N+6] .= [2.3013404932266894e-6, -6.600217715260397e-5, 0.22999341817950175] # lunar mantle angular velocity vector (mantle frame)
    # # TODO: update initial values of core angles/omegas
    # q0[6N+7:6N+9] .= [0.0025514083634858424, 0.4093771628786826, 2560.1154745041295] # lunar core Euler angles
    # q0[6N+10:6N+12] .= [0.006443733086714832, -0.0005981799500672227, 0.22968509873693774] # lunar core angular velocity vector (mantle frame)
    # q0[6N+13] = 9.930292723454279e-5 # DE430 TT-TDB at initial epoch
    return q0
end
