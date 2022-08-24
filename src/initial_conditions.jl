@doc raw"""
    initialcond(N, jd0::Float64=datetime2julian(DateTime(1969,6,28,0,0,0)))

Returns a vector with the initial conditions (`3N` positions [au] + `3N` velocities [au/day] + 
3 lunar mantle Euler angles [rad] + 3 mantle angular velocities [rad/day] +
3 lunar core Euler angles [rad] + 3 core angular velocities [rad/day] +
DE430 TT-TDB at initial epoch [days]) for the integration. Two possible values of `jd0` 
are supported:

- If `jd0 == datetime2julian(DateTime(1969,6,28,0,0,0))`, planets (+ Sun & Moon) and asteroids initial conditions are retrieved from files `ss11ic_1969Jun28.txt` and `ast343ic_1969Jun28.txt` respectively.

- If `jd0 == datetime2julian(DateTime(2000,1,1,12))`, planets (+ Sun & Moon) and asteroids initial conditions are retrieved from files `ss11ic.txt` and `ast343ic.txt` respectively.

For the initial conditions at 1969-Jun-28.0 see Tables 5 and 6 in page 47 (Sun + Planets + 
Moon positions and velocities), Table 7 in page 49 (lunar mantle and core libration 
angles/rates) and Table 13 in pages 60-74 (asteroids positions and velocities) of 
https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract.
"""
function initialcond(N, jd0::Float64=datetime2julian(DateTime(1969,6,28,0,0,0)))
    # Output from JPL Horizons

    q0 = Array{Float64}(undef, 6N+13)     # Initial conditions array
    dt0 = julian2datetime(jd0)            # Convert jd0 to DateTime
    dt0_1969 = DateTime(1969,6,28,0,0,0)  # 1969-Jun-28 (TDB)
    dt0_2000 = DateTime(2000,1,1,12)      # 2000-Jan-1.5 (TDB)

    # 1969-Jun-28 (TDB)
    if dt0 == dt0_1969
        println("Initial conditions: 1969-Jun-28.0 (TDB)")
        ss_ic_filename = "ss11ic_1969Jun28.txt"     # Planets (+ Sun & Moon) initial conditions file
        ast_ic_filename = "ast343ic_1969Jun28.txt"  # Asteroids initial conditions file
        # Initial conditions for lunar mantle and core libration angles/rates
        # See Table 7 in page 49 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract
        q0[6N+1:6N+3] .= [0.00512830031411853500, 0.38239278420173690000, 1.29416700274878300000]      # Lunar mantle Euler angles
        q0[6N+4:6N+6] .= [0.00004573724185991433, -0.00000218986174567295, 0.22994486018992250000]     # Lunar mantle angular velocity vector (mantle frame)
        q0[6N+7:6N+9] .= [-0.00241990927040684100, 0.41101946488652730000, -0.46309468558363680000]    # Lunar core Euler angles
        q0[6N+10:6N+12] .= [-0.00661836772247824400, -0.00107295445159005100, 0.22964879652299730000]  # Lunar core angular velocity vector (mantle frame)
        # DE430 TT-TDB at initial epoch
        q0[6N+13] = -0.00016266592104301078 
    # 2000-Jan-1.5 (TDB)
    elseif dt0 == dt0_2000
        println("Initial conditions: 2000-Jan-1.5 (TDB)")
        ss_ic_filename = "ss11ic.txt"               # Planets (+ Sun & Moon) initial conditions file
        ast_ic_filename = "ast343ic.txt"            # Asteroids initial conditions file
        # Initial conditions for lunar mantle and core libration angles/rates
        q0[6N+1:6N+3] .= [-0.05414766382529318, 0.42485573826608863, 2564.2582726674223]       # Lunar mantle Euler angles
        q0[6N+4:6N+6] .= [2.3013404932266894e-6, -6.600217715260397e-5, 0.22999341817950175]   # Lunar mantle angular velocity vector (mantle frame)
        # TODO: update initial values of core angles/omegas
        q0[6N+7:6N+9] .= [0.0025514083634858424, 0.4093771628786826, 2560.1154745041295]       # Lunar core Euler angles
        q0[6N+10:6N+12] .= [0.006443733086714832, -0.0005981799500672227, 0.22968509873693774] # Lunar core angular velocity vector (mantle frame)
        # DE430 TT-TDB at initial epoch
        q0[6N+13] = 9.930292723454279e-5 
    # Neither 1969-Jun-28 (TDB) or 2000-Jan-1.5 (TDB)
    else
        @assert(false, "Initial time must either correspond to $(dt0_1969) or $(dt0_2000).")
    end

    # Fill initial conditions for Sun, Moon, planets and Pluto
    # Order is the following:
    # Sun, Mercury, Venus, Earth-Moon barycenter, Moon (wrt geocenter), Mars, Jupiter, Saturn, Uranus, Neptune, Pluto
    smppic = readdlm( joinpath(dirname(pathof(PlanetaryEphemeris)), ss_ic_filename) )
    for i in 1:11
        ith_body_ind = nbodyind(N, i)             # Indices of the i-th body
        q0[ith_body_ind[1:3]] .= smppic[i, 1:3]   # Position vector
        q0[ith_body_ind[4:6]] .= smppic[i, 4:6]   # Velocity vector
    end
    # Fill initial conditions for 343 asteroids used in integration of JPL DE430 ephemeris
    ast343ic = readdlm( joinpath(dirname(pathof(PlanetaryEphemeris)), ast_ic_filename) )
    for i in 12:N
        ith_body_ind = nbodyind(N, i)                  # Indices of the i-th body
        q0[ith_body_ind[1:3]] .= ast343ic[i-11, 1:3]   # Position vector
        q0[ith_body_ind[4:6]] .= ast343ic[i-11, 4:6]   # Velocity vector
    end

    return q0
end
