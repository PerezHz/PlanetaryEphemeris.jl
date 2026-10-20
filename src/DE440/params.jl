"""
    DE440Constants{T <: Real}

Physical parameters of the JPL DE440 dynamical model. Parameters listed in the
technical comments file of DE440 keep their names there.

    DE440Constants{T}(filename)

Read the parameters in `filename` (see [`read_parameters`](@ref)), convert them to `T`
and derive the remaining quantities used by the model.

# Fields

- `N_ext::Int`: number of bodies in extended-body accelerations.
- `GM::Vector{T}`: mass parameters [au^3/day^2].
- `UJ_interaction::Matrix{Bool}`: whether the extended body `j` (column) interacts
    with the point mass `i` (row).
- `SUNPOLERA::T`, `SUNPOLEDEC::T`: right ascension and declination of the Sun's
    rotation pole [deg].
- `R_sun::T`: Sun's radius [au].
- `MOISUN::T`: Sun's polar moment of inertia divided by `M R^2`.
- `OMGSUN::T`: Sun's rotation rate [rad/day].
- `EPSSRPE::T`, `EPSSRPM::T`: solar radiation pressure factors of the Earth and
    the Moon.
- `RE_au::T`: Earth's radius [au].
- `J2E::T`, `J2EDOT::T`: Earth's second zonal harmonic and its rate of change [1/yr].
- `ROTEX::T`, `ROTEY::T`: rotations about the X and Y axes in the Earth orientation
    model [arcsec].
- `JSEM::Matrix{T}`: zonal harmonics `J_n` times the radius of the body [au] to the
    `n`-th power (row: body, column: degree).
- `n1SEM::Vector{Int}`: degree of the zonal harmonics expansion of each extended body.
- `n2M::Int`: degree of the Moon's tesseral harmonics expansion.
- `CM::Matrix{T}`, `SM::Matrix{T}`: Moon's tesseral harmonics `C_{nm}`, `S_{nm}` times
    the Moon's radius [au] to the `n`-th power.
- `R_moon::T`: Moon's radius [au].
- `K2M::T`: Moon's potential Love number.
- `TAUM::T`: time-lag of the lunar solid-body tide [day].
- `n_moon::T`: Moon's mean motion [rad/day].
- `ITM_und::Matrix{T}`: Moon's undistorted total moment of inertia.
- `I_c::Matrix{T}`: Moon's core moment of inertia.
- `k_ν::T`: friction coefficient between the Moon's core and mantle.
- `C_c_m_A_c::T`: numerical factor of the torque on the Moon's mantle due to
    the interaction between the core and mantle.
- `K2E0::T`, `K2E1::T`, `K2E2::T`: Earth's potential Love numbers for
    long-period, diurnal and semi-diurnal deformations.
- `TAUE0::T`, `TAUE1::T`, `TAUE2::T`: Earth's orbital time-lags for long-period,
    diurnal and semi-diurnal deformations [day].
- `TAUER0::T`, `TAUER1::T`, `TAUER2::T`: Earth's rotational time-lags for long-period,
    diurnal and semi-diurnal deformations [day].
- `tid_num_coeff::T`: overall numerical factor of the tidal acceleration of the
    Moon, `1.5 * (1 + GM[mo] / GM[ea])`.
"""
struct DE440Constants{T <: Real}
    N_ext::Int
    GM::Vector{T}
    UJ_interaction::Matrix{Bool}
    SUNPOLERA::T
    SUNPOLEDEC::T
    R_sun::T
    MOISUN::T
    OMGSUN::T
    EPSSRPE::T
    EPSSRPM::T
    RE_au::T
    J2E::T
    J2EDOT::T
    ROTEX::T
    ROTEY::T
    JSEM::Matrix{T}
    n1SEM::Vector{Int}
    n2M::Int
    CM::Matrix{T}
    SM::Matrix{T}
    R_moon::T
    K2M::T
    TAUM::T
    n_moon::T
    ITM_und::Matrix{T}
    I_c::Matrix{T}
    k_ν::T
    C_c_m_A_c::T
    K2E0::T
    K2E1::T
    K2E2::T
    TAUE0::T
    TAUE1::T
    TAUE2::T
    TAUER0::T
    TAUER1::T
    TAUER2::T
    tid_num_coeff::T

    function DE440Constants{T}(filename::AbstractString) where {T <: Real}
        # Parameters in file
        params = read_parameters(filename)
        # Return parameter `name` converted to `T`, keeping track of the used ones
        used = Set{String}()
        function param(name::String)
            haskey(params, name) || throw(ArgumentError("Missing parameter \
                $name in $filename"))
            push!(used, name)
            return T.(params[name])
        end

        # Mass parameters
        GM = param("GM")
        N = length(GM)

        # Model structure
        # See Section 3.2 in pages 5-6 of https://doi.org/10.3847/1538-3881/abd414
        # Number of bodies in extended-body accelerations
        N_ext = 11
        # Degree of zonal harmonics expansion for each body with extended body accelerations
        n1SEM = [2, 0, 0, 5, 6]   # Sun, Mercury, Venus, Earth, Moon
        # Degree of tesseral harmonics expansion for the Moon
        n2M = 6
        # Matrix of extended body interactions
        UJ_interaction = fill(false, N, N)
        # Sun's J_2 interacts with all other bodies
        UJ_interaction[2:N, su] .= true
        # Earth's grav potential interacts with Sun, Mercury, Venus, Moon, Mars, Jupiter and Saturn
        UJ_interaction[union(1:ea-1,ea+1:8), ea] .= true
        # Moon's grav potential interacts with Sun, Mercury, Venus, Earth, Mars, Jupiter and Saturn
        UJ_interaction[union(1:mo-1,mo+1:8), mo] .= true

        # Sun
        SUNPOLERA, SUNPOLEDEC = param("SUNPOLERA"), param("SUNPOLEDEC")
        # Solar radius in au
        R_sun = param("ASUN")/au
        # Second zonal harmonic coefficient J_2 * Suns's radius in au^2
        JS = [param("J2SUN")*R_sun^2]
        # Lense-Thirring and solar radiation pressure
        MOISUN, OMGSUN = param("MOISUN"), param("OMGSUN")
        EPSSRPE, EPSSRPM = param("EPSSRPE"), param("EPSSRPM")

        # Earth
        # Earth's radius in au
        RE_au = param("RE")/au
        J2E, J2EDOT = param("J2E"), param("J2EDOT")
        # Vector of zonal harmonic coefficients J_n * Earth's radius in au ^n
        JE = [J2E*RE_au^2, param("J3E")*RE_au^3, param("J4E")*RE_au^4, param("J5E")*RE_au^5]
        # Orientation offsets
        ROTEX, ROTEY = param("ROTEX"), param("ROTEY")
        # Love numbers, orbital and rotational time-lags
        K2E0, K2E1, K2E2 = param("K2E0"), param("K2E1"), param("K2E2")
        TAUE0, TAUE1, TAUE2 = param("TAUE0"), param("TAUE1"), param("TAUE2")
        TAUER0, TAUER1, TAUER2 = param("TAUER0"), param("TAUER1"), param("TAUER2")

        # Moon
        # Lunar radius in au and its powers
        R_moon = param("AM")/au
        R_moon_n = (R_moon, R_moon^2, R_moon^3, R_moon^4, R_moon^5, R_moon^6)
        K2M, TAUM = param("K2M"), param("TAUM")
        # Vector of zonal harmonic coefficients J_n * Moons's radius in au ^n
        JM = [param("J$(n)M")*R_moon_n[n] for n in 2:6]
        # Matrices of tesseral harmonic coefficients C_{nm}, S_{nm} * Moon's radius in au ^n
        CM = zeros(T, 6, 6)
        SM = zeros(T, 6, 6)
        for n in 2:n2M
            CM[n, 1:n] .= [param("C$(n)$(m)M") for m in 1:n]*R_moon_n[n]
            SM[n, 1:n] .= [param("S$(n)$(m)M") for m in 1:n]*R_moon_n[n]
        end
        # Matrix of zonal harmonic coefficients J_n * radius of the corresponding body ^n
        JSEM = zeros(T, 5, 6)
        JSEM[su, 2:2] = JS
        JSEM[ea, 2:5] = JE
        JSEM[mo, 2:6] = JM
        # Diagonal elements of undistorted lunar mantle moment of inertia
        # See equation (49) in page 8 of https://doi.org/10.3847/1538-3881/abd414
        J2M, LBET, LGAM = param("J2M"), param("LBET"), param("LGAM")
        A_T = (2(1-LBET*LGAM)/(2LBET-LGAM+LBET*LGAM))*J2M
        B_T = (2(1+LGAM)/(2LBET-LGAM+LBET*LGAM))*J2M
        C_T = (2(1+LBET)/(2LBET-LGAM+LBET*LGAM))*J2M
        # Friction coefficient between core and mantle
        k_ν = param("KVC")*(C_T*GM[mo]*(R_moon^2))
        # Diagonal elements of lunar core moment of inertia
        IFAC, COBLAT = param("IFAC"), param("COBLAT")
        A_c = IFAC*C_T*(1-COBLAT)
        B_c = IFAC*C_T*(1-COBLAT)
        C_c = IFAC*C_T
        # Numerical factor for the torque on the mantle due to the interaction between the core and mantle
        C_c_m_A_c = (C_c - A_c)*(GM[mo]*R_moon^2)
        # Lunar undistorted total moment of inertia
        ITM_und = diagm([A_T, B_T, C_T]*GM[mo]*R_moon^2)
        # Lunar core moment of inertia
        I_c = diagm([A_c, B_c, C_c]*GM[mo]*R_moon^2)
        # Lunar mean motion (rad/day)
        n_moon = sqrt((GM[ea]+GM[mo])/((param("AMOON")/au)^3))
        # Overall numerical factor of the tidal acceleration of the Moon
        tid_num_coeff = 1.5*(1.0 + GM[mo]/GM[ea])

        # Check that all the parameters in file were used
        unused = setdiff(keys(params), used)
        isempty(unused) || throw(ArgumentError("Unknown parameters in $filename: \
            $(join(sort!(collect(unused)), ", "))"))

        return new{T}(N_ext, GM, UJ_interaction, SUNPOLERA, SUNPOLEDEC, R_sun, MOISUN,
            OMGSUN, EPSSRPE, EPSSRPM, RE_au, J2E, J2EDOT, ROTEX, ROTEY, JSEM, n1SEM, n2M,
            CM, SM, R_moon, K2M, TAUM, n_moon, ITM_und, I_c, k_ν, C_c_m_A_c, K2E0, K2E1,
            K2E2, TAUE0, TAUE1, TAUE2, TAUER0, TAUER1, TAUER2, tid_num_coeff)
    end
end

# Print method for DE440Constants
show(io::IO, x::DE440Constants) = print(io, "DE440Constants")

"""
    DE440Params{T <: Real}

Parameters for the JPL DE440 planetary ephemeris model.

# Fields

- `N::Int`: total number of bodies.
- `N_bwd::Int`: number of bodies used to compute time-delayed
    tidal interactions.
- `jd0::T`: reference epoch [JDTDB].
- `consts::DE440Constants{T}`: physical parameters of the model.
"""
struct DE440Params{T <: Real}
    N::Int
    N_bwd::Int
    jd0::T
    consts::DE440Constants{T}
end

# Outer constructor
function DE440Params(jd0::T; consts::DE440Constants{T} = DE440Constants{T}(joinpath(
                     pkgdir(PlanetaryEphemeris), "data", "de440params.txt"))
                     ) where {T <: Real}
    # Total number of bodies
    N = length(consts.GM)
    # Number of bodies used to compute time-delayed tidal interactions
    N_bwd = 11

    return DE440Params{T}(N, N_bwd, jd0, consts)
end

# Print method for DE440Params
show(io::IO, x::DE440Params) = print(io, "DE440Params")
