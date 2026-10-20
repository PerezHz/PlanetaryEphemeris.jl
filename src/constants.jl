# PlanetaryEphemeris abbreviation
const PE = PlanetaryEphemeris

# Path to PlanetaryEphemeris src directory
const src_path = dirname(pathof(PlanetaryEphemeris))

# Integration parameters
const _order = 25
const abstol = 1.0E-20

# Important bodies indices
const su = 1 # Sun's index
const ea = 4 # Earth's index
const mo = 5 # Moon's index

# Physical constants in various units
# See Table 4 in page 47 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract
const au = 1.495978707E8                       # Astronomical unit value in km
const yr = 365.25                              # Days in a Julian year
const daysec = 86_400                          # Number of seconds in a day

# Parameters related to speed of light, c
const clightkms = 2.99792458E5                 # Speed of light, km/sec
const c_au_per_day = daysec*(clightkms/au)     # Speed of light in au per day
const c_au_per_sec = clightkms/au              # Speed of light in au per sec
const c_cm_per_sec = 100_000*clightkms         # Speed of light in cm per sec
const c_p2 = 29979.063823897606                # Speed of light^2 in au^2/day^2
const c_m2 = 3.3356611996764786e-5             # Speed of light^-2 in day^2/au^2
const c_m4 = 1.1126635639027125e-9             # Speed of light^-4 in day^4/au^4

"""
    nbodyind(N, idxs)

Return the indices of the positions and velocities of the `idxs`-th bodies
in a vector of the form:  `3N` positions + `3N` velocities + Lunar physical
librations + TT-TDB.
"""
nbodyind(N::Int, i::Int) = union(3i-2:3i, 3*(N+i)-2:3*(N+i))

function nbodyind(N::Int, idxs::AbstractVector{Int})
    jdxs = mapreduce(Base.Fix1(nbodyind, N), vcat, idxs)
    sort!(jdxs)
    return sort(jdxs)
end

const J2000 = 2451545.0                        # 2000 epoch julian date

# Numerical factors for recursion relations of Legendre polynomials up to degree 6
# See equations (175)-(178) in page 33 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract
const fact1_jsem = [(2n-1)/n for n in 1:6]  # (2n - 1) / n
const fact2_jsem = [(n-1)/n for n in 1:6]   # (n - 1) / n
const fact3_jsem = [n for n in 1:6]         # n
const fact4_jsem = fact3_jsem .+ 1          # n + 1
const fact5_jsem = fact3_jsem .+ 2          # n + 2

# Numerical factors for recursion relations of Associated Legendre polynomials
# See equations (180)-(183) in pages 33-34 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract
const lnm1 = [(2n-1)/(n-m) for n in 1:6, m in 1:6]       # (2n - 1) / (n - m)
const lnm2 = [-(n+m-1)/(n-m) for n in 1:6, m in 1:6]     # -(n + m - 1) / (n - m)
const lnm3 = [-n for n in 1:6]                           # -n
const lnm4 = [n+m for n in 1:6, m in 1:6]                # (n + m)
const lnm5 = [2n-1 for n in 1:6]                         # (2n - 1)
const lnm6 = [-(n+1) for n in 1:6]                       # -(n + 1)
const lnm7 = [m for m in 1:6]                            # m

# Standard value of nominal mean angular velocity of Earth (rad/day),
# See Explanatory Supplement to the Astronomical Almanac 2014, section 7.4.3.3,
# page 296: 7.2921151467e-5 rad/second
const ω_E = daysec*7.2921151467e-5 # (7.2921151467e-5 rad/sec)*daysec -> rad/day

const L_G = 6.969290134e-10 # rate of Terrestrial Time (TT) wrt Geocentric Coordinate Time (TCG)
const L_B = 1.550519768e-8  # rate of TDB wrt TCB
const TDB_0 = -65.5e-6daysec # initial offset (days)
const T_0 = 2443144.5003725 # Julian day
const one_plus_L_B_minus_L_G = 1 + L_B - L_G
