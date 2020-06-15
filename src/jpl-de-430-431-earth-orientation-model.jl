# JPL-DE-430/431 model for the orientation of the Earth

# From JPL DE430 documentation (Folkner et al., 2014):
# Only the long-term change of the Earth's orientation is modeled in the ephemeris
# integration. The Earth orientation model used for the DE 430 and 431 integration
# is based on the International Astronomical Union (IAU) 1976 precession model\* with
# an estimated linear correction and on a modified IAU 1980 nutation model \* including
# only terms with a period of 18.6 years.

export t2c_jpl_de430, c2t_jpl_de430, pole_rotation

# The mean longitude of the ascending node of the lunar orbit measured on the ecliptic
# plane from the mean equinox of date is calculated by

Ω(t) = deg2rad( (125+2/60+40.280/3600)-(1934+8/60+10.539/3600)*(t/36525)+(7.455/3600)*(t/36525)^2+(0.008/3600)*(t/36525)^3 )

# where `t` is the TDB time in Julian days from J2000.0. The nutations in longitude
# $\Delta \psi$ and obliquity $\Delta \epsilon$ are given by

#(these expressions coincide with Moyer, 2003, eq. 5-190, page 5-72)
Delta_psi(t) = deg2rad( (-17.1996/3600)*sin(Ω(t)) )
Delta_epsilon(t) = deg2rad( (9.2025/3600)*cos(Ω(t)) )

# The true pole of date unit vector $\vec p_\mathrm{d}$ is computed by rotating the
# Earth-fixed pole vector by the effect of the 18.6-year nutation term to give

function pole_date(t)
    # ϵ0 = ϵ̄(t)
    # Δϵ = Delta_epsilon(t)
    # Δψ = Delta_psi(t)
    # ϵ = ϵ0 + Δϵ
    # pdx = sin(Δψ)*sin(ϵ)
    # pdy = cos(Δψ)*sin(ϵ)*cos(ϵ0) - cos(ϵ)*sin(ϵ0)
    # pdz = cos(Δψ)*sin(ϵ)*sin(ϵ0) + cos(ϵ)*cos(ϵ0)
    Nut_mat = nutation_iau80(t)
    # inv(Nut_mat)*[0,0,1] == Nut_mat[3,:]
    return Nut_mat[3,:]
end

# where the mean obliquity $\bar \epsilon$ is given by

ϵ̄(t) = deg2rad( (84381.448/3600)-(46.815/3600)*(t/36525)-(0.00059/3600)*(t/36525)^2+(0.001813/3600)*(t/36525)^3 )

# The pole unit vector in the inertial frame $\vec p_\mathrm{E}$ is computed by
# precessing the pole of date with an estimated linear correction,

function pole_frame(t)
    p_d = pole_date(t)
    # inverse of precession matrix with linear corrections to precession angles
    pm = (Rz(Zeta(t))*Ry(-Theta(t))*Rz(zeta(t))*Rx(-phi_x(t))*Ry(-phi_y(t)))
    return pm*p_d

end

# where

phi_x(t) = deg2rad( (phi_x0 + (t/yr)*Dt_phi_x)/3600 )
phi_y(t) = deg2rad( (phi_y0 + (t/yr)*Dt_phi_y)/3600 )

phi_x0 = 5.6754203322893470E-03 #x-axis rotation at J2000.0 (arcseconds)
phi_y0 = -1.7022656914989530E-02 #y-axis rotation at J2000.0 (arcseconds)
Dt_phi_x = 2.7689915574483550E-04 #Negative obliquity rate correction (arcseconds/year)
Dt_phi_y = -1.2118591216559240E-03 #Precession rate correction times sine of obliquity (arcseconds/year)

# are estimated linear corrections with offsets and rates given by `phi_x0`, `phi_y0`,
# `Dt_phi_x` and `Dt_phi_y` and the precession angles are given by

# Zeta = t -> deg2rad( (2306.2181/3600)*(t/36525)+(0.30188/3600)*(t/36525)^2+(0.017998/3600)*(t/36525)^3 )
# Theta = t -> deg2rad( (2004.3109/3600)*(t/36525)-(0.42665/3600)*(t/36525)^2-(0.041833/3600)*(t/36525)^3 )
# zeta = t -> deg2rad( (2306.2181/3600)*(t/36525)+(1.09468/3600)*(t/36525)^2+(0.018203/3600)*(t/36525)^3 )

function Zeta(t)
    t_cy = t/(100yr)
    Zeta_arcsec = 0.017998*t_cy
    Zeta_arcsec = (0.30188 + Zeta_arcsec)*t_cy
    Zeta_arcsec = (2306.2181 + Zeta_arcsec)*t_cy
    return deg2rad(Zeta_arcsec/3600)
end

function Theta(t)
    t_cy = t/(100yr)
    Theta_arcsec = -0.041833*t_cy
    Theta_arcsec = (-0.42665 + Theta_arcsec)*t_cy
    Theta_arcsec = (2004.3109 + Theta_arcsec)*t_cy
    return deg2rad(Theta_arcsec/3600)
end

function zeta(t)
    t_cy = t/(100yr)
    zeta_arcsec = 0.018203*t_cy
    zeta_arcsec = (1.09468 + zeta_arcsec)*t_cy
    zeta_arcsec = (2306.2181 + zeta_arcsec)*t_cy
    return deg2rad(zeta_arcsec/3600)
end

# The rotation matrices are defined by

function Rx(alpha::T) where {T<:Number}
    res = Array{T}(undef, 3, 3)
    res[1, 1] = one(alpha)
    res[2, 1] = zero(alpha)
    res[3, 1] = zero(alpha)
    res[1, 2] = zero(alpha)
    res[2, 2] = cos(alpha)
    res[3, 2] = -sin(alpha)
    res[1, 3] = zero(alpha)
    res[2, 3] = sin(alpha)
    res[3, 3] = cos(alpha)
    return res
end

# Lieske et al. (1977, 1979) introduces the convention that this matrix has
# opposite signs outside the sine terms (wrt Rx, Rz)
function Ry(alpha::T) where {T<:Number}
    res = Array{T}(undef, 3, 3)
    res[1, 1] = cos(alpha)
    res[2, 1] = zero(alpha)
    res[3, 1] = sin(alpha)
    res[1, 2] = zero(alpha)
    res[2, 2] = one(alpha)
    res[3, 2] = zero(alpha)
    res[1, 3] = -sin(alpha)
    res[2, 3] = zero(alpha)
    res[3, 3] = cos(alpha)
    return res
end

function Rz(alpha::T) where {T<:Number}
    res = Array{T}(undef, 3, 3)
    res[1, 1] = cos(alpha)
    res[2, 1] = -sin(alpha)
    res[3, 1] = zero(alpha)
    res[1, 2] = sin(alpha)
    res[2, 2] = cos(alpha)
    res[3, 2] = zero(alpha)
    res[1, 3] = zero(alpha)
    res[2, 3] = zero(alpha)
    res[3, 3] = one(alpha)
    return res
end

# myatan(x, y) = y>=zero(x)?( x>=zero(x)?atan(y/x):(atan(y/x)+pi) ):( x>=zero(x)?(atan(y/x)+2pi):(atan(y/x)+pi) )
# myatan2(x, y) = y>=zero(x)?( x>=zero(x)?atan(y/x):(atan(y/x)-pi) ):( x>=zero(x)?(atan(y/x)):(atan(y/x)+pi) )

function pole_ra(t)
    pole_frame_t = pole_frame(t)
    return atan( pole_frame_t[2]/pole_frame_t[1] ) + π
end

function pole_dec(t)
    pole_frame_t = pole_frame(t)
    return atan(  pole_frame_t[3]/sqrt(pole_frame_t[1]^2+pole_frame_t[2]^2)  )
end

function pole_radec(t)
    pole_frame_t = pole_frame(t)
    pole_ra = atan( pole_frame_t[2]/pole_frame_t[1] ) + π
    pole_dec = atan(  pole_frame_t[3]/sqrt(pole_frame_t[1]^2+pole_frame_t[2]^2)  )
    return pole_ra, pole_dec
end

# rotation from inertial frame to frame with pole at right ascension α, declination δ and prime meridian at W
# taken from matrix A in Moyer (2003), eq. 6-3 (page 6-5)
# note that rotation matrices (Rx, Ry, Rz) convention is the same as in Folkner et al. (2014)
function pole_rotation(α::T, δ::T, W::T=zero(α)) where {T <: Number}
    return Rz(W)*Rx(π/2-δ)*Rz(π/2+α)
end

# rotation matrix from inertial frame to Earth pole at time t (days) since J2000.0
function earth_pole_rotation(t)
    # α_ep = Earth pole at time t (Julian days) since J2000.0
    # δ_ep = Earth pole at time t (Julian days) since J2000.0
    α_ep, δ_ep = pole_radec(t)
    return pole_rotation(α_ep, δ_ep)
end

# IAU 1980 nutation matrix (ESAA 1992)
# Δψ (rad): nutation in longitude
# ϵ (rad): mean obliquity
# ϵ0 (rad): true obliquity of date
# same as Eq. (5-152) in Moyer, 2003
function nutation_iau80(t)
    ϵ0 = ϵ̄(t)
    Δϵ = Delta_epsilon(t)
    Δψ = Delta_psi(t)
    ϵ = ϵ0 + Δϵ
    return Rx(-ϵ)*Rz(-Δψ)*Rx(ϵ0)
end

# terrestrial-to-celestial coordinate transformation, JPL DE 430/431 Earth orientation model
function t2c_jpl_de430(t)
    # inv_P_iau7680 = Rz(Zeta(t))*Ry(-Theta(t))*Rz(zeta(t))
    # corrections = Rx(-phi_x(t))*Ry(-phi_y(t))
    # # inv_N_iau80 = inv(nutation_iau80(t))
    # inv_N_iau80 = transpose(nutation_iau80(t))
    # return inv_P_iau7680*corrections*inv_N_iau80 #corrections*inv_P_iau7680*inv_N_iau80
    return transpose(c2t_jpl_de430(t))
end

# celestial-to-terrestrial coordinate transformation, JPL DE 430/431 Earth orientation model
function c2t_jpl_de430(t)
    P_iau7680 = Rz(-zeta(t))*Ry(Theta(t))*Rz(-Zeta(t)) # Moyer (2003), eq. 5-147 (page 5-59)
    corrections = Ry(phi_y(t))*Rx(phi_x(t))
    N_iau80 = nutation_iau80(t) # Moyer (2003), eq. 5-152 (page 5-60)
    return N_iau80*corrections*P_iau7680
end

function moon_omega(ϕ::Taylor1, θ::Taylor1, ψ::Taylor1)
    dϕ = differentiate(ϕ)
    dθ = differentiate(θ)
    dψ = differentiate(ψ)
    ωx = dϕ*sin(θ)*sin(ψ)+dθ*cos(ψ)
    ωy = dϕ*sin(θ)*cos(ψ)-dθ*sin(ψ)
    ωz = dψ+dϕ*cos(θ)
    return [ωx, ωy, ωz]
end

#first term of time-dependent part of lunar total moment of inertia (Folkner et al., 2014, eq. 41, 1st term)
function ITM1(x::T, dx::T, y::T, dy::T, z::T, dz::T) where {T <: Number}
    # evaluate lunar geocentric position at time t-τ_M
    xd = x-dx*τ_M # x(t-τ_M)
    yd = y-dy*τ_M # y(t-τ_M)
    zd = z-dz*τ_M # z(t-τ_M)
    rd2 = xd^2+yd^2+zd^2 # r(t-τ_M)^2
    rd5 = rd2^2.5 # r(t-τ_M)^5
    # compute corresponding matrix elements
    m = Matrix{T}(undef, 3, 3)
    m[1,1] = xd^2-rd2/3
    m[2,2] = yd^2-rd2/3
    m[3,3] = zd^2-rd2/3
    m[1,2] = xd*yd
    m[2,1] = m[1,2]
    m[1,3] = xd*zd
    m[3,1] = m[1,3]
    m[2,3] = yd*zd
    m[3,2] = m[2,3]
    return (-(k_2M*μ[ea]*R_moon^5)/rd5)*m
end

#second term of time-dependent part of lunar total moment of inertia (Folkner et al., 2014, eq. 41, 2nd term)
# note: Euler angles must be evaluated at time `t-τ_M`
function ITM2(ϕ::T, θ::T, ψ::T) where {T <: Number}
    m = Matrix{T}(undef, 3, 3)
    ω = moon_omega(ϕ, θ, ψ)
    ωx2 = ω[1]*ω[1]
    ωy2 = ω[2]*ω[2]
    ωz2 = ω[3]*ω[3]
    ω2 = ωx2+ωy2+ωz2
    aux = (ω2-n_moon^2)/3
    m[1,1] = ωx2-aux
    m[2,2] = ωy2-aux
    m[3,3] = ωz2-(ω2+2n_moon^2)/3
    m[1,2] = ω[1]*ω[2]
    m[2,1] = m[1,2]
    m[1,3] = ω[1]*ω[3]
    m[3,1] = m[1,3]
    m[2,3] = ω[2]*ω[3]
    m[3,2] = m[2,3]
    return ((k_2M*R_moon^5)/3)*m
end
