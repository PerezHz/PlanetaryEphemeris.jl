# JPL-DE-430/431 model for the orientation of the Earth

# From JPL DE430 documentation (Folkner et al., 2014):
# Only the long-term change of the Earth's orientation is modeled in the ephemeris
# integration. The Earth orientation model used for the DE 430 and 431 integration
# is based on the International Astronomical Union (IAU) 1976 precession model\* with
# an estimated linear correction and on a modified IAU 1980 nutation model \* including
# only terms with a period of 18.6 years.

export t2c_jpl_de430, c2t_jpl_de430, pole_rotation

@doc raw"""
    Ω(t)

Returns the longitude (in radians) of the mean ascending node of the lunar orbit on the 
ecliptic, measured from the mean equinox of date 
```math
\Omega(t) = 125^\circ 02' 40''.280 - \left(1934^\circ 8' 10''.539\right) T + 7''.455 T^2 + 0''.008 T^3,
```
where ``t = 36,525 T`` is the TDB time in Julian days from J2000.0.

See equation (5-64) in page 5-27 of https://doi.org/10.1002/0471728470.
"""
Ω(t) = deg2rad( (125+2/60+40.280/3600)-(1934+8/60+10.539/3600)*(t/36525)+(7.455/3600)*(t/36525)^2+(0.008/3600)*(t/36525)^3 )

@doc raw"""
    Delta_psi(t)

Returns the nutation in longitude (in radians) 
```math
\Delta\psi(t) = -17''.1996 \sin\Omega,
```
where ``t = 36,525 T`` is the TDB time in Julian days from J2000.0 and ``\Omega`` is the
longitude of the mean ascending node of the lunar orbit. 

See equation (5-190) in page 5-72 of https://doi.org/10.1002/0471728470.

See also [`Ω`](@ref).
"""
Delta_psi(t) = deg2rad( (-17.1996/3600)*sin(Ω(t)) )

@doc raw"""
    Delta_epsilon(t)

Returns the nutation in obliquity (in radians) 
```math
\Delta\epsilon(t) = 9''.2025 \cos\Omega,
```
where ``t = 36,525 T`` is the TDB time in Julian days from J2000.0 and ``\Omega`` is the
longitude of the mean ascending node of the lunar orbit.

See equation (5-190) in page 5-72 of https://doi.org/10.1002/0471728470.

See also [`Ω`](@ref).
"""
Delta_epsilon(t) = deg2rad( (9.2025/3600)*cos(Ω(t)) )

@doc raw"""
    pole_date(t)

Returns the true pole of date unit vector ``\mathbf{p}_\mathrm{d}``, computed by rotating the Earth-fixed pole vector by the effect of the 18.6-year nutation term. ``t`` is the TDB time in Julian days from J2000.0.

See equation (23) in page 11 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract.

See also [`nutation_iau80`](@ref).
"""
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

@doc raw"""
    ϵ̄(t)

Returns the mean obliquity (in radians) 
```math
\bar{\epsilon}(t) = 84,381''.448 - 46''.8150 T - 0''.00059 T^2 + 0''.001813 T^3,
```
where ``t = 36,525 T`` is the TDB time in Julian days from J2000.0.

See equation (5-153) in page 5-61 of https://doi.org/10.1002/0471728470.
"""
ϵ̄(t) = deg2rad( (84381.448/3600)-(46.815/3600)*(t/36525)-(0.00059/3600)*(t/36525)^2+(0.001813/3600)*(t/36525)^3 )

@doc raw"""
    pole_frame(t)
    
Returns the pole unit vector in the inertial frame ``\mathbf{p}_\mathrm{E}``, computed by 
precessing the pole of date with an estimated linear correction
```math
\mathbf{p}_\mathrm{E} = R_z(\zeta_A)R_y(-\theta_A)R_z(z_A)R_x(-\phi_x)R_y(-\phi_y)\mathbf{p}_\mathrm{d},
```
where ``t`` is the TDB time in Julian days from J2000.0, ``\mathbf{p}_d`` is the true pole 
of date unit vector, ``(\zeta_A, \theta_A, z_A)`` are the equatorial precession angles and 
``(\phi_x, \phi_y)`` are the linear corrections.

See equation (25) in page 11 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract.

See also [`pole_date`](@ref), [`Zeta`](@ref), [`Theta`](@ref), [`zeta`](@ref), [`phi_x`](@ref)
and [`phi_y`](@ref).
"""
function pole_frame(t)
    # True pole of date unit vector
    p_d = pole_date(t)
    # Inverse of precession matrix with linear corrections to precession angles
    pm = (Rz(Zeta(t))*Ry(-Theta(t))*Rz(zeta(t))*Rx(-phi_x(t))*Ry(-phi_y(t)))
    return pm*p_d

end

# where

@doc raw"""
    phi_x(t)

Returns the X-axis linear correction to precession (in radians)
```math
\phi_x = \phi_{x0} + 100T\times \frac{d\phi_x}{dt},
```
where ``t = 36,525 T`` is the TDB time in Julian days from J2000.0, ``\phi_{x0}`` is the 
X–axis rotation at J2000.0 (in arcseconds) and ``\frac{d\phi_x}{dt}`` is the negative 
obliquity rate correction (in arcseconds per year).

See equation (25) in page 11 and Table 10 in page 50 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract.

See also [`pole_frame`](@ref).
"""
phi_x(t) = deg2rad( (phi_x0 + (t/yr)*Dt_phi_x)/3600 )

phi_x0 = 5.6754203322893470E-03     # x-axis rotation at J2000.0 (arcseconds)
Dt_phi_x = 2.7689915574483550E-04   # Negative obliquity rate correction (arcseconds/year)

@doc raw"""
    phi_y(t)

Returns the Y-axis linear correction to precession (in radians)
```math
\phi_y = \phi_{y0} + 100T\times \frac{d\phi_y}{dt},
```
where ``t = 36,525 T`` is the TDB time in Julian days from J2000.0, ``\phi_{y0}`` is the 
Y–axis rotation at J2000.0 (in arcseconds) and ``\frac{d\phi_y}{dt}`` is precession rate 
correction times sine of obliquity (in arcseconds per year).

See equation (25) in page 11 and Table 10 in page 50 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract.

See also [`pole_frame`](@ref).
"""
phi_y(t) = deg2rad( (phi_y0 + (t/yr)*Dt_phi_y)/3600 )

phi_y0 = -1.7022656914989530E-02     # y-axis rotation at J2000.0 (arcseconds)
Dt_phi_y = -1.2118591216559240E-03   # Precession rate correction times sine of obliquity (arcseconds/year)

@doc raw"""
    Zeta(t)
    
Returns the ``\zeta_A`` equatorial precession angle (in radians)
```math
\zeta_A(t) = 2306''.2181 T + 0''.30188 T^2 + 0''.017998 T^3,
```
where ``t = 36,525 T`` is the TDB time in Julian days from J2000.0.

See equation (5-143) in page 5-58 of https://doi.org/10.1002/0471728470.

See also [`Theta`](@ref) and [`zeta`](@ref).
"""
function Zeta(t)
    t_cy = t/(100yr)
    Zeta_arcsec = 0.017998*t_cy
    Zeta_arcsec = (0.30188 + Zeta_arcsec)*t_cy
    Zeta_arcsec = (2306.2181 + Zeta_arcsec)*t_cy
    return deg2rad(Zeta_arcsec/3600)
end

@doc raw"""
    Theta(t)
    
Returns the ``\theta_A`` equatorial precession angle (in radians)
```math
\theta_A(t) = 2004''.3109 T - 0''.42665 T^2 - 0''.041833 T^3,
```
where ``t = 36,525 T`` is the TDB time in Julian days from J2000.0.

See equation (5-143) in page 5-58 of https://doi.org/10.1002/0471728470.

See also [`Zeta`](@ref) and [`zeta`](@ref).
"""
function Theta(t)
    t_cy = t/(100yr)
    Theta_arcsec = -0.041833*t_cy
    Theta_arcsec = (-0.42665 + Theta_arcsec)*t_cy
    Theta_arcsec = (2004.3109 + Theta_arcsec)*t_cy
    return deg2rad(Theta_arcsec/3600)
end

@doc raw"""
    zeta(t)
    
Returns the ``z_A`` equatorial precession angle (in radians)
```math
z_A(t) = 2306''.2181 T + 1''.09468 T^2 + 0''.018203 T^3,
```
where ``t = 36,525 T`` is the TDB time in Julian days from J2000.0.

See equation (5-143) in page 5-58 of https://doi.org/10.1002/0471728470.

See also [`Zeta`](@ref) and [`Theta`](@ref).
"""
function zeta(t)
    t_cy = t/(100yr)
    zeta_arcsec = 0.018203*t_cy
    zeta_arcsec = (1.09468 + zeta_arcsec)*t_cy
    zeta_arcsec = (2306.2181 + zeta_arcsec)*t_cy
    return deg2rad(zeta_arcsec/3600)
end

# The rotation matrices are defined by

@doc raw"""
    Rx(alpha::T) where {T<:Number}

Returns the rotation matrix around the x-axis
```math
R_x(\alpha) = 
\left[
\begin{array}{ccc}
    1 & 0 & 0 \\
    0 & \cos\alpha & \sin\alpha \\
    0 & -\sin\alpha & \cos\alpha \\    
\end{array}
\right],
```
where ``\alpha`` is an angle in radians. 

See equation (11) in page 9 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract.

See also [`Ry`](@ref) and [`Rz`](@ref).
"""
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

@doc raw"""
    Ry(alpha::T) where {T<:Number}

Returns the rotation matrix around the y-axis
```math
R_y(\alpha) = 
\left[
\begin{array}{ccc}
    \cos\alpha & 0 & -\sin\alpha \\
    0 & 1 & 0 \\
    \sin\alpha & 0 & \cos\alpha \\    
\end{array}
\right],
```
where ``\alpha`` is an angle in radians. 

See equation (12) in page 9 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract.

``\textbf{Note:}`` Lieske et al. (1977, 1979) introduces the convention that this matrix 
has opposite signs outside the sine terms (with respect to `Rx`, `Rz`). See equation (4) 
in page 3 of https://ui.adsabs.harvard.edu/abs/1977A%26A....58....1L/abstract and equation 
(3) in page 283 of https://ui.adsabs.harvard.edu/abs/1979A%26A....73..282L/abstract.

See also [`Rx`](@ref) and [`Rz`](@ref).
"""
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

@doc raw"""
    Rz(alpha::T) where {T<:Number}

Returns the rotation matrix around the z-axis
```math
R_z(\alpha) = 
\left[
\begin{array}{ccc}
    \cos\alpha & \sin\alpha & 0 \\
    -\sin\alpha & \cos\alpha & 0 \\
    0 & 0 & 1 \\    
\end{array}
\right],
```
where ``\alpha`` is an angle in radians. 

See equation (13) in page 9 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract.

See also [`Rx`](@ref) and [`Ry`](@ref).
"""
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

@doc raw"""
    pole_ra(t)

Returns the right ascension (in radians) of the pole unit vector in the inertial frame ``\mathbf{p}_\mathrm{E}(t)``
```math
\alpha = \arctan\left(\frac{p_\mathrm{Ey}(t)}{p_\mathrm{Ex}(t)}\right) + \pi,
```
where ``t`` is the TDB time in Julian days from J2000.0.

See also [`pole_frame`](@ref), [`pole_dec`](@ref) and [`pole_radec`](@ref).
"""
function pole_ra(t)
    pole_frame_t = pole_frame(t)
    return atan( pole_frame_t[2]/pole_frame_t[1] ) + π
end

@doc raw"""
    pole_dec(t)

Returns the declination (in radians) of the pole unit vector in the inertial frame ``\mathbf{p}_\mathrm{E}(t)``
```math
\delta = \arctan\left(\frac{p_\mathrm{Ez}(t)}{\sqrt{p_\mathrm{Ex}^2(t) + p_\mathrm{Ey}^2(t)}}\right),
```
where ``t`` is the TDB time in Julian days from J2000.0.

See also [`pole_frame`](@ref), [`pole_ra`](@ref) and [`pole_radec`](@ref).
"""
function pole_dec(t)
    pole_frame_t = pole_frame(t)
    return atan(  pole_frame_t[3]/sqrt(pole_frame_t[1]^2+pole_frame_t[2]^2)  )
end

@doc raw"""
    pole_radec(t)

Returns the right ascension and declination (both in radians) of the pole unit vector in the inertial frame ``\mathbf{p}_\mathrm{E}(t)``
```math
\begin{align*}
    \alpha & = \arctan\left(\frac{p_\mathrm{Ey}(t)}{p_\mathrm{Ex}(t)}\right) + \pi \\
    \delta & = \arctan\left(\frac{p_\mathrm{Ez}(t)}{\sqrt{p_\mathrm{Ex}^2(t) + p_\mathrm{Ey}^2(t)}}\right),
\end{align*}
```
where ``t`` is the TDB time in Julian days from J2000.0.

See also [`pole_frame`](@ref), [`pole_ra`](@ref) and [`pole_dec`](@ref).
"""
function pole_radec(t)
    pole_frame_t = pole_frame(t)
    pole_ra = atan( pole_frame_t[2]/pole_frame_t[1] ) + π
    pole_dec = atan(  pole_frame_t[3]/sqrt(pole_frame_t[1]^2+pole_frame_t[2]^2)  )
    return pole_ra, pole_dec
end

@doc raw"""
    pole_rotation(α::T, δ::T, W::T=zero(α)) where {T <: Number}

Returns the rotation matrix from the inertial frame to the frame with pole at right 
ascension ``\alpha``, declination ``\delta`` and prime meridian at ``W``
```math
A = R_z(W + \Delta W)R_x\left(\frac{\pi}{2} - \delta - \Delta\delta\right)R_z\left(\alpha + \Delta\alpha + \frac{\pi}{2}\right).
```

See equation (6-3) in page 6-5 of https://doi.org/10.1002/0471728470.

``\textbf{Note:}`` Rotation matrices ``(R_x, R_y, R_z)`` convention is the same as 
equations (11), (12) and (13) in page 9 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract.

See also [`Rx`](@ref) and [`Rz`](@ref).
"""
function pole_rotation(α::T, δ::T, W::T=zero(α)) where {T <: Number}
    return Rz(W)*Rx(π/2-δ)*Rz(π/2+α)
end

@doc raw"""
    earth_pole_rotation(t)

Returns the rotation matrix from inertial frame to Earth pole at time t (days) since J2000.0.

See also [`pole_radec`](@ref) and [`pole_rotation`](@ref).
"""
function earth_pole_rotation(t)
    # α_ep = Earth pole at time t (Julian days) since J2000.0
    # δ_ep = Earth pole at time t (Julian days) since J2000.0
    α_ep, δ_ep = pole_radec(t)
    return pole_rotation(α_ep, δ_ep)
end

@doc raw"""
    nutation_iau80(t)

Returns the IAU 1980 nutation matrix (Explanatory Supplement to the Astronomical Almanac 1992)
```math
N = R_x(-\bar{\epsilon}(t) - \Delta\epsilon(t))R_z(-\Delta\psi(t))R_x(\bar{\epsilon}(t)),
```
where ``t`` is the TDB time in Julian days from J2000.0, ``\bar{\epsilon}`` is the mean 
obliquity, ``\Delta \psi`` is the nutation in longitude and ``\Delta\epsilon`` is the 
nutation in obliquity. All angles are in radians. 

See equation (5-152) in page 5-60 of https://doi.org/10.1002/0471728470.

See also [`ϵ̄`](@ref), [`Delta_epsilon`](@ref) and [`Delta_psi`](@ref).
"""
function nutation_iau80(t)
    ϵ0 = ϵ̄(t)                # Mean obliquity (rad)
    Δϵ = Delta_epsilon(t)    # Nutation in obliquity (rad)
    Δψ = Delta_psi(t)        # Nutation in longitude (rad)
    ϵ = ϵ0 + Δϵ
    return Rx(-ϵ)*Rz(-Δψ)*Rx(ϵ0)
end

@doc raw"""
    t2c_jpl_de430(t)

Returns the matrix for terrestrial-to-celestial coordinate transformation, according to JPL DE 430/431 Earth orientation model
```math
\texttt{t2c_jpl_de430}(t) = \texttt{c2t_jpl_de430}^T(t),
```
where ``t`` is the TDB time in Julian days from J2000.0.

See also [`c2t_jpl_de430`](@ref).

"""
function t2c_jpl_de430(t)
    # inv_P_iau7680 = Rz(Zeta(t))*Ry(-Theta(t))*Rz(zeta(t))
    # corrections = Rx(-phi_x(t))*Ry(-phi_y(t))
    # # inv_N_iau80 = inv(nutation_iau80(t))
    # inv_N_iau80 = transpose(nutation_iau80(t))
    # return inv_P_iau7680*corrections*inv_N_iau80 #corrections*inv_P_iau7680*inv_N_iau80
    return transpose(c2t_jpl_de430(t))
end

@doc raw"""
    c2t_jpl_de430(t)

Returns the matrix for celestial-to-terrestrial coordinate transformation, according to 
JPL DE 430/431 Earth orientation model
```math
\texttt{c2t_jpl_de430}(t) = A\times C\times N,
```
where ``t`` is the TDB time in Julian days from J2000.0, ``A = R_z(-z_A(t))R_y(\theta_A(t))R_z(-\zeta_A(t))``
is the precession matrix, ``C = R_y(\phi_y(t))R_x(\phi_x(t))`` is the linear corrections matrix
and ``N`` is the nutation matrix.

See equation (5-147) in page (5-59) and equation (5-152) in page (5-60) of https://doi.org/10.1002/0471728470.
Also see equation (25) in page 11 and Table 10 in page 50 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract.

See also [`Rx`](@ref), [`Ry`](@ref), [`Rz`](@ref) and [`nutation_iau80`](@ref).
"""
function c2t_jpl_de430(t)
    # Precession matrix, see equation (5-147) in page (5-59) of https://doi.org/10.1002/0471728470
    P_iau7680 = Rz(-zeta(t))*Ry(Theta(t))*Rz(-Zeta(t))
    # Linear corrections to precession, see equation (25) in page 11 and Table 10 in page 50 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract
    corrections = Ry(phi_y(t))*Rx(phi_x(t))
    # Nutation matrix, see equation (5-152) in page 5-60 of https://doi.org/10.1002/0471728470
    N_iau80 = nutation_iau80(t)
    return N_iau80*corrections*P_iau7680
end

@doc raw"""
    moon_omega(ϕ::Taylor1, θ::Taylor1, ψ::Taylor1)
    
Returns the Moon's angular velocity, computed by differentiating the Euler angles 
``(\phi, \theta, \psi)``
```math
\begin{align*}
    \omega_x & = \dot{\phi}\sin\theta\sin\psi + \dot{\theta}\cos\psi \\
    \omega_y & = \dot{\phi}\sin\theta\cos\psi - \dot{\theta}\sin\psi \\
    \omega_z & = \dot{\phi}\cos\theta + \dot{\psi}.
\end{align*}
```
"""
function moon_omega(ϕ::Taylor1, θ::Taylor1, ψ::Taylor1)
    dϕ = differentiate(ϕ)
    dθ = differentiate(θ)
    dψ = differentiate(ψ)
    ωx = dϕ*sin(θ)*sin(ψ)+dθ*cos(ψ)
    ωy = dϕ*sin(θ)*cos(ψ)-dθ*sin(ψ)
    ωz = dψ+dϕ*cos(θ)
    return [ωx, ωy, ωz]
end

@doc raw"""
    ITM1(x::T, y::T, z::T) where {T <: Number}

Returns the first term of the time-dependent part of lunar total moment of intertia
```math
-\frac{k_{2,M} m_E R_M^5}{r^5}
\left[
\begin{array}{ccc}
    x^2 - \frac{1}{3}r^2 & xy & xz \\
    xy & y^2 - \frac{1}{3}r^2 & yz \\
    xz & yz & z^2 - \frac{1}{3}r^2 \\
\end{array}
\right],
```
where ``\mathbf{r} = (x, y, z)`` is the position of the Moon relative to Earth referred to 
the mantle frame, ``k_{2,M}`` is the lunar potential Love number, ``m_E`` is the mass of
the Earth and ``R_M`` is the equatorial radius of the Moon.

See the first term of equation (41) in page 17 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract.

See also [`ITM2`](@ref) and [`ITM`](@ref).
"""
function ITM1(x::T, y::T, z::T) where {T <: Number}
    r2 = x^2+y^2+z^2 # r^2
    r5 = r2^2.5 # r^5
    # Compute corresponding matrix elements
    m = Matrix{T}(undef, 3, 3)
    # Diagonal
    m[1,1] = x^2-r2/3
    m[2,2] = y^2-r2/3
    m[3,3] = z^2-r2/3
    # Off-diagonal (the matrix is symmetric)
    m[1,2] = x*y
    m[2,1] = m[1,2]
    m[1,3] = x*z
    m[3,1] = m[1,3]
    m[2,3] = y*z
    m[3,2] = m[2,3]
    return (-(k_2M*μ[ea]*R_moon^5)/r5)*m
end

@doc raw"""
    ITM2(ωx::T, ωy::T, ωz::T) where {T <: Number}

Returns the second term of the time-dependent part of lunar total moment of intertia
```math
\frac{k_{2,M} R_M^5}{3G}
\left[
\begin{array}{ccc}
    \omega_{m,x}^2 - \frac{1}{3}(\omega_m^2 - n^2) & \omega_{m,x}\omega_{m,y} & \omega_{m,x}\omega_{m,z} \\
    \omega_{m,x}\omega_{m,y} & \omega_{m,y}^2 - \frac{1}{3}(\omega_m^2 - n^2) & \omega_{m,y}\omega_{m,z} \\
    \omega_{m,x}\omega_{m,z} & \omega_{m,y}\omega_{m,z} & \omega_{m,z}^2 - \frac{1}{3}(\omega_m^2 + 2n^2) \\
\end{array}
\right],
```
where ``\mathbf{\omega}_m = (\omega_{m,x}, \omega_{m,y}, \omega_{m,z})`` is the angular 
velocity of the mantle in the mantle frame, ``k_{2,M}`` is the lunar potential Love number, 
``R_M`` is the equatorial radius of the Moon and ``n`` is the lunar mean motion.

``\textbf{Note:}`` Euler angles must be evaluated at time ``t-τ_M``.

See the second term of equation (41) in page 17 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract.

See also [`ITM1`](@ref) and [`ITM`](@ref).
"""
function ITM2(ωx::T, ωy::T, ωz::T) where {T <: Number}
    # Compute corresponding matrix elements
    m = Matrix{T}(undef, 3, 3)
    ωx2 = ωx*ωx       
    ωy2 = ωy*ωy
    ωz2 = ωz*ωz
    # Angular speed
    ω2 = ωx2+ωy2+ωz2
    aux = (ω2-n_moon^2)/3
    # Diagonal
    m[1,1] = ωx2-aux
    m[2,2] = ωy2-aux
    m[3,3] = ωz2-(ω2+2n_moon^2)/3
    # Off diagonal (the matrix is symmetric)
    m[1,2] = ωx*ωy
    m[2,1] = m[1,2]
    m[1,3] = ωx*ωz
    m[3,1] = m[1,3]
    m[2,3] = ωy*ωz
    m[3,2] = m[2,3]
    return ((k_2M*R_moon^5)/3)*m
end

@doc raw"""
    ITM(q::Vector{T}, eulang::Vector{T}, ω_m::Vector{T}) where {T <: Number}

Returns lunar mantle inertia tensor
```math
\mathbf{I}_m(t) = \tilde{\mathbf{I}}_m 
-
\frac{k_{2,M} m_E R_M^5}{r^5}
\left[
\begin{array}{ccc}
    x^2 - \frac{1}{3}r^2 & xy & xz \\
    xy & y^2 - \frac{1}{3}r^2 & yz \\
    xz & yz & z^2 - \frac{1}{3}r^2 \\
\end{array}
\right]
+
\frac{k_{2,M} R_M^5}{3G}
\left[
\begin{array}{ccc}
    \omega_{m,x}^2 - \frac{1}{3}(\omega_m^2 - n^2) & \omega_{m,x}\omega_{m,y} & \omega_{m,x}\omega_{m,z} \\
    \omega_{m,x}\omega_{m,y} & \omega_{m,y}^2 - \frac{1}{3}(\omega_m^2 - n^2) & \omega_{m,y}\omega_{m,z} \\
    \omega_{m,x}\omega_{m,z} & \omega_{m,y}\omega_{m,z} & \omega_{m,z}^2 - \frac{1}{3}(\omega_m^2 + 2n^2) \\
\end{array}
\right],
```
where ``\mathbf{r} = (x, y, z)`` is the position of the Moon relative to Earth referred to the
mantle frame, ``k_{2,M}`` is the lunar potential Love number, ``m_E`` is the mass of the Earth,
``R_M`` is the equatorial radius of the Moon, ``\mathbf{\omega}_m = (\omega_{m,x}, \omega_{m,y}, \omega_{m,z})``
is the angular velocity of the mantle in the mantle frame and ``n`` is the lunar mean motion.

``\textbf{Note:}`` Euler angles must be evaluated at time ``t-τ_M``.

See equation (41) in page 17 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract.

See also [`ITM1`](@ref) and [`ITM2`](@ref).
"""
function ITM(q::Vector{T}, eulang::Vector{T}, ω_m::Vector{T}) where {T <: Number}
    # Coordinate transformation matrix: inertial frame -> lunar mantle frame
    M_del_mo = pole_rotation(eulang[1] - (pi/2), (pi/2) - eulang[2], eulang[3])
    # Transform delayed geocentric position of Moon (space-fixed->lunar mantle frame)
    X_me_del_τ_M = q[3mo-2:3mo] .- q[3ea-2:3ea]
    X_me_del_τ_M_lm = M_del_mo*X_me_del_τ_M
    # First term of the time-dependent part of lunar total moment of intertia
    itm1 = ITM1(X_me_del_τ_M_lm[1], X_me_del_τ_M_lm[2], X_me_del_τ_M_lm[3])
    # Second term of the time-dependent part of lunar total moment of intertia
    itm2 = ITM2(ω_m[1], ω_m[2], ω_m[3])
    # Lunar mantle intertia tensor
    return (ITM_und-I_c)*one(q[1]) + itm1 + itm2
end