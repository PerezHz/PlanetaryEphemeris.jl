# This file is part of the PlanetaryEphemeris.jl package; MIT licensed

# JPL-DE-430/431 model for the orientation of the Earth

# From JPL DE430 documentation (Folkner et al., 2014):
# Only the long-term change of the Earth's orientation is modeled in the ephemeris
# integration. The Earth orientation model used for the DE 430 and 431 integration
# is based on the International Astronomical Union (IAU) 1976 precession model\* with
# an estimated linear correction and on a modified IAU 1980 nutation model \* including
# only terms with a period of 18.6 years.

@doc raw"""
    Ω(t)

Return the longitude (in radians) of the mean ascending node of the lunar orbit on the
ecliptic, measured from the mean equinox of date
```math
\Omega(t) = 125^\circ 02' 40''.280 - \left(1934^\circ 8' 10''.539\right) T + 7''.455 T^2 + 0''.008 T^3,
```
where ``t = 36,525 T`` is the TDB time in Julian days from J2000.0.

See equation (5-64) in page 5-27 of https://doi.org/10.1002/0471728470.
"""
Ω(t) = deg2rad( evaluate(Taylor1(coeffs_Ω_A, get_order(t)), t/36525) )

const coeffs_Ω_A = [125+2/60+40.280/3600, -(1934+8/60+10.539/3600), 7.455/3600, 0.008/3600]

@doc raw"""
    Ω!(res, t, auxs)

Inplace version of ``Ω(t)``; see [`Ω`](@ref) for details.
"""
function Ω!(res::Taylor1{T}, t::Taylor1{T}, auxs) where {T<:Real}
    t_cy = auxs[1]
    auxhorner = auxs[2] # aux for _horner!
    Ω_A = auxs[3]  # Ω_A
    for i in eachindex(coeffs_Ω_A)
        Ω_A.coeffs[i] = coeffs_Ω_A[i]
    end
    for ord in eachindex(t)
        res[ord] = zero(t[ord])
        t_cy[ord] = t[ord]/36525 # t/36525
        auxhorner[ord] = zero(t[ord])
    end
    TaylorSeries._horner!(res, Ω_A, t_cy, auxhorner) # Ω_A(t_cy)
    for ord in eachindex(t)
        res[ord] = res[ord]*(convert(float(T), pi) / 180)#deg2rad(res[ord])
    end
    return nothing
end

@doc raw"""
    Delta_psi(t)

Return the nutation in longitude (in radians)
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
    Delta_psi!(res, t, auxs)

Inplace version of ``Delta_psi(t)``; see [`Delta_psi`](@ref) for details.
"""
function Delta_psi!(deltapsi::Taylor1{T}, t::Taylor1{T}, auxs) where {T<:Real}
    numfactor = -17.1996/3600
    omega_t_cy = auxs[1]
    Ω!(omega_t_cy, t, view(auxs, 2:4))
    sinO = auxs[5]
    cosO = auxs[6]
    for ord in eachindex(t)
        TaylorSeries.sincos!(sinO, cosO, omega_t_cy, ord)
        deltapsi[ord] = numfactor*sinO[ord]*(convert(float(T), pi) / 180)#deg2rad( numfactor*sinO[ord] )
    end
    return nothing
end

@doc raw"""
    Delta_epsilon(t)

Return the nutation in obliquity (in radians)
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
    Delta_epsilon!(res, t, auxs)

Inplace version of ``Delta_epsilon(t)``; see [`Delta_epsilon`](@ref) for details.
"""
function Delta_epsilon!(deltaepsilon::Taylor1{T}, t::Taylor1{T}, auxs) where {T<:Real}
    numfactor = 9.2025/3600
    omega_t_cy = auxs[1]
    Ω!(omega_t_cy, t, view(auxs, 2:4))
    sinO = auxs[5]
    cosO = auxs[6]
    for ord in eachindex(t)
        TaylorSeries.sincos!(sinO, cosO, omega_t_cy, ord)
        deltaepsilon[ord] = numfactor*cosO[ord]*(convert(float(T), pi) / 180)#deg2rad( numfactor*cosO[ord] )
    end
    return nothing
end

@doc raw"""
    pole_date(t)

Return the true pole of date unit vector ``\mathbf{p}_\mathrm{d}``, computed by rotating the Earth-fixed pole vector by the effect of the 18.6-year nutation term. ``t`` is the TDB time in Julian days from J2000.0.

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

Return the mean obliquity (in radians)
```math
\bar{\epsilon}(t) = 84,381''.448 - 46''.8150 T - 0''.00059 T^2 + 0''.001813 T^3,
```
where ``t = 36,525 T`` is the TDB time in Julian days from J2000.0.

See equation (5-153) in page 5-61 of https://doi.org/10.1002/0471728470.
"""
ϵ̄(t) = deg2rad( evaluate(Taylor1(coeffs_ϵ̄, get_order(t)), t/36525) )

const coeffs_ϵ̄ = [84381.448/3600, -46.815/3600, -0.00059/3600, 0.001813/3600]

@doc raw"""
    ϵ̄!(res, t, auxs)

Inplace version of ``ϵ̄(t)``; see [`ϵ̄`](@ref) for details.
"""
function ϵ̄!(res::Taylor1{T}, t::Taylor1{T}, auxs) where {T<:Real}
    t_cy = auxs[1]
    auxhorner = auxs[2] # aux for _horner!
    ϵ̄_A = auxs[3]  # ϵ̄_A
    for i in eachindex(coeffs_ϵ̄)
        ϵ̄_A.coeffs[i] = coeffs_ϵ̄[i]
    end
    for ord in eachindex(t)
        res[ord] = zero(t[ord])
        t_cy[ord] = t[ord]/36525 # t/(36525)
        auxhorner[ord] = zero(t[ord])
    end
    TaylorSeries._horner!(res, ϵ̄_A, t_cy, auxhorner) # ϵ̄_A(t_cy)
    for ord in eachindex(t)
        res[ord] = res[ord]*(convert(float(T), pi) / 180)#deg2rad(res[ord])
    end
    return nothing
end

@doc raw"""
    pole_frame(t)

Return the pole unit vector in the inertial frame ``\mathbf{p}_\mathrm{E}``, computed by
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

Return the X-axis linear correction to precession (in radians)
```math
\phi_x = \phi_{x0} + 100T\times \frac{d\phi_x}{dt},
```
where ``t = 36,525 T`` is the TDB time in Julian days from J2000.0, ``\phi_{x0}`` is the
X–axis rotation at J2000.0 (in arcseconds) and ``\frac{d\phi_x}{dt}`` is the negative
obliquity rate correction (in arcseconds per year).

See equation (25) in page 11 and Table 10 in page 50 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract.

See also [`pole_frame`](@ref).
"""
function phi_x(t)
    phix = zero(t)
    for ord in eachindex(t)
        phix[ord] = (t[ord]/yr)*Dt_phi_x/3600
    end
    phix[0] += phi_x0
    for ord in eachindex(t)
        phix[ord] = deg2rad( phix[ord] )
    end
    return phix
end

@doc raw"""
    phi_x!(res, t)

Inplace version of ``phi_x(t)``; see [`phi_x`](@ref) for details.
"""
function phi_x!(phix::Taylor1{T}, t::Taylor1{T}) where {T<:Real}
    for ord in eachindex(t)
        phix[ord] = (t[ord]/yr)*Dt_phi_x/3600
    end
    phix[0] += phi_x0
    for ord in eachindex(t)
        phix[ord] = phix[ord]*(convert(float(T), pi) / 180)#deg2rad( phix[ord] )
    end
    return nothing
end


const phi_x0 = 5.6754203322893470E-03     # x-axis rotation at J2000.0 (arcseconds)
const Dt_phi_x = 2.7689915574483550E-04   # Negative obliquity rate correction (arcseconds/year)

@doc raw"""
    phi_y(t)

Return the Y-axis linear correction to precession (in radians)
```math
\phi_y = \phi_{y0} + 100T\times \frac{d\phi_y}{dt},
```
where ``t = 36,525 T`` is the TDB time in Julian days from J2000.0, ``\phi_{y0}`` is the
Y–axis rotation at J2000.0 (in arcseconds) and ``\frac{d\phi_y}{dt}`` is precession rate
correction times sine of obliquity (in arcseconds per year).

See equation (25) in page 11 and Table 10 in page 50 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract.

See also [`pole_frame`](@ref).
"""
function phi_y(t)
    phiy = zero(t)
    for ord in eachindex(t)
        phiy[ord] = (t[ord]/yr)*Dt_phi_y/3600
    end
    phiy[0] += phi_y0
    for ord in eachindex(t)
        phiy[ord] = deg2rad( phiy[ord] )
    end
    return phiy
end

@doc raw"""
    phi_y!(res, t)

Inplace version of ``phi_y(t)``; see [`phi_y`](@ref) for details.
"""
function phi_y!(phiy::Taylor1{T}, t::Taylor1{T}) where {T<:Real}
    for ord in eachindex(t)
        phiy[ord] = (t[ord]/yr)*Dt_phi_y/3600
    end
    phiy[0] += phi_y0
    for ord in eachindex(t)
        phiy[ord] = phiy[ord]*(convert(float(T), pi) / 180)#deg2rad( phiy[ord] )
    end
    return nothing
end

const phi_y0 = -1.7022656914989530E-02     # y-axis rotation at J2000.0 (arcseconds)
const Dt_phi_y = -1.2118591216559240E-03   # Precession rate correction times sine of obliquity (arcseconds/year)

@doc raw"""
    Zeta(t)

Return the ``\zeta_A`` equatorial precession angle (in radians)
```math
\zeta_A(t) = 2306''.2181 T + 0''.30188 T^2 + 0''.017998 T^3,
```
where ``t = 36,525 T`` is the TDB time in Julian days from J2000.0.

See equation (5-143) in page 5-58 of https://doi.org/10.1002/0471728470.

See also [`Theta`](@ref) and [`zeta`](@ref).
"""
function Zeta(t)
    t_cy = t/(100yr)
    # Zeta_arcsec = 0.017998*t_cy
    # Zeta_arcsec = (0.30188 + Zeta_arcsec)*t_cy
    # Zeta_arcsec = (2306.2181 + Zeta_arcsec)*t_cy
    ZetaA = Taylor1(coeffs_Zeta, get_order(t))
    Zeta_arcsec = ZetaA(t_cy)
    return deg2rad(Zeta_arcsec/3600)
end

@doc raw"""
    Zeta!(res, t, auxs)

Inplace version of ``Zeta(t)``; see [`Zeta`](@ref) for details.
"""
function Zeta!(res::Taylor1{T}, t::Taylor1{T}, auxs) where {T<:Real}
    t_cy = auxs[1]
    auxhorner = auxs[2] # aux for _horner!
    Zeta_A = auxs[3] # Zeta_A
    for i in eachindex(coeffs_Zeta)
        Zeta_A.coeffs[i] = coeffs_Zeta[i]
    end
    hundyrs = 100yr
    for ord in eachindex(t)
        res[ord] = zero(t[ord])
        t_cy[ord] = t[ord]/hundyrs # t/(100yr)
        auxhorner[ord] = zero(t[ord])
    end
    TaylorSeries._horner!(res, Zeta_A, t_cy, auxhorner) # Zeta_A(t_cy)
    for ord in eachindex(t)
        res[ord] = (res[ord]/3600)*(convert(float(T), pi) / 180)#deg2rad(res[ord]/3600)
    end
    return nothing
end

const coeffs_Zeta = [0.0, 2306.2181, 0.30188, 0.017998]

@doc raw"""
    Theta(t)

Return the ``\theta_A`` equatorial precession angle (in radians)
```math
\theta_A(t) = 2004''.3109 T - 0''.42665 T^2 - 0''.041833 T^3,
```
where ``t = 36,525 T`` is the TDB time in Julian days from J2000.0.

See equation (5-143) in page 5-58 of https://doi.org/10.1002/0471728470.

See also [`Zeta`](@ref) and [`zeta`](@ref).
"""
function Theta(t)
    t_cy = t/(100yr)
    # Theta_arcsec = -0.041833*t_cy
    # Theta_arcsec = (-0.42665 + Theta_arcsec)*t_cy
    # Theta_arcsec = (2004.3109 + Theta_arcsec)*t_cy
    theta_A = Taylor1(coeffs_Theta, get_order(t))
    Theta_arcsec = theta_A(t_cy)
    return deg2rad(Theta_arcsec/3600)
end

@doc raw"""
    Theta!(res, t, auxs)

Inplace version of ``Theta(t)``; see [`Theta`](@ref) for details.
"""
function Theta!(res::Taylor1{T}, t::Taylor1{T}, auxs) where {T<:Real}
    t_cy = auxs[1]
    auxhorner = auxs[2]
    theta_A = auxs[3]  # Theta_A
    for i in eachindex(coeffs_Theta)
        theta_A.coeffs[i] = coeffs_Theta[i]
    end
    hundyrs = 100yr
    for ord in eachindex(t)
        res[ord] = zero(t[ord])
        t_cy[ord] = t[ord]/hundyrs # t/(100yr)
        auxhorner[ord] = zero(t[ord]) # aux for _horner!
    end
    TaylorSeries._horner!(res, theta_A, t_cy, auxhorner) # Theta_A(t_cy)
    for ord in eachindex(t)
        res[ord] = (res[ord]/3600)*(convert(float(T), pi) / 180)#deg2rad(res[ord]/3600)
    end
    return nothing
end

const coeffs_Theta = [0.0, 2004.3109, -0.42665, -0.041833]


@doc raw"""
    zeta(t)

Return the ``z_A`` equatorial precession angle (in radians)
```math
z_A(t) = 2306''.2181 T + 1''.09468 T^2 + 0''.018203 T^3,
```
where ``t = 36,525 T`` is the TDB time in Julian days from J2000.0.

See equation (5-143) in page 5-58 of https://doi.org/10.1002/0471728470.

See also [`Zeta`](@ref) and [`Theta`](@ref).
"""
function zeta(t)
    t_cy = t/(100yr)
    # zeta_arcsec = 0.018203*t_cy
    # zeta_arcsec = (1.09468 + zeta_arcsec)*t_cy
    # zeta_arcsec = (2306.2181 + zeta_arcsec)*t_cy
    zeta_A = Taylor1(coeffs_zeta, get_order(t))
    zeta_arcsec = zeta_A(t_cy)
    return deg2rad(zeta_arcsec/3600)
end

@doc raw"""
    zeta!(res, t, auxs)

Inplace version of ``zeta(t)``; see [`zeta`](@ref) for details.
"""
function zeta!(res::Taylor1{T}, t::Taylor1{T}, auxs) where {T<:Real}
    t_cy = auxs[1]
    auxhorner = auxs[2] # aux for _horner!
    zeta_A = auxs[3] # zeta_A
    for i in eachindex(coeffs_zeta)
        zeta_A.coeffs[i] = coeffs_zeta[i]
    end
    hundyrs = 100yr
    for ord in eachindex(t)
        res[ord] = zero(t[ord])
        t_cy[ord] = t[ord]/hundyrs # t/(100yr)
        auxhorner[ord] = zero(t[ord])
    end
    TaylorSeries._horner!(res, zeta_A, t_cy, auxhorner) # zeta_A(t_cy)
    for ord in eachindex(t)
        res[ord] = (res[ord]/3600)*(convert(float(T), pi) / 180)#deg2rad(res[ord]/3600)
    end
    return nothing
end

const coeffs_zeta = [0.0, 2306.2181, 1.09468, 0.018203]


# The rotation matrices are defined by

@doc raw"""
    Rx(alpha::T) where {T<:Number}

Return the rotation matrix around the x-axis
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
    # Allocate memory
    res = Array{T}(undef, 3, 3)
    # Local variables
    one_alpha = one(alpha)
    zero_alpha = zero(alpha)
    sin_alpha, cos_alpha = sincos(alpha)
    # Matrix elements
    res[1, 1] = one_alpha
    res[2, 1] = zero_alpha
    res[3, 1] = zero_alpha
    res[1, 2] = zero_alpha
    res[2, 2] = cos_alpha
    res[3, 2] = -sin_alpha
    res[1, 3] = zero_alpha
    res[2, 3] = sin_alpha
    res[3, 3] = cos_alpha
    return res
end

@doc raw"""
    Rx!(res, alpha, auxsR)

Inplace version of ``Rx(alpha)``; see [`Rx`](@ref) for details.
"""
function Rx!(res::Matrix{Taylor1{T}}, alpha::Taylor1{T}, auxsR) where {T<:Real}
    # Local variables
    one_alpha = auxsR[1] # one(alpha)
    zero_alpha = auxsR[2] # zero(alpha)
    sin_alpha = auxsR[3] # sin(alpha)
    cos_alpha = auxsR[4] # cos(alpha)
    for ord in eachindex(alpha)
        TaylorSeries.sincos!(sin_alpha, cos_alpha, alpha, ord)
        # Matrix elements
        res[1, 1][ord] = one_alpha[ord]
        res[2, 1][ord] = zero_alpha[ord]
        res[3, 1][ord] = zero_alpha[ord]
        res[1, 2][ord] = zero_alpha[ord]
        res[2, 2][ord] = cos_alpha[ord]
        res[3, 2][ord] = -sin_alpha[ord]
        res[1, 3][ord] = zero_alpha[ord]
        res[2, 3][ord] = sin_alpha[ord]
        res[3, 3][ord] = cos_alpha[ord]
    end
    return nothing
end

@doc raw"""
    Ry(alpha::T) where {T<:Number}

Return the rotation matrix around the y-axis
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
    # Allocate memory
    res = Array{T}(undef, 3, 3)
    # Local variables
    one_alpha = one(alpha)
    zero_alpha = zero(alpha)
    sin_alpha, cos_alpha = sincos(alpha)
    # Matrix elements
    res[1, 1] = cos_alpha
    res[2, 1] = zero_alpha
    res[3, 1] = sin_alpha
    res[1, 2] = zero_alpha
    res[2, 2] = one_alpha
    res[3, 2] = zero_alpha
    res[1, 3] = -sin_alpha
    res[2, 3] = zero_alpha
    res[3, 3] = cos_alpha
    return res
end

@doc raw"""
    Ry!(res, alpha, auxsR)

Inplace version of ``Ry(alpha)``; see [`Ry`](@ref) for details.
"""
function Ry!(res::Matrix{Taylor1{T}}, alpha::Taylor1{T}, auxsR) where {T<:Real}
    # Local variables
    one_alpha = auxsR[1] # one(alpha)
    zero_alpha = auxsR[2] # zero(alpha)
    sin_alpha = auxsR[3] # sin(alpha)
    cos_alpha = auxsR[4] # cos(alpha)
    for ord in eachindex(alpha)
        TaylorSeries.sincos!(sin_alpha, cos_alpha, alpha, ord)
        # Matrix elements
        res[1, 1][ord] = cos_alpha[ord]
        res[2, 1][ord] = zero_alpha[ord]
        res[3, 1][ord] = sin_alpha[ord]
        res[1, 2][ord] = zero_alpha[ord]
        res[2, 2][ord] = one_alpha[ord]
        res[3, 2][ord] = zero_alpha[ord]
        res[1, 3][ord] = -sin_alpha[ord]
        res[2, 3][ord] = zero_alpha[ord]
        res[3, 3][ord] = cos_alpha[ord]
    end
    return nothing
end

@doc raw"""
    Rz(alpha::T) where {T<:Number}

Return the rotation matrix around the z-axis
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
    # Allocate memory
    res = Array{T}(undef, 3, 3)
    # Local variables
    one_alpha = one(alpha)
    zero_alpha = zero(alpha)
    sin_alpha, cos_alpha = sincos(alpha)
    # Matrix elements
    res[1, 1] = cos_alpha
    res[2, 1] = -sin_alpha
    res[3, 1] = zero_alpha
    res[1, 2] = sin_alpha
    res[2, 2] = cos_alpha
    res[3, 2] = zero_alpha
    res[1, 3] = zero_alpha
    res[2, 3] = zero_alpha
    res[3, 3] = one_alpha
    return res
end

@doc raw"""
    Rz!(res, alpha, auxsR)

Inplace version of ``Rz(alpha)``; see [`Rz`](@ref) for details.
"""
function Rz!(res::Matrix{Taylor1{T}}, alpha::Taylor1{T}, auxsR) where {T<:Real}
    # Local variables
    one_alpha = auxsR[1] # one(alpha)
    zero_alpha = auxsR[2] # zero(alpha)
    sin_alpha = auxsR[3] # sin(alpha)
    cos_alpha = auxsR[4] # cos(alpha)
    for ord in eachindex(alpha)
        TaylorSeries.sincos!(sin_alpha, cos_alpha, alpha, ord)
        # Matrix elements
        res[1, 1][ord] = cos_alpha[ord]
        res[2, 1][ord] = -sin_alpha[ord]
        res[3, 1][ord] = zero_alpha[ord]
        res[1, 2][ord] = sin_alpha[ord]
        res[2, 2][ord] = cos_alpha[ord]
        res[3, 2][ord] = zero_alpha[ord]
        res[1, 3][ord] = zero_alpha[ord]
        res[2, 3][ord] = zero_alpha[ord]
        res[3, 3][ord] = one_alpha[ord]
    end
    return nothing
end

# myatan(x, y) = y>=zero(x)?( x>=zero(x)?atan(y/x):(atan(y/x)+pi) ):( x>=zero(x)?(atan(y/x)+2pi):(atan(y/x)+pi) )
# myatan2(x, y) = y>=zero(x)?( x>=zero(x)?atan(y/x):(atan(y/x)-pi) ):( x>=zero(x)?(atan(y/x)):(atan(y/x)+pi) )

@doc raw"""
    pole_ra(t)

Return the right ascension (in radians) of the pole unit vector in the inertial frame ``\mathbf{p}_\mathrm{E}(t)``
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

Return the declination (in radians) of the pole unit vector in the inertial frame ``\mathbf{p}_\mathrm{E}(t)``
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

Return the right ascension and declination (both in radians) of the pole unit vector in the inertial frame ``\mathbf{p}_\mathrm{E}(t)``
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

Return the rotation matrix from the inertial frame to the frame with pole at right
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

Return the rotation matrix from inertial frame to Earth pole at time t (days) since J2000.0.

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

Return the IAU 1980 nutation matrix (Explanatory Supplement to the Astronomical Almanac 1992)
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

function nutation_iau80!(res, t::Taylor1, auxNut, auxsR, auxMat)
    ϵ0 = auxNut[1]        # Mean obliquity (rad)
    Δϵ = auxNut[2]        # Nutation in obliquity (rad)
    Δψ = auxNut[3]        # Nutation in longitude (rad)
    nΔψ = auxNut[4]       # Negative nutation in longitude (rad)
    nϵ = auxNut[5]        # ϵ = -(ϵ0 + Δϵ)
    # omega_t_cy = auxNut[6] # Ω(t)
    ϵ̄!(ϵ0, t, view(auxNut, 7:9))
    Delta_epsilon!(Δϵ, t, view(auxNut, 6:11))
    Delta_psi!(Δψ, t, view(auxNut, 6:11))
    for ind in eachindex(t)
        nΔψ[ind] = -Δψ[ind]
        nϵ[ind] = -(ϵ0[ind] + Δϵ[ind]) # -ϵ
    end
    Rxnϵ = auxMat[1]
    RznΔψ = auxMat[2]
    Rxϵ0 = auxMat[3]
    Rx!(Rxnϵ, nϵ, auxsR)
    Rz!(RznΔψ, nΔψ, auxsR)
    Rx!(Rxϵ0, ϵ0, auxsR)
    res = Rxnϵ*RznΔψ*Rxϵ0
    return nothing
end

@doc raw"""
    t2c_jpl_de430(t)

Return the matrix for terrestrial-to-celestial coordinate transformation, according to JPL DE 430/431 Earth orientation model
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

# The following methods are used by NEOs.jl for propagations involving
# extended body interactions
function t2c_jpl_de430(dsj2k::Taylor1{T}, zero_q::Taylor1{T}) where {T <: Real}
    M = t2c_jpl_de430(dsj2k)
    return M .+ zero_q
end

function t2c_jpl_de430(dsj2k::Taylor1{T}, zero_q::Taylor1{TaylorN{T}}) where {T <: Real}
    M = t2c_jpl_de430(dsj2k)
    return M .+ zero_q
end

function t2c_jpl_de430(dsj2k::Taylor1{T}, zero_q::Taylor1{Taylor1{T}}) where {T <: Real}
    M = t2c_jpl_de430(dsj2k)
    one_q = one(zero_q.coeffs[1])
    _M_ = @. Taylor1(getfield(M, :coeffs) * one_q)
    return _M_
end

t2c_jpl_de430!(dsj2k::Taylor1{T}, zero_q::Taylor1{<:Number}, ::Nothing) where {T <: Real} =
    t2c_jpl_de430(dsj2k, zero_q)

function t2c_jpl_de430!(dsj2k::Taylor1{T}, zero_q::Taylor1{T},
                        rotatBuf::RetAlloc{Taylor1{T}}) where {T <: Real}
    c2t_jpl_de430!(dsj2k, rotatBuf)
    resdagg = rotatBuf.v2[13]
    resdagg = transpose(rotatBuf.v2[12])
    return resdagg .+ zero_q
end

function t2c_jpl_de430!(dsj2k::Taylor1{T}, zero_q::Taylor1{TaylorN{T}},
                        rotatBuf::RetAlloc{Taylor1{T}}) where {T <: Real}
    c2t_jpl_de430!(dsj2k, rotatBuf)
    resdagg = rotatBuf.v2[13]
    resdagg = transpose(rotatBuf.v2[12])
    return resdagg .+ zero_q
end

function t2c_jpl_de430!(dsj2k::Taylor1{T}, zero_q::Taylor1{Taylor1{T}},
                        rotatBuf::RetAlloc{Taylor1{T}}) where {T <: Real}
    c2t_jpl_de430!(dsj2k, rotatBuf)
    resdagg = rotatBuf.v2[13]
    resdagg = transpose(rotatBuf.v2[12])
    one_q = one(zero_q.coeffs[1])
    return @. Taylor1(getfield(resdagg, :coeffs) * one_q)
end


@doc raw"""
    c2t_jpl_de430(t)

Return the matrix for celestial-to-terrestrial coordinate transformation, according to
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
    c2t_jpl_de430!(t, rotatBuf)

Inplace version of ``c2t_jpl_de430(t)``; see [`c2t_jpl_de430`](@ref) for details.
"""
function c2t_jpl_de430!(t::Taylor1{T}, rotatBuf::RetAlloc{Taylor1{T}}) where {T<:Real}
    # Preliminaries
    zt = rotatBuf.v0[1]
    nzt = rotatBuf.v0[2]
    Tt = rotatBuf.v0[3]
    Zt = rotatBuf.v0[4]
    nZt = rotatBuf.v0[5]
    phiy = rotatBuf.v0[6]
    phix = rotatBuf.v0[7]
    auxs_angles = rotatBuf.v1[1]
    auxsR = rotatBuf.v1[2]
    auxNut = rotatBuf.v1[3]
    Rzz = rotatBuf.v2[1]
    RyT = rotatBuf.v2[2]
    RzZ = rotatBuf.v2[3]
    P_iau7680 = rotatBuf.v2[4]
    Ryphiy = rotatBuf.v2[5]
    Rxphix = rotatBuf.v2[6]
    corrections = rotatBuf.v2[7]
    N_iau80 = rotatBuf.v2[8]
    # Rxnϵ = rotatBuf.v2[9]
    # RznΔψ = rotatBuf.v2[10]
    # Rxϵ0 = rotatBuf.v2[11]
    res = rotatBuf.v2[12]
    Raux = rotatBuf.v2[13]
    # Angles
    zeta!(zt, t, auxs_angles)  # zeta(t)
    Theta!(Tt, t, auxs_angles) # Theta(t)
    Zeta!(Zt, t, auxs_angles)
    for ord in eachindex(t)
        nzt[ord] = -zt[ord]
        nZt[ord] = -Zt[ord]
    end
    # Precession matrix, see equation (5-147) in page (5-59) of https://doi.org/10.1002/0471728470
    # P_iau7680 = Rz(-zeta(t))*Ry(Theta(t))*Rz(-Zeta(t))
    Rz!(Rzz, nzt, auxsR)
    Ry!(RyT, Tt, auxsR)
    Rz!(RzZ, nZt, auxsR)
    # Raux = Rzz*RyT
    # P_iau7680 = Raux*RzZ
    mul!(Raux, Rzz, RyT)
    mul!(P_iau7680, Raux, RzZ)
    # Linear corrections to precession, see equation (25) in page 11 and Table 10 in page 50 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract
    # corrections = Ry(phi_y(t))*Rx(phi_x(t))
    phi_y!(phiy, t)
    phi_x!(phix, t)
    Ry!(Ryphiy, phiy, auxsR)
    Rx!(Rxphix, phix, auxsR)
    # corrections = Ryphiy*Rxphix
    mul!(corrections, Ryphiy, Rxphix)
    # Nutation matrix, see equation (5-152) in page 5-60 of https://doi.org/10.1002/0471728470
    # N_iau80 = nutation_iau80(t)
    nutation_iau80!(N_iau80, t, auxNut, auxsR, view(rotatBuf.v2, 9:11))
    # Raux = N_iau80*corrections
    # res = Raux*P_iau7680
    mul!(Raux,  N_iau80, corrections)
    mul!(res,  Raux, P_iau7680)
    return nothing
end

@doc raw"""
    allocate_c2t_jpl_de430(t)

Returns ``c2t_jpl_de430(t)`` (see details in [`c2t_jpl_de430`](@ref)) and a `RetAlloc` object with the
objects that can be used to recicle memory.
"""
function allocate_c2t_jpl_de430(t::Taylor1{T}) where {T<:Real}
    t_cy = t/(100yr)
    # Angles zeta(t), Theta(t), Zeta(t), negatives, and allocations
    zt = zero(t) # zt = zeta(t)
    auxhorner = zero(t)
    angle_A = zero(t)
    auxs_angles = [t_cy, auxhorner, angle_A] # auxs for angles
    zeta!(zt, t, auxs_angles) # zt = zeta(t)
    nzt = -zt
    Tt = zero(t) # Tt = Theta(t)
    Theta!(Tt, t, auxs_angles) # Theta(t)
    Zt = zero(t) # Zt = Zeta(t)
    Zeta!(Zt, t, auxs_angles)
    nZt = -Zt
    # Allocations for calculations involving Rx, Ry and Rz
    one_alpha = one(nzt)
    zero_alpha = zero(nzt)
    sin_alpha, cos_alpha = sincos(nzt)
    auxsR = [one_alpha, zero_alpha, sin_alpha, cos_alpha]
    Rzz = Rz(nzt)
    RyT = Ry(Tt)
    RzZ = Rz(nZt)
    # Precession matrix, see equation (5-147) in page (5-59) of https://doi.org/10.1002/0471728470
    # P_iau7680 = Rz(-zeta(t))*Ry(Theta(t))*Rz(-Zeta(t))
    Raux = Rzz*RyT
    P_iau7680 = Raux*RzZ
    # Angles and allocations for linear corrections of precession
    phiy = phi_y(t)
    phix = phi_x(t)
    Ryphiy = Ry(phiy)
    Rxphix = Rx(phix)
    # Linear corrections to precession, see equation (25) in page 11 and Table 10 in page 50 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract
    # corrections = Ry(phi_y(t))*Rx(phi_x(t))
    corrections = Ryphiy*Rxphix
    # Nutation matrix, see equation (5-152) in page 5-60 of https://doi.org/10.1002/0471728470
    # N_iau80 = nutation_iau80(t)
    ϵ0 = zero(t)        # Mean obliquity (rad)
    Δϵ = zero(t)        # Nutation in obliquity (rad)
    Δψ = zero(t)        # Nutation in longitude (rad)
    nΔψ = zero(t)       # Negative nutation in longitude (rad)
    nϵ = zero(t)        # ϵ = -(ϵ0 + Δϵ)
    omega_t_cy = zero(t) # Ω(t)
    auxNut = [ϵ0, Δϵ, Δψ, nΔψ, nϵ, omega_t_cy, t_cy, auxhorner, angle_A, sin_alpha, cos_alpha]
    ϵ̄!(ϵ0, t, auxs_angles)
    Delta_epsilon!(Δϵ, t, view(auxNut, 6:11))
    Delta_psi!(Δψ, t, view(auxNut, 6:11))
    nΔψ = -Δψ           # Negative of nutation in longitude (rad)
    nϵ = -(ϵ0 + Δϵ)
    Rxnϵ = Rx(nϵ)
    RznΔψ = Rz(nΔψ)
    Rxϵ0 = Rx(ϵ0)
    Raux = Rxnϵ*RznΔψ
    N_iau80 = Raux*Rxϵ0
    # c2t_jpl_de430
    Raux = N_iau80*corrections
    res = Raux*P_iau7680
    # resdagg = Matrix(transpose(res))
    # Returned RetAlloc{Taylor1{T}} object
    rotatBuf = RetAlloc{Taylor1{T}}(
        [zt, nzt, Tt, Zt, nZt, phiy, phix],
        [auxs_angles, auxsR, auxNut],
        [Rzz, RyT, RzZ, P_iau7680, Ryphiy, Rxphix, corrections, N_iau80, Rxnϵ, RznΔψ, Rxϵ0, res, Raux],
        [Array{Taylor1{T}}(undef, 0, 0, 0)],
        [Array{Taylor1{T}}(undef, 0, 0, 0, 0)])
    return rotatBuf
end

@doc raw"""
    moon_omega(ϕ::Taylor1, θ::Taylor1, ψ::Taylor1)

Return the Moon's angular velocity, computed by differentiating the Euler angles
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
    sinθ, cosθ = sincos(θ)
    sinψ, cosψ = sincos(ψ)
    dϕ = ordpres_differentiate(ϕ)
    dθ = ordpres_differentiate(θ)
    dψ = ordpres_differentiate(ψ)
    ωx = dϕ*sinθ*sinψ + dθ*cosψ
    ωy = dϕ*sinθ*cosψ - dθ*sinψ
    ωz = dψ + dϕ*cosθ
    return [ωx, ωy, ωz]
end

@doc raw"""
    ITM1(x::T, y::T, z::T) where {T <: Number}

Return the first term of the time-dependent part of lunar total moment of intertia
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

Return the second term of the time-dependent part of lunar total moment of intertia
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

Return lunar mantle inertia tensor
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