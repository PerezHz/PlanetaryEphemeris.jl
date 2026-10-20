# This file is part of the PlanetaryEphemeris.jl package; MIT licensed

# Orientation models common to the JPL DE430 and DE440 ephemerides: rotation matrices,
# IAU 1980 nutation restricted to the 18.6-year terms, mean obliquity, and the
# orientation, angular velocity and moment of inertia of the Moon

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
Ω(t) = deg2rad( evaluate(Taylor1(coeffs_Ω_A, TS.order(t)), t/36525) )

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
    Delta_psi(t, Δψ_A)

Return the nutation in longitude (in radians)
```math
\Delta\psi(t) = \Delta\psi_A \sin\Omega,
```
where ``t = 36,525 T`` is the TDB time in Julian days from J2000.0, ``\Delta\psi_A`` is
the amplitude of the 18.6-year term (in arcseconds; ``-17''.1996`` in the IAU 1980 model)
and ``\Omega`` is the longitude of the mean ascending node of the lunar orbit.

See equation (5-190) in page 5-72 of https://doi.org/10.1002/0471728470.

See also [`Ω`](@ref).
"""
Delta_psi(t, Δψ_A) = deg2rad( (Δψ_A/3600)*sin(Ω(t)) )

@doc raw"""
    Delta_psi!(res, t, auxs, Δψ_A)

Inplace version of ``Delta_psi(t, Δψ_A)``; see [`Delta_psi`](@ref) for details.
"""
function Delta_psi!(deltapsi::Taylor1{T}, t::Taylor1{T}, auxs, Δψ_A) where {T<:Real}
    numfactor = Δψ_A/3600
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
    Delta_epsilon(t, Δϵ_A)

Return the nutation in obliquity (in radians)
```math
\Delta\epsilon(t) = \Delta\epsilon_A \cos\Omega,
```
where ``t = 36,525 T`` is the TDB time in Julian days from J2000.0, ``\Delta\epsilon_A`` is
the amplitude of the 18.6-year term (in arcseconds; ``9''.2025`` in the IAU 1980 model)
and ``\Omega`` is the longitude of the mean ascending node of the lunar orbit.

See equation (5-190) in page 5-72 of https://doi.org/10.1002/0471728470.

See also [`Ω`](@ref).
"""
Delta_epsilon(t, Δϵ_A) = deg2rad( (Δϵ_A/3600)*cos(Ω(t)) )

@doc raw"""
    Delta_epsilon!(res, t, auxs, Δϵ_A)

Inplace version of ``Delta_epsilon(t, Δϵ_A)``; see [`Delta_epsilon`](@ref) for details.
"""
function Delta_epsilon!(deltaepsilon::Taylor1{T}, t::Taylor1{T}, auxs, Δϵ_A) where {T<:Real}
    numfactor = Δϵ_A/3600
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
    pole_date(t, Δψ_A, Δϵ_A)

Return the true pole of date unit vector ``\mathbf{p}_\mathrm{d}``, computed by rotating the Earth-fixed pole vector by the effect of the 18.6-year nutation term with amplitudes ``\Delta\psi_A`` and ``\Delta\epsilon_A`` (in arcseconds). ``t`` is the TDB time in Julian days from J2000.0.

See equation (23) in page 11 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract.

See also [`nutation_iau80`](@ref).
"""
function pole_date(t, Δψ_A, Δϵ_A)
    # ϵ0 = ϵ̄(t)
    # Δϵ = Delta_epsilon(t)
    # Δψ = Delta_psi(t)
    # ϵ = ϵ0 + Δϵ
    # pdx = sin(Δψ)*sin(ϵ)
    # pdy = cos(Δψ)*sin(ϵ)*cos(ϵ0) - cos(ϵ)*sin(ϵ0)
    # pdz = cos(Δψ)*sin(ϵ)*sin(ϵ0) + cos(ϵ)*cos(ϵ0)
    Nut_mat = nutation_iau80(t, Δψ_A, Δϵ_A)
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
ϵ̄(t) = deg2rad( evaluate(Taylor1(coeffs_ϵ̄, TS.order(t)), t/36525) )

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
    nutation_iau80(t, Δψ_A, Δϵ_A)

Return the IAU 1980 nutation matrix (Explanatory Supplement to the Astronomical Almanac 1992),
including only the 18.6-year terms with amplitudes ``\Delta\psi_A`` and ``\Delta\epsilon_A``
(in arcseconds)
```math
N = R_x(-\bar{\epsilon}(t) - \Delta\epsilon(t))R_z(-\Delta\psi(t))R_x(\bar{\epsilon}(t)),
```
where ``t`` is the TDB time in Julian days from J2000.0, ``\bar{\epsilon}`` is the mean
obliquity, ``\Delta \psi`` is the nutation in longitude and ``\Delta\epsilon`` is the
nutation in obliquity. All angles are in radians.

See equation (5-152) in page 5-60 of https://doi.org/10.1002/0471728470.

See also [`ϵ̄`](@ref), [`Delta_epsilon`](@ref) and [`Delta_psi`](@ref).
"""
function nutation_iau80(t, Δψ_A, Δϵ_A)
    ϵ0 = ϵ̄(t)                     # Mean obliquity (rad)
    Δϵ = Delta_epsilon(t, Δϵ_A)   # Nutation in obliquity (rad)
    Δψ = Delta_psi(t, Δψ_A)       # Nutation in longitude (rad)
    ϵ = ϵ0 + Δϵ
    return Rx(-ϵ)*Rz(-Δψ)*Rx(ϵ0)
end

@doc raw"""
    nutation_iau80!(res, t, auxNut, auxsR, auxMat, Δψ_A, Δϵ_A)

Inplace version of ``nutation_iau80(t, Δψ_A, Δϵ_A)``; see [`nutation_iau80`](@ref) for details.
"""
function nutation_iau80!(res, t::Taylor1, auxNut, auxsR, auxMat, Δψ_A, Δϵ_A)
    ϵ0 = auxNut[1]        # Mean obliquity (rad)
    Δϵ = auxNut[2]        # Nutation in obliquity (rad)
    Δψ = auxNut[3]        # Nutation in longitude (rad)
    nΔψ = auxNut[4]       # Negative nutation in longitude (rad)
    nϵ = auxNut[5]        # ϵ = -(ϵ0 + Δϵ)
    # omega_t_cy = auxNut[6] # Ω(t)
    ϵ̄!(ϵ0, t, view(auxNut, 7:9))
    Delta_epsilon!(Δϵ, t, view(auxNut, 6:11), Δϵ_A)
    Delta_psi!(Δψ, t, view(auxNut, 6:11), Δψ_A)
    for ind in eachindex(t)
        nΔψ[ind] = -Δψ[ind]
        nϵ[ind] = -(ϵ0[ind] + Δϵ[ind]) # -ϵ
    end
    Rxnϵ = auxMat[1]
    RznΔψ = auxMat[2]
    Rxϵ0 = auxMat[3]
    Raux =  auxMat[5]
    Rx!(Rxnϵ, nϵ, auxsR)
    Rz!(RznΔψ, nΔψ, auxsR)
    Rx!(Rxϵ0, ϵ0, auxsR)
    matmul!(Raux, RznΔψ, Rxϵ0)
    matmul!(res, Rxnϵ, Raux)
    return nothing
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
    ITM1(x::T, y::T, z::T, consts) where {T <: Number}

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

The physical parameters of the model are taken from `consts` (see [`DE430Constants`](@ref)).

See the first term of equation (41) in page 17 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract.

See also [`ITM2`](@ref) and [`ITM`](@ref).
"""
function ITM1(x::T, y::T, z::T, consts) where {T <: Number}
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
    return (-consts.ITM1_coeff/r5)*m
end

@doc raw"""
    ITM2(ωx::T, ωy::T, ωz::T, consts) where {T <: Number}

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

The physical parameters of the model are taken from `consts` (see [`DE430Constants`](@ref)).

See the second term of equation (41) in page 17 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract.

See also [`ITM1`](@ref) and [`ITM`](@ref).
"""
function ITM2(ωx::T, ωy::T, ωz::T, consts) where {T <: Number}
    # Compute corresponding matrix elements
    m = Matrix{T}(undef, 3, 3)
    ωx2 = ωx*ωx
    ωy2 = ωy*ωy
    ωz2 = ωz*ωz
    # Angular speed
    ω2 = ωx2+ωy2+ωz2
    aux = (ω2-consts.n_moon_p2)/3
    # Diagonal
    m[1,1] = ωx2-aux
    m[2,2] = ωy2-aux
    m[3,3] = ωz2-(ω2+2*consts.n_moon_p2)/3
    # Off diagonal (the matrix is symmetric)
    m[1,2] = ωx*ωy
    m[2,1] = m[1,2]
    m[1,3] = ωx*ωz
    m[3,1] = m[1,3]
    m[2,3] = ωy*ωz
    m[3,2] = m[2,3]
    return consts.ITM2_coeff*m
end

@doc raw"""
    ITM(q::Vector{T}, eulang::Vector{T}, ω_m::Vector{T}, consts) where {T <: Number}

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

The physical parameters of the model are taken from `consts` (see [`DE430Constants`](@ref)).

See equation (41) in page 17 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract.

See also [`ITM1`](@ref) and [`ITM2`](@ref).
"""
function ITM(q::Vector{T}, eulang::Vector{T}, ω_m::Vector{T}, consts) where {T <: Number}
    # Coordinate transformation matrix: inertial frame -> lunar mantle frame
    M_del_mo = pole_rotation(eulang[1] - (pi/2), (pi/2) - eulang[2], eulang[3])
    # Transform delayed geocentric position of Moon (space-fixed->lunar mantle frame)
    X_me_del_τ_M = q[3mo-2:3mo] .- q[3ea-2:3ea]
    X_me_del_τ_M_lm = M_del_mo*X_me_del_τ_M
    # First term of the time-dependent part of lunar total moment of intertia
    itm1 = ITM1(X_me_del_τ_M_lm[1], X_me_del_τ_M_lm[2], X_me_del_τ_M_lm[3], consts)
    # Second term of the time-dependent part of lunar total moment of intertia
    itm2 = ITM2(ω_m[1], ω_m[2], ω_m[3], consts)
    # Lunar mantle intertia tensor
    return consts.I_m_und*one(q[1]) + itm1 + itm2
end
