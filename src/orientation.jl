# This file is part of the PlanetaryEphemeris.jl package; MIT licensed

# Orientation models: functions shared by the JPL DE430/431 and DE440/441 models

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
A = R_z(W)R_x\left(\frac{\pi}{2} - \delta\right)R_z\left(\alpha + \frac{\pi}{2}\right).
```

See equation (6-3) in page 6-5 of https://doi.org/10.1002/0471728470.

``\textbf{Note:}`` Rotation matrices ``(R_x, R_y, R_z)`` convention is the same as
equations (11), (12) and (13) in page 9 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract.

See also [`Rx`](@ref) and [`Rz`](@ref).
"""
function pole_rotation(α::T, δ::T, W::T=zero(α)) where {T <: Number}
    return Rz(W)*Rx(π/2-δ)*Rz(π/2+α)
end

# IAU 1980 nutation model including only the 18.6-year term, and mean obliquity of the
# ecliptic. Both are used by the DE430/431 and DE440/441 Earth orientation models; see
# equations (21)-(24) of Folkner et al. (2014) and (20)-(24) of Park et al. (2021).

# Folkner et al. (2014), eq. (21), prints the rate as 1934°08'10''.549: a typo confirmed by
# JPL (see the README)
const coeffs_Ω_A = [125+2/60+40.280/3600, -(1934+8/60+10.539/3600), 7.455/3600, 0.008/3600]

@doc raw"""
    Ω(t::Number)

Return the longitude (in radians) of the mean ascending node of the lunar orbit on the
ecliptic, measured from the mean equinox of date
```math
\Omega(t) = 125^\circ 02' 40''.280 - \left(1934^\circ 8' 10''.539\right) T + 7''.455 T^2 + 0''.008 T^3,
```
where ``t = 36,525 T`` is the TDB time in Julian days from J2000.0.

See equation (5-64) in page 5-27 of https://doi.org/10.1002/0471728470.
"""
Ω(t::Real) = deg2rad(evalpoly(t/36525, coeffs_Ω_A))

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

# Folkner et al. (2014), eq. (22), prints -17''.206262: a typo confirmed by JPL (see the
# README)
const coeff_Δψ = -17.1996 # arcseconds

@doc raw"""
    Delta_psi(t::Number)

Return the nutation in longitude (in radians)
```math
\Delta\psi(t) = -17''.1996 \sin\Omega,
```
where ``t = 36,525 T`` is the TDB time in Julian days from J2000.0 and ``\Omega`` is the
longitude of the mean ascending node of the lunar orbit.

See equation (5-190) in page 5-72 of https://doi.org/10.1002/0471728470.

See also [`Ω`](@ref).
"""
Delta_psi(t::Real) = deg2rad((coeff_Δψ/3600)*sin(Ω(t)))

@doc raw"""
    Delta_psi!(res, t, auxs)

Inplace version of ``Delta_psi(t)``; see [`Delta_psi`](@ref) for details.
"""
function Delta_psi!(deltapsi::Taylor1{T}, t::Taylor1{T}, auxs) where {T<:Real}
    numfactor = coeff_Δψ/3600
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

# Folkner et al. (2014), eq. (22), prints 9''.205348: a typo confirmed by JPL (see the
# README)
const coeff_Δϵ = 9.2025 # arcseconds

@doc raw"""
    Delta_epsilon(t::Number)

Return the nutation in obliquity (in radians)
```math
\Delta\epsilon(t) = 9''.2025 \cos\Omega,
```
where ``t = 36,525 T`` is the TDB time in Julian days from J2000.0 and ``\Omega`` is the
longitude of the mean ascending node of the lunar orbit.

See equation (5-190) in page 5-72 of https://doi.org/10.1002/0471728470.

See also [`Ω`](@ref).
"""
Delta_epsilon(t::Real) = deg2rad((coeff_Δϵ/3600)*cos(Ω(t)))

@doc raw"""
    Delta_epsilon!(res, t, auxs)

Inplace version of ``Delta_epsilon(t)``; see [`Delta_epsilon`](@ref) for details.
"""
function Delta_epsilon!(deltaepsilon::Taylor1{T}, t::Taylor1{T}, auxs) where {T<:Real}
    numfactor = coeff_Δϵ/3600
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

const coeffs_ϵ̄ = [84381.448/3600, -46.815/3600, -0.00059/3600, 0.001813/3600]

@doc raw"""
    ϵ̄(t::Number)

Return the mean obliquity (in radians)
```math
\bar{\epsilon}(t) = 84,381''.448 - 46''.8150 T - 0''.00059 T^2 + 0''.001813 T^3,
```
where ``t = 36,525 T`` is the TDB time in Julian days from J2000.0.

See equation (5-153) in page 5-61 of https://doi.org/10.1002/0471728470.
"""
ϵ̄(t::Real) = deg2rad(evalpoly(t/36525, coeffs_ϵ̄))

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
    nutation_iau80(t::Number)

Return the nutation matrix of the JPL DE430/431 Earth orientation model, i.e. the IAU 1980
nutation model with only the 18.6-year term
```math
N = R_x(-\bar{\epsilon}(t) - \Delta\epsilon(t))R_z(-\Delta\psi(t))R_x(\bar{\epsilon}(t)),
```
where ``t`` is the TDB time in Julian days from J2000.0, ``\bar{\epsilon}`` is the mean
obliquity, ``\Delta \psi`` is the nutation in longitude and ``\Delta\epsilon`` is the
nutation in obliquity. All angles are in radians.

See equation (5-152) in page 5-60 of https://doi.org/10.1002/0471728470.

See also [`ϵ̄`](@ref), [`Delta_epsilon`](@ref) and [`Delta_psi`](@ref).
"""
function nutation_iau80(t::Real)
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
    Raux =  auxMat[5]
    Rx!(Rxnϵ, nϵ, auxsR)
    Rz!(RznΔψ, nΔψ, auxsR)
    Rx!(Rxϵ0, ϵ0, auxsR)
    matmul!(Raux, RznΔψ, Rxϵ0)
    matmul!(res, Rxnϵ, Raux)
    return nothing
end

# Allocating Taylor1 methods: fall back to the in-place methods

# Auxiliaries for Ω!, ϵ̄!, Zeta!, Theta! and zeta!
_allocate_angle_auxs(t::Taylor1{T}) where {T<:Real} = [zero(t), zero(t), Taylor1(zeros(T, 4))]

# Auxiliaries for Delta_psi! and Delta_epsilon!
_allocate_nutation_auxs(t::Taylor1{T}) where {T<:Real} =
    [zero(t); _allocate_angle_auxs(t); zero(t); zero(t)]

for (f, f!, auxs) in ((:Ω, :Ω!, :_allocate_angle_auxs), (:ϵ̄, :ϵ̄!, :_allocate_angle_auxs),
        (:Delta_psi, :Delta_psi!, :_allocate_nutation_auxs),
        (:Delta_epsilon, :Delta_epsilon!, :_allocate_nutation_auxs))
    @eval function $f(t::Taylor1{T}) where {T<:Real}
        res = zero(t)
        $(f!)(res, t, $(auxs)(t))
        return res
    end
end

# Auxiliaries for nutation_iau80! (same layout as in allocate_c2t_jpl_de430)
function _allocate_nutation_iau80(t::Taylor1{T}) where {T<:Real}
    sin_alpha, cos_alpha = zero(t), zero(t)
    auxNut = [zero(t), zero(t), zero(t), zero(t), zero(t), zero(t), zero(t), zero(t),
        Taylor1(zeros(T, 4)), sin_alpha, cos_alpha]
    auxsR = [one(t), zero(t), sin_alpha, cos_alpha]
    auxMat = [[zero(t) for _ in 1:3, _ in 1:3] for _ in 1:5]
    return auxNut, auxsR, auxMat
end

function nutation_iau80(t::Taylor1{T}) where {T<:Real}
    auxNut, auxsR, auxMat = _allocate_nutation_iau80(t)
    N_iau80 = [zero(t) for _ in 1:3, _ in 1:3]
    nutation_iau80!(N_iau80, t, auxNut, auxsR, auxMat)
    return N_iau80
end

# Lunar Euler angles

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
