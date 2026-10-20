# This file is part of the PlanetaryEphemeris.jl package; MIT licensed

# JPL-DE-430/431 model for the orientation of the Earth

# From JPL DE430 documentation (Folkner et al., 2014):
# Only the long-term change of the Earth's orientation is modeled in the ephemeris
# integration. The Earth orientation model used for the DE 430 and 431 integration
# is based on the International Astronomical Union (IAU) 1976 precession model\* with
# an estimated linear correction and on a modified IAU 1980 nutation model \* including
# only terms with a period of 18.6 years.

@doc raw"""
    pole_frame(t, consts::DE430Constants)

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
function pole_frame(t, consts::DE430Constants)
    # True pole of date unit vector
    p_d = pole_date(t, consts.nut_Δψ, consts.nut_Δϵ)
    # Inverse of precession matrix with linear corrections to precession angles
    pm = (Rz(Zeta(t, consts))*Ry(-Theta(t, consts))*Rz(zeta(t, consts))*Rx(-phi_x(t, consts))*
        Ry(-phi_y(t, consts)))
    return pm*p_d
end

# where

@doc raw"""
    phi_x(t, consts::DE430Constants)

Return the X-axis linear correction to precession (in radians)
```math
\phi_x = \phi_{x0} + 100T\times \frac{d\phi_x}{dt},
```
where ``t = 36,525 T`` is the TDB time in Julian days from J2000.0, ``\phi_{x0}`` is the
X–axis rotation at J2000.0 (in arcseconds, `ROTEX`) and ``\frac{d\phi_x}{dt}`` is the negative
obliquity rate correction (in arcseconds per year, `DROTEX`).

See equation (25) in page 11 and Table 10 in page 50 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract.

See also [`pole_frame`](@ref).
"""
function phi_x(t, consts::DE430Constants)
    phix = zero(t)
    for ord in eachindex(t)
        phix[ord] = (t[ord]/yr)*consts.DROTEX/3600
    end
    phix[0] += consts.ROTEX/3600
    for ord in eachindex(t)
        phix[ord] = deg2rad( phix[ord] )
    end
    return phix
end

@doc raw"""
    phi_x!(res, t, consts::DE430Constants)

Inplace version of ``phi_x(t, consts)``; see [`phi_x`](@ref) for details.
"""
function phi_x!(phix::Taylor1{T}, t::Taylor1{T}, consts::DE430Constants) where {T<:Real}
    for ord in eachindex(t)
        phix[ord] = (t[ord]/yr)*consts.DROTEX/3600
    end
    phix[0] += consts.ROTEX/3600
    for ord in eachindex(t)
        phix[ord] = phix[ord]*(convert(float(T), pi) / 180)#deg2rad( phix[ord] )
    end
    return nothing
end

@doc raw"""
    phi_y(t, consts::DE430Constants)

Return the Y-axis linear correction to precession (in radians)
```math
\phi_y = \phi_{y0} + 100T\times \frac{d\phi_y}{dt},
```
where ``t = 36,525 T`` is the TDB time in Julian days from J2000.0, ``\phi_{y0}`` is the
Y–axis rotation at J2000.0 (in arcseconds, `ROTEY`) and ``\frac{d\phi_y}{dt}`` is precession rate
correction times sine of obliquity (in arcseconds per year, `DROTEY`).

See equation (25) in page 11 and Table 10 in page 50 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract.

See also [`pole_frame`](@ref).
"""
function phi_y(t, consts::DE430Constants)
    phiy = zero(t)
    for ord in eachindex(t)
        phiy[ord] = (t[ord]/yr)*consts.DROTEY/3600
    end
    phiy[0] += consts.ROTEY/3600
    for ord in eachindex(t)
        phiy[ord] = deg2rad( phiy[ord] )
    end
    return phiy
end

@doc raw"""
    phi_y!(res, t, consts::DE430Constants)

Inplace version of ``phi_y(t, consts)``; see [`phi_y`](@ref) for details.
"""
function phi_y!(phiy::Taylor1{T}, t::Taylor1{T}, consts::DE430Constants) where {T<:Real}
    for ord in eachindex(t)
        phiy[ord] = (t[ord]/yr)*consts.DROTEY/3600
    end
    phiy[0] += consts.ROTEY/3600
    for ord in eachindex(t)
        phiy[ord] = phiy[ord]*(convert(float(T), pi) / 180)#deg2rad( phiy[ord] )
    end
    return nothing
end

@doc raw"""
    Zeta(t, consts::DE430Constants)

Return the ``\zeta_A`` equatorial precession angle (in radians)
```math
\zeta_A(t) = 2306''.2181 T + 0''.30188 T^2 + 0''.017998 T^3,
```
where ``t = 36,525 T`` is the TDB time in Julian days from J2000.0.

See equation (5-143) in page 5-58 of https://doi.org/10.1002/0471728470.

See also [`Theta`](@ref) and [`zeta`](@ref).
"""
function Zeta(t, consts::DE430Constants)
    t_cy = t/(100yr)
    # Zeta_arcsec = 0.017998*t_cy
    # Zeta_arcsec = (0.30188 + Zeta_arcsec)*t_cy
    # Zeta_arcsec = (2306.2181 + Zeta_arcsec)*t_cy
    ZetaA = Taylor1(consts.coeffs_Zeta, TS.order(t))
    Zeta_arcsec = ZetaA(t_cy)
    return deg2rad(Zeta_arcsec/3600)
end

@doc raw"""
    Zeta!(res, t, auxs, consts::DE430Constants)

Inplace version of ``Zeta(t, consts)``; see [`Zeta`](@ref) for details.
"""
function Zeta!(res::Taylor1{T}, t::Taylor1{T}, auxs, consts::DE430Constants) where {T<:Real}
    t_cy = auxs[1]
    auxhorner = auxs[2] # aux for _horner!
    Zeta_A = auxs[3] # Zeta_A
    for i in eachindex(consts.coeffs_Zeta)
        Zeta_A.coeffs[i] = consts.coeffs_Zeta[i]
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

@doc raw"""
    Theta(t, consts::DE430Constants)

Return the ``\theta_A`` equatorial precession angle (in radians)
```math
\theta_A(t) = 2004''.3109 T - 0''.42665 T^2 - 0''.041833 T^3,
```
where ``t = 36,525 T`` is the TDB time in Julian days from J2000.0.

See equation (5-143) in page 5-58 of https://doi.org/10.1002/0471728470.

See also [`Zeta`](@ref) and [`zeta`](@ref).
"""
function Theta(t, consts::DE430Constants)
    t_cy = t/(100yr)
    # Theta_arcsec = -0.041833*t_cy
    # Theta_arcsec = (-0.42665 + Theta_arcsec)*t_cy
    # Theta_arcsec = (2004.3109 + Theta_arcsec)*t_cy
    theta_A = Taylor1(consts.coeffs_Theta, TS.order(t))
    Theta_arcsec = theta_A(t_cy)
    return deg2rad(Theta_arcsec/3600)
end

@doc raw"""
    Theta!(res, t, auxs, consts::DE430Constants)

Inplace version of ``Theta(t, consts)``; see [`Theta`](@ref) for details.
"""
function Theta!(res::Taylor1{T}, t::Taylor1{T}, auxs, consts::DE430Constants) where {T<:Real}
    t_cy = auxs[1]
    auxhorner = auxs[2]
    theta_A = auxs[3]  # Theta_A
    for i in eachindex(consts.coeffs_Theta)
        theta_A.coeffs[i] = consts.coeffs_Theta[i]
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

@doc raw"""
    zeta(t, consts::DE430Constants)

Return the ``z_A`` equatorial precession angle (in radians)
```math
z_A(t) = 2306''.2181 T + 1''.09468 T^2 + 0''.018203 T^3,
```
where ``t = 36,525 T`` is the TDB time in Julian days from J2000.0.

See equation (5-143) in page 5-58 of https://doi.org/10.1002/0471728470.

See also [`Zeta`](@ref) and [`Theta`](@ref).
"""
function zeta(t, consts::DE430Constants)
    t_cy = t/(100yr)
    # zeta_arcsec = 0.018203*t_cy
    # zeta_arcsec = (1.09468 + zeta_arcsec)*t_cy
    # zeta_arcsec = (2306.2181 + zeta_arcsec)*t_cy
    zeta_A = Taylor1(consts.coeffs_zeta, TS.order(t))
    zeta_arcsec = zeta_A(t_cy)
    return deg2rad(zeta_arcsec/3600)
end

@doc raw"""
    zeta!(res, t, auxs, consts::DE430Constants)

Inplace version of ``zeta(t, consts)``; see [`zeta`](@ref) for details.
"""
function zeta!(res::Taylor1{T}, t::Taylor1{T}, auxs, consts::DE430Constants) where {T<:Real}
    t_cy = auxs[1]
    auxhorner = auxs[2] # aux for _horner!
    zeta_A = auxs[3] # zeta_A
    for i in eachindex(consts.coeffs_zeta)
        zeta_A.coeffs[i] = consts.coeffs_zeta[i]
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

@doc raw"""
    pole_ra(t, consts::DE430Constants)

Return the right ascension (in radians) of the pole unit vector in the inertial frame ``\mathbf{p}_\mathrm{E}(t)``
```math
\alpha = \arctan\left(\frac{p_\mathrm{Ey}(t)}{p_\mathrm{Ex}(t)}\right) + \pi,
```
where ``t`` is the TDB time in Julian days from J2000.0.

See also [`pole_frame`](@ref), [`pole_dec`](@ref) and [`pole_radec`](@ref).
"""
function pole_ra(t, consts::DE430Constants)
    pole_frame_t = pole_frame(t, consts)
    return atan( pole_frame_t[2]/pole_frame_t[1] ) + π
end

@doc raw"""
    pole_dec(t, consts::DE430Constants)

Return the declination (in radians) of the pole unit vector in the inertial frame ``\mathbf{p}_\mathrm{E}(t)``
```math
\delta = \arctan\left(\frac{p_\mathrm{Ez}(t)}{\sqrt{p_\mathrm{Ex}^2(t) + p_\mathrm{Ey}^2(t)}}\right),
```
where ``t`` is the TDB time in Julian days from J2000.0.

See also [`pole_frame`](@ref), [`pole_ra`](@ref) and [`pole_radec`](@ref).
"""
function pole_dec(t, consts::DE430Constants)
    pole_frame_t = pole_frame(t, consts)
    return atan(  pole_frame_t[3]/sqrt(pole_frame_t[1]^2+pole_frame_t[2]^2)  )
end

@doc raw"""
    pole_radec(t, consts::DE430Constants)

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
function pole_radec(t, consts::DE430Constants)
    pole_frame_t = pole_frame(t, consts)
    pole_ra = atan( pole_frame_t[2]/pole_frame_t[1] ) + π
    pole_dec = atan(  pole_frame_t[3]/sqrt(pole_frame_t[1]^2+pole_frame_t[2]^2)  )
    return pole_ra, pole_dec
end

@doc raw"""
    earth_pole_rotation(t, consts::DE430Constants)

Return the rotation matrix from inertial frame to Earth pole at time t (days) since J2000.0.

See also [`pole_radec`](@ref) and [`pole_rotation`](@ref).
"""
function earth_pole_rotation(t, consts::DE430Constants)
    # α_ep = Earth pole at time t (Julian days) since J2000.0
    # δ_ep = Earth pole at time t (Julian days) since J2000.0
    α_ep, δ_ep = pole_radec(t, consts)
    return pole_rotation(α_ep, δ_ep)
end

@doc raw"""
    t2c_jpl_de430(t, consts::DE430Constants)

Return the matrix for terrestrial-to-celestial coordinate transformation, according to JPL DE 430/431 Earth orientation model
```math
\texttt{t2c_jpl_de430}(t) = \texttt{c2t_jpl_de430}^T(t),
```
where ``t`` is the TDB time in Julian days from J2000.0.

See also [`c2t_jpl_de430`](@ref).

"""
function t2c_jpl_de430(t, consts::DE430Constants)
    # inv_P_iau7680 = Rz(Zeta(t))*Ry(-Theta(t))*Rz(zeta(t))
    # corrections = Rx(-phi_x(t))*Ry(-phi_y(t))
    # # inv_N_iau80 = inv(nutation_iau80(t))
    # inv_N_iau80 = transpose(nutation_iau80(t))
    # return inv_P_iau7680*corrections*inv_N_iau80 #corrections*inv_P_iau7680*inv_N_iau80
    return transpose(c2t_jpl_de430(t, consts))
end

# The following methods are used by NEOs.jl for propagations involving
# extended body interactions
function t2c_jpl_de430(dsj2k::Taylor1{T}, zero_q::Taylor1{T},
        consts::DE430Constants) where {T <: Real}
    M = t2c_jpl_de430(dsj2k, consts)
    return M .+ zero_q
end

function t2c_jpl_de430(dsj2k::Taylor1{T}, zero_q::Taylor1{TaylorN{T}},
        consts::DE430Constants) where {T <: Real}
    M = t2c_jpl_de430(dsj2k, consts)
    return M .+ zero_q
end

function t2c_jpl_de430(dsj2k::Taylor1{T}, zero_q::Taylor1{Taylor1{T}},
        consts::DE430Constants) where {T <: Real}
    M = t2c_jpl_de430(dsj2k, consts)
    one_q = one(zero_q.coeffs[1])
    _M_ = @. Taylor1(getfield(M, :coeffs) * one_q)
    return _M_
end

@doc raw"""
    t2c_jpl_de430!(M, ea, t, zero_q, rotatBuf::Union{Nothing,RetAlloc{Taylor1{T}}}, consts::DE430Constants)

Inplace version of ``t2c_jpl_de430(t, consts)``; see [`t2c_jpl_de430`](@ref) for details.
"""
t2c_jpl_de430!(::Array{Taylor1{S},3}, ::Int, dsj2k::Taylor1{T},
        zero_q::Taylor1{S}, ::Nothing, consts::DE430Constants) where {T <: Real, S <: Number} =
    t2c_jpl_de430(dsj2k, zero_q, consts)

function t2c_jpl_de430!(M::Array{Taylor1{T},3}, ea::Int, dsj2k::Taylor1{T},
        ::Taylor1{T}, rotatBuf::RetAlloc{Taylor1{T}}, consts::DE430Constants) where {T <: Real}
    c2t_jpl_de430!(dsj2k, rotatBuf, consts)
    for j in axes(M, 2)
        for i in axes(M, 1)
            for k in eachindex(M[i, j, ea])
                M[i, j, ea][k] = rotatBuf.v2[12][j, i][k]
            end
        end
    end
    return M
end

function t2c_jpl_de430!(M::Array{Taylor1{TaylorN{T}},3}, ea::Int, dsj2k::Taylor1{T},
        ::Taylor1{TaylorN{T}}, rotatBuf::RetAlloc{Taylor1{T}}, consts::DE430Constants) where {T <: Real}
    c2t_jpl_de430!(dsj2k, rotatBuf, consts)
    for j in axes(M, 2)
        for i in axes(M, 1)
            for ord in eachindex(M[i, j, ea])
                M[i, j, ea][ord][0][1] = rotatBuf.v2[12][j, i][ord]
            end
        end
    end
    return M
end

function t2c_jpl_de430!(M::Array{Taylor1{Taylor1{T}},3}, ea::Int, dsj2k::Taylor1{T},
        ::Taylor1{Taylor1{T}}, rotatBuf::RetAlloc{Taylor1{T}}, consts::DE430Constants) where {T <: Real}
    c2t_jpl_de430!(dsj2k, rotatBuf, consts)
    for j in axes(M, 2)
        for i in axes(M, 1)
            for ord in eachindex(M[i, j, ea])
                M[i, j, ea][ord][0] = rotatBuf.v2[12][j, i][ord]
            end
        end
    end
    return M
end


@doc raw"""
    c2t_jpl_de430(t, consts::DE430Constants)

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
function c2t_jpl_de430(t, consts::DE430Constants)
    # Precession matrix, see equation (5-147) in page (5-59) of https://doi.org/10.1002/0471728470
    P_iau7680 = Rz(-zeta(t, consts))*Ry(Theta(t, consts))*Rz(-Zeta(t, consts))
    # Linear corrections to precession, see equation (25) in page 11 and Table 10 in page 50 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract
    corrections = Ry(phi_y(t, consts))*Rx(phi_x(t, consts))
    # Nutation matrix, see equations (22)-(23) in page 11 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract
    N_iau80 = nutation_iau80(t, consts.nut_Δψ, consts.nut_Δϵ)
    return N_iau80*corrections*P_iau7680
end

@doc raw"""
    c2t_jpl_de430!(t, rotatBuf, consts::DE430Constants)

Inplace version of ``c2t_jpl_de430(t, consts)``; see [`c2t_jpl_de430`](@ref) for details.
"""
function c2t_jpl_de430!(t::Taylor1{T}, rotatBuf::RetAlloc{Taylor1{T}},
        consts::DE430Constants) where {T<:Real}
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
    zeta!(zt, t, auxs_angles, consts)  # zeta(t)
    Theta!(Tt, t, auxs_angles, consts) # Theta(t)
    Zeta!(Zt, t, auxs_angles, consts)
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
    matmul!(Raux, Rzz, RyT)
    matmul!(P_iau7680, Raux, RzZ)
    # Linear corrections to precession, see equation (25) in page 11 and Table 10 in page 50 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract
    # corrections = Ry(phi_y(t))*Rx(phi_x(t))
    phi_y!(phiy, t, consts)
    phi_x!(phix, t, consts)
    Ry!(Ryphiy, phiy, auxsR)
    Rx!(Rxphix, phix, auxsR)
    # corrections = Ryphiy*Rxphix
    matmul!(corrections, Ryphiy, Rxphix)
    # Nutation matrix, see equation (5-152) in page 5-60 of https://doi.org/10.1002/0471728470
    # N_iau80 = nutation_iau80(t)
    nutation_iau80!(N_iau80, t, auxNut, auxsR, view(rotatBuf.v2, 9:13), consts.nut_Δψ,
        consts.nut_Δϵ)
    # Raux = N_iau80*corrections
    # res = Raux*P_iau7680
    matmul!(Raux,  N_iau80, corrections)
    matmul!(res,  Raux, P_iau7680)
    return nothing
end

@doc raw"""
    allocate_c2t_jpl_de430(t)

Returns a `RetAlloc{Taylor1{T}}` object with objects that can be used to recicle memory for
[`c2t_jpl_de430!`](@ref)) and [`t2c_jpl_de430!`](@ref)) methods. Minimum computations
are performed here.
"""
function allocate_c2t_jpl_de430(t::Taylor1{T}) where {T<:Real}
    t_cy = t/(100yr)
    # Angles zeta(t), Theta(t), Zeta(t), negatives, and allocations
    zt = zero(t) # zt = zeta(t)
    auxhorner = zero(t)
    angle_A = zero(t)
    auxs_angles = [t_cy, auxhorner, angle_A] # auxs for angles
    nzt = -zt
    Tt = zero(t) # Tt = Theta(t)
    Zt = zero(t) # Zt = Zeta(t)
    nZt = -Zt
    # Allocations for calculations involving Rx, Ry and Rz
    one_alpha = one(nzt)
    zero_alpha = zero(nzt)
    sin_alpha, cos_alpha = sincos(nzt)
    auxsR = [one_alpha, zero_alpha, sin_alpha, cos_alpha]
    Rzz = Rz(nzt)
    RyT = zero.(Rzz)
    RzZ = zero.(Rzz)
    # Allocations for precession matrix
    Raux = zero.(Rzz)
    P_iau7680 = zero.(Rzz)
    # Angles and allocations for linear corrections of precession
    phiy = zero(t)
    phix = zero(t)
    Ryphiy = zero.(Rzz)
    Rxphix = zero.(Rzz)
    # Allocations for linear corrections to precession
    corrections = zero.(Rzz)
    # Allocations for nutation matrix
    ϵ0 = zero(t)        # Mean obliquity (rad)
    Δϵ = zero(t)        # Nutation in obliquity (rad)
    Δψ = zero(t)        # Nutation in longitude (rad)
    nΔψ = zero(t)       # Negative nutation in longitude (rad)
    nϵ = zero(t)        # ϵ = -(ϵ0 + Δϵ)
    omega_t_cy = zero(t) # Ω(t)
    auxNut = [ϵ0, Δϵ, Δψ, nΔψ, nϵ, omega_t_cy, t_cy, auxhorner, angle_A, sin_alpha, cos_alpha]
    Rxnϵ = zero.(Rzz)
    RznΔψ = zero.(Rzz)
    Rxϵ0 = zero.(Rzz)
    res = zero.(Rzz) # c2t_jpl_de430
    N_iau80 = zero.(Rzz)
    # Returned RetAlloc{Taylor1{T}} object
    rotatBuf = RetAlloc{Taylor1{T}}(
        [zt, nzt, Tt, Zt, nZt, phiy, phix],
        [auxs_angles, auxsR, auxNut],
        [Rzz, RyT, RzZ, P_iau7680, Ryphiy, Rxphix, corrections, N_iau80, Rxnϵ, RznΔψ, Rxϵ0, res, Raux],
        [Array{Taylor1{T}}(undef, 0, 0, 0)],
        [Array{Taylor1{T}}(undef, 0, 0, 0, 0)])
    return rotatBuf
end
