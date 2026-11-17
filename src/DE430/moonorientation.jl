# This file is part of the PlanetaryEphemeris.jl package; MIT licensed

# JPL DE430/431 model for the orientation of the Moon: moment of inertia of the lunar
# mantle; see Section IV.B of Folkner et al. (2014)

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
