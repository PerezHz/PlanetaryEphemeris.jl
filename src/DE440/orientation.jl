# This file is part of the PlanetaryEphemeris.jl package; MIT licensed

# JPL-DE-440/441 model for the orientation of the Earth

# From JPL DE440 documentation (Park et al., 2021):
# Only the long-term change of the Earth orientation is modeled in the ephemeris
# integration. The Earth orientation model used for DE440 and DE441 is based on a
# long-term precession model (Vondrák et al. 2011) and a modified nutation model based
# on the IAU 1980 precession model including only terms with a period of 18.6 years.

# Sum of the periodic and polynomial terms of the long-term precession model, in
# arcseconds, at `t_cy` Julian centuries from J2000.0
function _ltp_terms(t_cy, pol::AbstractMatrix, per::AbstractMatrix)
    w = 2π*t_cy
    a = zero(t_cy)
    b = zero(t_cy)
    # Periodic terms
    for i in axes(per, 1)
        s, c = sincos(w/per[i, 1])
        a += c*per[i, 2] + s*per[i, 4]
        b += c*per[i, 3] + s*per[i, 5]
    end
    # Polynomial terms
    tn = one(t_cy)
    for i in axes(pol, 2)
        a += pol[1, i]*tn
        b += pol[2, i]*tn
        tn *= t_cy
    end
    return a, b
end

@doc raw"""
    ltp_ecliptic_pole(t, consts::DE440Constants)

Return the ecliptic pole unit vector ``\mathbf{k}`` in the J2000.0 mean equator and
equinox frame, according to the long-term precession model of Vondrák et al. (2011),
where ``t`` is the TDB time in Julian days from J2000.0.

See https://doi.org/10.1051/0004-6361/201117274 and function `eraLtpecl` of ERFA
(https://github.com/liberfa/erfa).

See also [`ltp_equator_pole`](@ref) and [`ltp_precession`](@ref).
"""
function ltp_ecliptic_pole(t, consts::DE440Constants)
    # P_A and Q_A (radians)
    p, q = _ltp_terms(t/36525, consts.ltp_pqpol, consts.ltp_pqper)
    p, q = deg2rad(p/3600), deg2rad(q/3600)
    # Ecliptic pole vector
    z = sqrt(1 - p^2 - q^2)
    s, c = sincos(deg2rad(consts.ltp_eps0/3600))
    return [p, -q*c - z*s, -q*s + z*c]
end

@doc raw"""
    ltp_equator_pole(t, consts::DE440Constants)

Return the mean equator pole unit vector ``\mathbf{n}`` in the J2000.0 mean equator and
equinox frame, according to the long-term precession model of Vondrák et al. (2011),
where ``t`` is the TDB time in Julian days from J2000.0.

See https://doi.org/10.1051/0004-6361/201117274 and function `eraLtpequ` of ERFA
(https://github.com/liberfa/erfa).

See also [`ltp_ecliptic_pole`](@ref) and [`ltp_precession`](@ref).
"""
function ltp_equator_pole(t, consts::DE440Constants)
    # X_A and Y_A (radians)
    x, y = _ltp_terms(t/36525, consts.ltp_xypol, consts.ltp_xyper)
    x, y = deg2rad(x/3600), deg2rad(y/3600)
    # Equator pole vector
    return [x, y, sqrt(1 - x^2 - y^2)]
end

@doc raw"""
    ltp_precession(t, consts::DE440Constants)

Return the long-term precession matrix from the J2000.0 mean equator and equinox frame
to the mean equator and equinox of date
```math
V^T = \left[
\begin{array}{c}
    (\mathbf{n} \times \mathbf{k}) / |\mathbf{n} \times \mathbf{k}| \\
    \mathbf{n} \times (\mathbf{n} \times \mathbf{k}) / |\mathbf{n} \times (\mathbf{n} \times \mathbf{k})| \\
    \mathbf{n} / |\mathbf{n}| \\
\end{array}
\right],
```
where ``t`` is the TDB time in Julian days from J2000.0, ``\mathbf{k}`` is the ecliptic
pole and ``\mathbf{n}`` is the mean equator pole, both in the J2000.0 frame.

See equation (26) in page 5 of https://doi.org/10.3847/1538-3881/abd414 and function
`eraLtp` of ERFA (https://github.com/liberfa/erfa).

See also [`ltp_ecliptic_pole`](@ref) and [`ltp_equator_pole`](@ref).
"""
function ltp_precession(t, consts::DE440Constants)
    n = ltp_equator_pole(t, consts)
    k = ltp_ecliptic_pole(t, consts)
    # Equinox of date
    eqx = cross(n, k)
    eqx = eqx / sqrt(eqx[1]^2 + eqx[2]^2 + eqx[3]^2)
    v = cross(n, eqx)
    return permutedims(hcat(eqx, v, n))
end

@doc raw"""
    c2t_jpl_de440(t, consts::DE440Constants)

Return the matrix for celestial-to-terrestrial coordinate transformation, according to
the JPL DE440/441 Earth orientation model
```math
\texttt{c2t_jpl_de440}(t) = N \times R_y(\Phi_y) R_x(\Phi_x) \times V^T,
```
where ``t`` is the TDB time in Julian days from J2000.0 (a `Taylor1`), ``V^T`` is the
long-term precession matrix, ``(\Phi_x, \Phi_y)`` are the offsets of the EME2000 pole from
the ICRF pole and ``N`` is the nutation matrix. The Earth pole unit vector in the
inertial frame is the third row of this matrix,
``\mathbf{p}_E = V R_x(-\Phi_x) R_y(-\Phi_y) \mathbf{p}_d``.

See equations (20)-(26) in page 5 of https://doi.org/10.3847/1538-3881/abd414.

See also [`ltp_precession`](@ref), [`nutation_iau80`](@ref) and [`t2c_jpl_de440`](@ref).
"""
function c2t_jpl_de440(t, consts::DE440Constants)
    # Long-term precession matrix
    P_ltp = ltp_precession(t, consts)
    # Offsets of the EME2000 pole from the ICRF pole
    offsets = Ry(deg2rad(consts.ROTEY/3600))*Rx(deg2rad(consts.ROTEX/3600))
    # Nutation matrix, see equations (20)-(24) in page 5 of https://doi.org/10.3847/1538-3881/abd414
    N_iau80 = nutation_iau80(t, consts.nut_Δψ, consts.nut_Δϵ)
    return N_iau80*offsets*P_ltp
end

@doc raw"""
    t2c_jpl_de440(t, consts::DE440Constants)

Return the matrix for terrestrial-to-celestial coordinate transformation, according to
the JPL DE440/441 Earth orientation model
```math
\texttt{t2c_jpl_de440}(t) = \texttt{c2t_jpl_de440}^T(t),
```
where ``t`` is the TDB time in Julian days from J2000.0 (a `Taylor1`).

See also [`c2t_jpl_de440`](@ref).
"""
t2c_jpl_de440(t, consts::DE440Constants) = transpose(c2t_jpl_de440(t, consts))
