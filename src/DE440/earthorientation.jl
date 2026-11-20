# This file is part of the PlanetaryEphemeris.jl package; MIT licensed

# JPL-DE-440/441 model for the orientation of the Earth

# From JPL DE440 documentation (Park et al., 2021):
# Only the long-term change of the Earth orientation is modeled in the ephemeris
# integration. The Earth orientation model used for DE440 and DE441 is based on a
# long-term precession model (Vondrák et al. 2011) and a modified nutation model based
# on the IAU 1980 precession model including only terms with a period of 18.6 years.

# The IAU 1980 nutation and the mean obliquity are in src/orientation.jl.

# Long-term precession model of Vondrák et al. (2011), see
# https://doi.org/10.1051/0004-6361/201117274. Coefficients in arcseconds, as in the
# Fortran code of its Appendix A except for the cosine coefficient of Q_A with period
# 882 cy: 198.296701 as in SOFA, instead of 198.296071. This moves the equinox by less
# than 2 mas, not the pole.

# Obliquity at J2000.0
const ltp_eps0 = 84381.406
# Polynomial coefficients of P_A (first row) and Q_A (second row) in powers of T
const ltp_pqpol = [ 5851.607687 -0.1189000 -0.00028913  0.000000101;
                   -1600.886300  1.1689818 -0.00000020 -0.000000437]
# Periodic terms of P_A and Q_A: period [cy], cos coefficients of P_A and Q_A, and sin
# coefficients of P_A and Q_A
const ltp_pqper = [ 708.15 -5486.751211 -684.661560   667.666730 -5523.863691;
                   2309.00   -17.127623 2446.283880 -2354.886252  -549.747450;
                   1620.00  -617.517403  399.671049  -428.152441  -310.998056;
                    492.20   413.442940 -356.652376   376.202861   421.535876;
                   1183.00    78.614193 -186.387003   184.778874   -36.776172;
                    622.00  -180.732815 -316.800070   335.321713  -145.278396;
                    882.00   -87.676083  198.296701  -185.138669   -34.744450;
                    547.00    46.140315  101.135679  -120.972830    22.885731]
# Polynomial coefficients of X_A (first row) and Y_A (second row) in powers of T
const ltp_xypol = [  5453.282155  0.4252841 -0.00037173 -0.000000152;
                   -73750.930350 -0.7675452 -0.00018725  0.000000231]
# Periodic terms of X_A and Y_A: period [cy], cos coefficients of X_A and Y_A, and sin
# coefficients of X_A and Y_A
const ltp_xyper = [ 256.75  -819.940624 75004.344875 81491.287984  1558.515853;
                    708.15 -8444.676815   624.033993   787.163481  7774.939698;
                    274.20  2600.009459  1251.136893  1251.296102 -2219.534038;
                    241.45  2755.175630 -1102.212834 -1257.950837 -2523.969396;
                   2309.00  -167.659835 -2660.664980 -2966.799730   247.850422;
                    492.20   871.855056   699.291817   639.744522  -846.485643;
                    396.10    44.769698   153.167220   131.600209 -1393.124055;
                    288.90  -512.313065  -950.865637  -445.040117   368.526116;
                    231.10  -819.415595   499.754645   584.522874   749.045012;
                   1610.00  -538.071099  -145.188210   -89.756563   444.704518;
                    620.00  -189.793622   558.116553   524.429630   235.934465;
                    157.87  -402.922932   -23.923029   -13.549067   374.049623;
                    220.30   179.516345  -165.405086  -210.157124  -171.330180;
                   1200.00    -9.814756     9.344131   -44.919798   -22.899655]

# Periodic plus polynomial terms (arcseconds) of a pair of long-term precession
# parameters at `T` Julian centuries from J2000.0
function _ltp_terms(T, pol::AbstractMatrix, per::AbstractMatrix)
    a, b = zero(T), zero(T)
    for i in axes(per, 1)
        s, c = sincos(2π*T/per[i, 1])
        a += c*per[i, 2] + s*per[i, 4]
        b += c*per[i, 3] + s*per[i, 5]
    end
    Tn = one(T)
    for j in axes(pol, 2)
        a += pol[1, j]*Tn
        b += pol[2, j]*Tn
        Tn *= T
    end
    return a, b
end

@doc raw"""
    ltp_ecliptic_pole(t::Real)

Return the ecliptic pole unit vector ``\mathbf{k}`` with respect to the J2000.0 mean
equator and equinox, according to the long-term precession model of Vondrák et al. (2011)
```math
\mathbf{k} = \left(P_A, -Q_A\cos\epsilon_0 - Z\sin\epsilon_0,
-Q_A\sin\epsilon_0 + Z\cos\epsilon_0\right), \quad Z = \sqrt{1 - P_A^2 - Q_A^2},
```
where ``t`` is the TDB time in Julian days from J2000.0 and ``\epsilon_0 = 84381''.406``.

See equation (8), Table 1 and Appendix A.1 of https://doi.org/10.1051/0004-6361/201117274.

See also [`ltp_equator_pole`](@ref) and [`ltp_precession`](@ref).
"""
ltp_ecliptic_pole(t::Real) = _ltp_ecliptic_pole(t)

# Generic methods, also used to test the in-place Taylor1 methods below
function _ltp_ecliptic_pole(t)
    p, q = _ltp_terms(t/36525, ltp_pqpol, ltp_pqper)
    p, q = deg2rad(p/3600), deg2rad(q/3600)
    z = sqrt(1 - p^2 - q^2)
    s, c = sincos(deg2rad(ltp_eps0/3600))
    return [p, -q*c - z*s, -q*s + z*c]
end

@doc raw"""
    ltp_equator_pole(t::Real)

Return the mean equator pole unit vector ``\mathbf{n} = (X_A, Y_A, \sqrt{1 - X_A^2 -
Y_A^2})`` with respect to the J2000.0 mean equator and equinox, according to the
long-term precession model of Vondrák et al. (2011), where ``t`` is the TDB time in
Julian days from J2000.0.

See equation (9), Table 2 and Appendix A.2 of https://doi.org/10.1051/0004-6361/201117274.

See also [`ltp_ecliptic_pole`](@ref) and [`ltp_precession`](@ref).
"""
ltp_equator_pole(t::Real) = _ltp_equator_pole(t)

function _ltp_equator_pole(t)
    x, y = _ltp_terms(t/36525, ltp_xypol, ltp_xyper)
    x, y = deg2rad(x/3600), deg2rad(y/3600)
    return [x, y, sqrt(1 - x^2 - y^2)]
end

@doc raw"""
    ltp_precession(t::Real)

Return the long-term precession matrix from the J2000.0 mean equator and equinox to the
mean equator and equinox of date, according to Vondrák et al. (2011)
```math
P = \left(
\begin{array}{c}
    \mathbf{e}^T \\
    (\mathbf{n} \times \mathbf{e})^T \\
    \mathbf{n}^T
\end{array}
\right), \quad
\mathbf{e} = \frac{\mathbf{n} \times \mathbf{k}}{|\mathbf{n} \times \mathbf{k}|},
```
where ``t`` is the TDB time in Julian days from J2000.0, and ``\mathbf{k}`` and
``\mathbf{n}`` are the ecliptic and mean equator poles.

``\textbf{Note:}`` this is the matrix ``V`` of equation (26) of Park et al. (2021). Since it
maps J2000.0 to date, the Earth pole in their equation (25) is
``\mathbf{p}_E = V^T R_x(-\Phi_x) R_y(-\Phi_y) \mathbf{p}_d``.

See equation (23) and Appendix A.3 of https://doi.org/10.1051/0004-6361/201117274, and
equations (25)-(26) of https://doi.org/10.3847/1538-3881/abd414.

See also [`ltp_ecliptic_pole`](@ref) and [`ltp_equator_pole`](@ref).
"""
ltp_precession(t::Real) = _ltp_precession(t)

function _ltp_precession(t)
    n = _ltp_equator_pole(t)
    k = _ltp_ecliptic_pole(t)
    # Equinox of date
    e = cross(n, k)
    e = e/sqrt(e[1]^2 + e[2]^2 + e[3]^2)
    return permutedims(hcat(e, cross(n, e), n))
end

# In-place Taylor1 methods

# Auxiliaries for _ltp_terms!: T, 2πT/period, sin, cos, polynomial (4 coefficients),
# its value and an auxiliary for _horner!
_allocate_ltp_terms_auxs(t::Taylor1{T}) where {T<:Real} =
    [zero(t), zero(t), zero(t), zero(t), Taylor1(zeros(T, 4)), zero(t), zero(t)]

# In-place version of _ltp_terms, with `t` the TDB time in Julian days from J2000.0
function _ltp_terms!(a::Taylor1{T}, b::Taylor1{T}, t::Taylor1{T}, pol::AbstractMatrix,
        per::AbstractMatrix, auxs) where {T<:Real}
    t_cy, arg, s, c, poly, pval, aux = auxs
    for ord in eachindex(t)
        t_cy[ord] = t[ord]/36525
        a[ord] = zero(t[ord])
        b[ord] = zero(t[ord])
    end
    # Periodic terms
    for i in axes(per, 1)
        for ord in eachindex(t)
            arg[ord] = (2π*t_cy[ord])/per[i, 1]
        end
        for ord in eachindex(t)
            TaylorSeries.sincos!(s, c, arg, ord)
            a[ord] += c[ord]*per[i, 2] + s[ord]*per[i, 4]
            b[ord] += c[ord]*per[i, 3] + s[ord]*per[i, 5]
        end
    end
    # Polynomial terms
    for (res, row) in ((a, 1), (b, 2))
        for j in axes(pol, 2)
            poly.coeffs[j] = pol[row, j]
        end
        for ord in eachindex(t)
            pval[ord] = zero(t[ord])
        end
        TaylorSeries._horner!(pval, poly, t_cy, aux)
        for ord in eachindex(t)
            res[ord] += pval[ord]
        end
    end
    return nothing
end

# Auxiliaries for ltp_ecliptic_pole! and ltp_equator_pole!: two direction cosines, their
# squares, 1 minus the sum of the squares, its square root, an auxiliary for sqrt! and
# the auxiliaries for _ltp_terms!
_allocate_ltp_pole_auxs(t::Taylor1{T}) where {T<:Real} =
    [[zero(t) for _ in 1:7]; _allocate_ltp_terms_auxs(t)]

# w ← 1 - u^2 - v^2 and z ← sqrt(w)
function _ltp_unit_z!(z, u, v, auxs)
    u2, v2, w, aux = auxs[3], auxs[4], auxs[5], auxs[7]
    TaylorSeries.mul!(u2, u, u)
    TaylorSeries.mul!(v2, v, v)
    for ord in eachindex(w)
        w[ord] = ((ord == 0 ? one(w[ord]) : zero(w[ord])) - u2[ord]) - v2[ord]
    end
    for ord in eachindex(w)
        TaylorSeries.sqrt!(z, w, aux, ord)
    end
    return nothing
end

@doc raw"""
    ltp_ecliptic_pole!(k, t, auxs)

In-place version of ``ltp_ecliptic_pole(t)``; see [`ltp_ecliptic_pole`](@ref) for details.
"""
function ltp_ecliptic_pole!(k::AbstractVector{Taylor1{T}}, t::Taylor1{T},
        auxs) where {T<:Real}
    p, q, z = auxs[1], auxs[2], auxs[6]
    _ltp_terms!(p, q, t, ltp_pqpol, ltp_pqper, view(auxs, 8:14))
    for ord in eachindex(t)
        p[ord] = deg2rad(p[ord]/3600)
        q[ord] = deg2rad(q[ord]/3600)
    end
    _ltp_unit_z!(z, p, q, auxs)
    s, c = sincos(deg2rad(ltp_eps0/3600))
    for ord in eachindex(t)
        k[1][ord] = p[ord]
        k[2][ord] = -q[ord]*c - z[ord]*s
        k[3][ord] = -q[ord]*s + z[ord]*c
    end
    return nothing
end

@doc raw"""
    ltp_equator_pole!(n, t, auxs)

In-place version of ``ltp_equator_pole(t)``; see [`ltp_equator_pole`](@ref) for details.
"""
function ltp_equator_pole!(n::AbstractVector{Taylor1{T}}, t::Taylor1{T},
        auxs) where {T<:Real}
    x, y = auxs[1], auxs[2]
    _ltp_terms!(x, y, t, ltp_xypol, ltp_xyper, view(auxs, 8:14))
    for ord in eachindex(t)
        x[ord] = deg2rad(x[ord]/3600)
        y[ord] = deg2rad(y[ord]/3600)
        n[1][ord] = x[ord]
        n[2][ord] = y[ord]
    end
    _ltp_unit_z!(n[3], x, y, auxs)
    return nothing
end

# Auxiliaries for ltp_precession!: equator pole (3), ecliptic pole (3), n × k (3), two
# products, |n × k|^2, |n × k|, an auxiliary for sqrt! and the auxiliaries for
# ltp_ecliptic_pole! and ltp_equator_pole!
_allocate_ltp_precession_auxs(t::Taylor1{T}) where {T<:Real} =
    [[zero(t) for _ in 1:14]; _allocate_ltp_pole_auxs(t)]

# c ← a1*b2 - a2*b1, with auxiliaries m1 and m2
function _ltp_cross!(c, a1, b2, a2, b1, m1, m2)
    TaylorSeries.mul!(m1, a1, b2)
    TaylorSeries.mul!(m2, a2, b1)
    for ord in eachindex(c)
        c[ord] = m1[ord] - m2[ord]
    end
    return nothing
end

@doc raw"""
    ltp_precession!(P, t, auxs)

In-place version of ``ltp_precession(t)``; see [`ltp_precession`](@ref) for details.
"""
function ltp_precession!(P::AbstractMatrix{Taylor1{T}}, t::Taylor1{T},
        auxs) where {T<:Real}
    n, k, e = view(auxs, 1:3), view(auxs, 4:6), view(auxs, 7:9)
    m1, m2, r2, r, aux = auxs[10], auxs[11], auxs[12], auxs[13], auxs[14]
    pole_auxs = view(auxs, 15:28)
    ltp_equator_pole!(n, t, pole_auxs)
    ltp_ecliptic_pole!(k, t, pole_auxs)
    # Equinox of date: e = n × k / |n × k|
    _ltp_cross!(e[1], n[2], k[3], n[3], k[2], m1, m2)
    _ltp_cross!(e[2], n[3], k[1], n[1], k[3], m1, m2)
    _ltp_cross!(e[3], n[1], k[2], n[2], k[1], m1, m2)
    TaylorSeries.mul!(m1, e[1], e[1])
    TaylorSeries.mul!(m2, e[2], e[2])
    for ord in eachindex(t)
        r2[ord] = m1[ord] + m2[ord]
    end
    TaylorSeries.mul!(m1, e[3], e[3])
    for ord in eachindex(t)
        r2[ord] += m1[ord]
    end
    for ord in eachindex(t)
        TaylorSeries.sqrt!(r, r2, aux, ord)
    end
    for j in 1:3, ord in eachindex(t)
        TaylorSeries.div!(P[1, j], e[j], r, ord)
    end
    # Second row: n × e
    _ltp_cross!(P[2, 1], n[2], P[1, 3], n[3], P[1, 2], m1, m2)
    _ltp_cross!(P[2, 2], n[3], P[1, 1], n[1], P[1, 3], m1, m2)
    _ltp_cross!(P[2, 3], n[1], P[1, 2], n[2], P[1, 1], m1, m2)
    # Third row: n
    for j in 1:3, ord in eachindex(t)
        P[3, j][ord] = n[j][ord]
    end
    return nothing
end

# Allocating Taylor1 methods: fall back to the in-place methods

for (f, f!) in ((:ltp_ecliptic_pole, :ltp_ecliptic_pole!),
        (:ltp_equator_pole, :ltp_equator_pole!))
    @eval function $f(t::Taylor1{T}) where {T<:Real}
        res = [zero(t) for _ in 1:3]
        $(f!)(res, t, _allocate_ltp_pole_auxs(t))
        return res
    end
end

function ltp_precession(t::Taylor1{T}) where {T<:Real}
    P = [zero(t) for _ in 1:3, _ in 1:3]
    ltp_precession!(P, t, _allocate_ltp_precession_auxs(t))
    return P
end

# Offsets of the EME2000 pole from the ICRF pole (ROTEX and ROTEY in the DE440 technical
# comments, converted from arcseconds to degrees); see equation (25) of Park et al. (2021)
const phi_x_de440 = 5.3439916044148049E-03/3600   # x-axis rotation (degrees)
const phi_y_de440 = -1.7128824840317532E-02/3600  # y-axis rotation (degrees)

@doc raw"""
    c2t_jpl_de440(t::Real)

Return the rotation matrix from the inertial (ICRF) frame (celestial) to the Earth's true
equator of date frame (terrestrial), according to the JPL DE440/441 Earth orientation
model
```math
\texttt{c2t_jpl_de440}(t) = N\times C\times P,
```
where ``t`` is the TDB time in Julian days from J2000.0, ``P`` is the long-term
precession matrix, ``C = R_y(\Phi_y)R_x(\Phi_x)`` contains the constant offsets of the
EME2000 pole from the ICRF pole and ``N`` is the nutation matrix. The third row is the
Earth's pole ``\mathbf{p}_\mathrm{E}`` in the inertial frame.

``\textbf{Note:}`` the Earth's rotation is not modeled: only the direction of the pole is
needed for the Earth's zonal harmonics.

See equations (20)-(26) in page 5 of https://doi.org/10.3847/1538-3881/abd414.

See also [`ltp_precession`](@ref), [`nutation_iau80`](@ref) and [`t2c_jpl_de440`](@ref).
"""
function c2t_jpl_de440(t::Real)
    # Long-term precession matrix
    P_ltp = ltp_precession(t)
    # Offsets of the EME2000 pole from the ICRF pole
    offsets = Ry(deg2rad(phi_y_de440))*Rx(deg2rad(phi_x_de440))
    # Nutation matrix
    N_iau80 = nutation_iau80(t)
    return N_iau80*offsets*P_ltp
end

@doc raw"""
    t2c_jpl_de440(t)

Return the rotation matrix from the Earth's true equator of date frame (terrestrial) to
the inertial (ICRF) frame (celestial), according to the JPL DE440/441 Earth orientation
model
```math
\texttt{t2c_jpl_de440}(t) = \texttt{c2t_jpl_de440}^T(t),
```
where ``t`` is the TDB time in Julian days from J2000.0.

See also [`c2t_jpl_de440`](@ref).
"""
t2c_jpl_de440(t) = transpose(c2t_jpl_de440(t))
