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
function ltp_ecliptic_pole(t::Real)
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
function ltp_equator_pole(t::Real)
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
function ltp_precession(t::Real)
    n = ltp_equator_pole(t)
    k = ltp_ecliptic_pole(t)
    # Equinox of date
    e = cross(n, k)
    e = e/norm(e)
    return permutedims(hcat(e, cross(n, e), n))
end
