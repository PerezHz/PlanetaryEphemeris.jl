# Solar System ( JPL DE430/431) dynamical model
# Bodies considered in the model are: the Sun, the eight planets, the Moon and
# the 343 main-belt asteroids included in the JPL DE430 ephemeris.
# Effects considered are:
# - Post-Newtonian point-mass accelerations between all bodies,
# - Figure-effects (oblateness) of the Earth (J2 and J3)
# - J2 effect of the Sun
# - J2 and J3 effect of the Moon
# - Kinematic model for the precession and nutation of the Earth's orientation (IAU 1976/1980 Earth orientation model)
# - Kinematic model for the Moons's orientation (Seidelmann et al., 2006)
# - Tidal secular acceleration of Moon due to rides raised on Earth by both the Moon and the Sun

@doc """
    ordpres_differentiate(a::Taylor1)

Returns the derivative of `a`, but preserving the order/degree of `a`. In comparison,
`TaylorSeries.differentiate` returns the returns a `Taylor1` object with one order/degree
less than the one of `a`.

See also [`TaylorSeries.differentiate`](@ref).
"""
function ordpres_differentiate(a::Taylor1{T}) where {T}
    res = zero(a)
    for ord in eachindex(res)
        TaylorSeries.differentiate!(res, a, ord)
    end
    return res
end

@doc raw"""
    special_eval(x::Vector{Taylor1{T}}, t::Taylor1{T}) where {T <: Number}

Evaluate each element of `x` at time `t`.
"""
function special_eval(x::Vector{Taylor1{T}}, t::Taylor1{T}) where {T <: Number}
    res = Vector{Taylor1{T}}(undef, length(x))
    for i in eachindex(res)
        res[i] = x[i](t)
    end
    return res
end


@doc raw"""
    NBP_pN_A_J23E_J23M_J2S!(dq, q, params, t)

Solar System (JPL DE430/431) dynamical model. Bodies considered in the model are: the Sun,
the eight planets, the Moon and the 343 main-belt asteroids included in the JPL DE430
ephemeris. Effects considered are:

- Post-Newtonian point-mass accelerations between all bodies: see equation (35) in page 7 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract
```math
\begin{align*}
    \mathbf{a}_i & = \sum_{j\neq i}\frac{\mu_j(\mathbf{r}_j - \mathbf{r}_i)}{r_{ij}^3}
    \left\lbrace
    1 - \frac{4}{c^2}\sum_{l\neq i}\frac{\mu_l}{r_{il}} - \frac{1}{c^2}\sum_{k\neq j}\frac{\mu_k}{r_{jk}} + \left(\frac{\dot{s}_i}{c}\right)^2 + 2\left(\frac{\dot{s}_j}{c}\right)^2
    - \frac{4}{c^2}\mathbf{v}_i\cdot\mathbf{v}_j - \frac{3}{2c^2}\left[\frac{(\mathbf{r}_i - \mathbf{r}_j)\cdot\mathbf{v}_j}{r_{ij}}\right]^2 + \frac{1}{2c^2}(\mathbf{r}_j - \mathbf{r}_i)\cdot\mathbf{a}_j
    \right\rbrace \\
    & \hspace{0.5cm} + \frac{1}{c^2}\sum_{j\neq i} \frac{\mu_j}{r_{ij}^3}[(\mathbf{r}_i - \mathbf{r}_j)\cdot(4\mathbf{v}_i - 3\mathbf{v}_j)](\mathbf{v}_i - \mathbf{v}_j)
    + \frac{7}{2c^2}\sum_{j\neq i}\frac{\mu_j\mathbf{a}_j}{r_{ij}},
\end{align*}
```
where ``\mathbf{v}_i = \dot{\mathbf{r}}_i``, ``\dot{s}_i = ||\mathbf{v}_i||^2`` and
``\mathbf{a}_i = \ddot{\mathbf{r}}_i``.

- Figure-effects (oblateness) of the Earth (``J_2`` and ``J_3``),

- ``J_2`` effect of the Sun and

- ``J_2`` and ``J_3`` effect of the Moon: see equation (28) in page 13 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract and equations (173) and (174) in page 33 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract
```math
\left[
\begin{array}{c}
    \ddot{\xi} \\
    \ddot{\eta} \\
    \ddot{\zeta} \\
\end{array}
\right]
=
-\frac{GM}{r^2}
\left\lbrace
\sum_{n=2}^{n_1} J_n\left(\frac{R}{r}\right)^n
\left[
\begin{array}{c}
    (n+1)P_n(\sin\phi) \\
    0 \\
    -\cos\phi P_n'(\sin\phi) \\
\end{array}
\right]
+\sum_{n=2}^{n_2} \left(\frac{R}{r}\right)^n\sum_{m=1}^n
\left[
\begin{array}{c}
    -(n+1)P_n^m(\sin\phi)[+C_{nm}\cos m\lambda + S_{nm}\sin m\lambda] \\
    m\sec\phi P_n^m(\sin\phi)[-C_{nm}\sin m\lambda + S_{nm}\cos m\lambda] \\
    \cos\phi P'\,_n^m(\sin\phi)[+C_{nm}\cos m\lambda + S_{nm}\sin m\lambda] \\
\end{array}
\right]
\right\rbrace,
```
where ``r`` is the center-of-mass separation between the two bodies; ``n_1`` and ``n_2`` are the maximum degrees
of the zonal and tesseral expansions, respectively; ``P_n(\sin\phi)`` is the Legendre polynomial of degree ``n``,
``P_n^m(\sin\phi)`` is the associated Legendre function of degree ``n`` and order ``m``, ``J_n`` is the zonal
harmonic coefficient for the extended body; ``C_{nm}``, ``S_{nm}`` are the tesseral harmonic coefficients for
the extended body, ``R`` is the equatorial radius of the extended body, ``\phi`` is the latitude of the point
mass relative to the body-fixed coordinate system in which the harmonics are expressed; and ``\lambda`` is the
east longitude of the point mass in the same body-fixed coordinate system. The primes denote differentiation
with respect to the argument ``\sin\phi``. The accelerations are transformed into the inertial frame by
application of the appropriate rotation matrix.

- Kinematic model for the precession and nutation of the Earth's orientation (IAU 1976/1980 Earth orientation model): see [`c2t_jpl_de430`](@ref).

- Kinematic model for the Moon's orientation (Seidelmann et al., 2006): see equations (14)-(15) in page 9 and equations (34)-(35) in page 16 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract.

``\textbf{Lunar mantle}``
```math
\begin{align*}
    \dot{\phi}_m & = (\omega_{m,x}\sin\psi_m + \omega_{m,y}\cos\psi_m)/\sin\theta_m \\
    \dot{\theta}_m & = \omega_{m,x}\cos\psi_m - \omega_{m,y}\sin\psi_m \\
    \dot{\psi}_m & = \omega_{m,z} - \dot{\phi_m}\cos\theta_m
\end{align*}
```
and
```math
    \dot{\mathbf{\omega}}_m = \mathbf{I}_m^{-1}\left\lbrace
        \sum_{A\neq M} \mathbf{N}_{M,fig M- pmA} + \mathbf{N}_{M,figM-figE} - \dot{\mathbf{I}}_m\mathbf{\omega}_m - \mathbf{\omega}_m \times \mathbf{I}_m\mathbf{\omega}_m + \mathbf{N}_{cmb}
    \right\rbrace,
```
where ``(\phi_m, \theta_m, \psi_m)`` are the lunar mantle Euler angles, ``\mathbf{\omega}_m`` is the angular
velocity of the mantle expressed in the mantle frame, ``\mathbf{I}_m`` is the mantle moment of inertia,
``\mathbf{N}_{M,figM-pmA}`` is the torque on the lunar
mantle from the point mass of body ``A`` , ``\mathbf{N}_{M, figM - figE}`` is the torque on the mantle due to
the extended figure of the Moon interacting with the extended figure of the Earth, and ``\mathbf{N}_{cmb}``
is the torque due to interaction between the mantle and core.

``\textbf{Lunar core}``

```math
\begin{align*}
    \dot{\phi}_c & = \omega_{c,z}^\dagger - \dot{\psi}_c\cos\theta_c \\
    \dot{\theta}_c & = \omega_{c,x}^\dagger \\
    \dot{\psi}_c & = -\omega_{c,y}^\dagger / \sin\theta_c
\end{align*}
```
and
```math
    \dot{\mathbf{\omega}}_c = \mathbf{I}_c^{-1}\left\lbrace
        -\mathbf{\omega}_m \times \mathbf{I}_c\mathbf{\omega}_c - \mathbf{N}_{cmb}
    \right\rbrace,
```
where ``(\phi_c, \theta_c, \psi_c)`` are the lunar core Euler angles, ``\mathbf{\omega}_c^\dagger`` is the
angular velocity of the core expressed in a frame define by the intersection of the core equator with the
inertial ``XY`` plane, ``\mathbf{I}_c`` is the core moment of inertia; ``\mathbf{\omega}_m`` and ``\mathbf{\omega}_c``
are the mantle and core angular velocities in the mantle frame; and ``\mathbf{N}_{cmb}`` is the torque due
to interaction between the mantle and core.

""" NBP_pN_A_J23E_J23M_J2S!

function NBP_pN_A_J23E_J23M_J2S!(dq, q, params, t)
    # N: number of bodies
    # jd0: initial Julian date
    local N, jd0 = params
    local S = eltype(q)   # Type of positions/velocities components

    local zero_q_1 = zero(q[1])                  # Zero of type S
    local one_t = one(t)                         # One of the same type as time t
    local dsj2k = t+(jd0-J2000)                  # Days since J2000.0 (TDB)
    # Matrix elements of lunar mantle moment of inertia at time t-τ_M (without tidal distortion)
    # See equations (36) to (41) in pages 16-17 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract
    # ITM(q_del_τ_M, eulang_del_τ_M)
    local I_m_t = (ITM_und-I_c).*one_t           # Undistorted moment of inertia of the mantle, see equation (40)
    local dI_m_t = ordpres_differentiate.(I_m_t) # Time-derivative of lunar mantle I at time t-τ_M
    local inv_I_m_t = inv(I_m_t)                 # Inverse of lunar mantle I matrix at time t-τ_M
    local I_c_t = I_c.*one_t                     # Lunar core I matrix, see equation (39)
    local inv_I_c_t = inv(I_c_t)                 # Inverse of lunar core I matrix
    local I_M_t = I_m_t+I_c_t                    # Total I matrix (mantle + core)

    #=
    Point-mass accelerations
    See equation (35) in page 7 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract
    =#

    # Note: All the following arrays are declared here in order to help @taylorize work

    # Difference between two positions (\mathbf{r}_i - \mathbf{r}_j)
    X = Array{S}(undef, N, N)         # X-axis component
    Y = Array{S}(undef, N, N)         # Y-axis component
    Z = Array{S}(undef, N, N)         # Z-axis component

    # Distance between two positions r_{ij} = ||\mathbf{r}_i - \mathbf{r}_j||
    r_p2 = Array{S}(undef, N, N)      # r_{ij}^2
    r_p1d2 = Array{S}(undef, N, N)    # sqrt(r_p2) <-> r_{ij}
    r_p3d2 = Array{S}(undef, N, N)    # r_p2^1.5 <-> r_{ij}^3
    r_p7d2 = Array{S}(undef, N, N)    # r_p2^3.5 <-> r_{ij}^7

    # Newtonian accelerations \mathbf{a}_{i} = \sum_{i\neq j} mu_i * (\mathbf{r_i} - \mathbf{r_j}) / r_{ij}^3
    newtonX = Array{S}(undef, N)      # X-axis component
    newtonY = Array{S}(undef, N)      # Y-axis component
    newtonZ = Array{S}(undef, N)      # Z-axis component
    # Newtonian coefficient, i.e., mass parameter / distance^3 -> \mu_i / r_{ij}^3
    newtonianCoeff = Array{S}(undef, N, N)

    # Post-Newtonian stuff

    # Difference between two velocities (\mathbf{v}_i - \mathbf{v}_j)
    U = Array{S}(undef, N, N)         # X-axis component
    V = Array{S}(undef, N, N)         # Y-axis component
    W = Array{S}(undef, N, N)         # Z-axis component

    # Weighted difference between two velocities (4\mathbf{v}_i - 3\mathbf{v}_j)
    _4U_m_3X = Array{S}(undef, N, N)  # X-axis component
    _4V_m_3Y = Array{S}(undef, N, N)  # Y-axis component
    _4W_m_3Z = Array{S}(undef, N, N)  # Z-axis component

    # Product of velocity components
    UU = Array{S}(undef, N, N)        # v_{ix}v_{jx}
    VV = Array{S}(undef, N, N)        # v_{iy}v_{jy}
    WW = Array{S}(undef, N, N)        # v_{iz}v_{jz}

    # Newtonian potential of 1 body \mu_i / r_{ij}
    newtonian1b_Potential = Array{S}(undef, N, N)
    # Newtonian potential of N bodies
    # \sum_{i\neq l} \frac{\mu_i}{r_{il}}
    newtonianNb_Potential = Array{S}(undef, N)

    # Newtonian coefficient * difference between two positions, i.e.,
    # \mu_i * (\mathbf{r_i} - \mathbf{r_j}) / r_{ij}^3
    newton_acc_X = Array{S}(undef, N, N)   # X-axis component
    newton_acc_Y = Array{S}(undef, N, N)   # Y-axis component
    newton_acc_Z = Array{S}(undef, N, N)   # Z-axis component

    # Combinations of velocities
    v2 = Array{S}(undef, N)                # Velocity magnitude squared ||\mathbf{v}_i||^2
    _2v2 = Array{S}(undef, N, N)           # 2 * ||\mathbf{v_i}||^2
    vi_dot_vj = Array{S}(undef, N, N)      # Dot product of two velocities \mathbf{v}_i\cdot\mathbf{v}_j

    # Second term of equation (35) in page 7 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract

    # Second term without (\mathbf{v}_i - \mathbf{v}_j)
    pn2 = Array{S}(undef, N, N)            # \mu_i * [(\mathbf{r_i} - \mathbf{r_j})\cdot(4\mathbf{v_i} - 3\mathbf{v_j})]
    # Full second term
    U_t_pn2 = Array{S}(undef, N, N)        # X-axis component
    V_t_pn2 = Array{S}(undef, N, N)        # Y-axis component
    W_t_pn2 = Array{S}(undef, N, N)        # Z-axis component

    # Third term of equation (35) in page 7 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract

    # Third term without newtonian accelerations \mathbf{a}_i
    pn3 = Array{S}(undef, N, N)
    # Full third term of equation (35)
    pNX_t_pn3 = Array{S}(undef, N, N)      # X-axis component
    pNY_t_pn3 = Array{S}(undef, N, N)      # Y-axis component
    pNZ_t_pn3 = Array{S}(undef, N, N)      # Z-axis component

    # First term of equation (35) in page 7 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract

    _4ϕj = Array{S}(undef, N, N)            # 4*\sum term inside {}
    ϕi_plus_4ϕj = Array{S}(undef, N, N)     # 4*\sum + \sum terms inside {}
    sj2_plus_2si2 = Array{S}(undef, N, N)   # \dot{s}_j^2 + 2\dot{s}_i^2 inside {}
    sj2_plus_2si2_minus_4vivj = Array{S}(undef, N, N)  # \dot{s}_j^2 + 2\dot{s}_i^2 - 4<, > terms inside {}
    ϕs_and_vs = Array{S}(undef, N, N)       # -4\sum - \sum + \dot{s}_j^2 + 2\dot{s}_i^2  - 4<, > terms inside {}
    pn1t1_7 = Array{S}(undef, N, N)         # Everything inside the {} in the first term except for the term with accelerations (last)
    # Last term inside the {}
    pNX_t_X = Array{S}(undef, N, N)     # X-axis component
    pNY_t_Y = Array{S}(undef, N, N)     # Y-axis component
    pNZ_t_Z = Array{S}(undef, N, N)     # Z-axis component
    # Everything inside the {} in the first term
    pn1 = Array{S}(undef, N, N)
    # Full first term
    X_t_pn1 = Array{S}(undef, N, N)     # X-axis component
    Y_t_pn1 = Array{S}(undef, N, N)     # Y-axis component
    Z_t_pn1 = Array{S}(undef, N, N)     # Z-axis component

    # Temporary post-Newtonian accelerations
    pntempX = Array{S}(undef, N)        # X-axis component
    pntempY = Array{S}(undef, N)        # Y-axis component
    pntempZ = Array{S}(undef, N)        # Z-axis component
    # Full post-Newtonian accelerations
    postNewtonX = Array{S}(undef, N)    # X-axis component
    postNewtonY = Array{S}(undef, N)    # Y-axis component
    postNewtonZ = Array{S}(undef, N)    # Z-axis component

    #=
    Extended body accelerations
    See equation (28) in page 13 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract
    and equations (173) and (174) in page 33 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract
    =#

    # (J_n, C_{mn}, S_{mn}) acceleration auxiliaries

    # Auxiliaries to compute body-fixed frame coordinates
    X_bf_1 = Array{S}(undef, N_ext, N_ext)
    Y_bf_1 = Array{S}(undef, N_ext, N_ext)
    Z_bf_1 = Array{S}(undef, N_ext, N_ext)
    X_bf_2 = Array{S}(undef, N_ext, N_ext)
    Y_bf_2 = Array{S}(undef, N_ext, N_ext)
    Z_bf_2 = Array{S}(undef, N_ext, N_ext)
    X_bf_3 = Array{S}(undef, N_ext, N_ext)
    Y_bf_3 = Array{S}(undef, N_ext, N_ext)
    Z_bf_3 = Array{S}(undef, N_ext, N_ext)
    # Body-fixed frame coordinates
    X_bf = Array{S}(undef, N_ext, N_ext)
    Y_bf = Array{S}(undef, N_ext, N_ext)
    Z_bf = Array{S}(undef, N_ext, N_ext)

    # Extended body accelerations (without mass parameter) in the inertial frame
    F_JCS_x = Array{S}(undef, N_ext, N_ext)
    F_JCS_y = Array{S}(undef, N_ext, N_ext)
    F_JCS_z = Array{S}(undef, N_ext, N_ext)
    # Temporary arrays for the sum of full extended body accelerations
    temp_accX_j = Array{S}(undef, N_ext, N_ext)
    temp_accY_j = Array{S}(undef, N_ext, N_ext)
    temp_accZ_j = Array{S}(undef, N_ext, N_ext)
    temp_accX_i = Array{S}(undef, N_ext, N_ext)
    temp_accY_i = Array{S}(undef, N_ext, N_ext)
    temp_accZ_i = Array{S}(undef, N_ext, N_ext)

    # Trigonometric functions of latitude ϕ and longitude λ in the body-fixed coordinate system
    # See equations (165)-(168) in page 33 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract
    sin_ϕ = Array{S}(undef, N_ext, N_ext)
    cos_ϕ = Array{S}(undef, N_ext, N_ext)
    sin_λ = Array{S}(undef, N_ext, N_ext)
    cos_λ = Array{S}(undef, N_ext, N_ext)

    # Distances
    r_xy = Array{S}(undef, N_ext, N_ext)  # X-Y projection magnitude in body-fixed frame sqrt(x_b^2 + y_b^2)
    r_p4 = Array{S}(undef, N_ext, N_ext)  # r_{ij}^4
    # Legendre polynomials
    P_n = Array{S}(undef, N_ext, N_ext, maximum(n1SEM)+1)   # Vector of Legendre polynomials
    dP_n = Array{S}(undef, N_ext, N_ext, maximum(n1SEM)+1)  # Vector of d/d(sin ϕ) of Legendre polynomials

    # Temporary arrays for the sum of accelerations due to zonal harmonics J_n
    temp_fjξ = Array{S}(undef, N_ext, N_ext, maximum(n1SEM)+1)       # ξ-axis component
    temp_fjζ = Array{S}(undef, N_ext, N_ext, maximum(n1SEM)+1)       # ζ-axis component
    temp_rn = Array{S}(undef, N_ext, N_ext, maximum(n1SEM)+1)        # r_{ij}^{n+2}
    # Temporary arrays for the vector sum in equation (173) in page 33 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract
    temp_CS_ξ = Array{S}(undef, N_ext, N_ext, n1SEM[mo], n1SEM[mo])  # ξ-axis component
    temp_CS_η = Array{S}(undef, N_ext, N_ext, n1SEM[mo], n1SEM[mo])  # η-axis component
    temp_CS_ζ = Array{S}(undef, N_ext, N_ext, n1SEM[mo], n1SEM[mo])  # ζ-axis component
    # Accelerations due to lunar tesseral harmonics beyond C_{21} and S_{21}
    F_CS_ξ_36 = Array{S}(undef, N_ext, N_ext)  # ξ-axis component
    F_CS_η_36 = Array{S}(undef, N_ext, N_ext)  # η-axis component
    F_CS_ζ_36 = Array{S}(undef, N_ext, N_ext)  # ζ-axis component
    # Accelerations due to third zonal harmonic and beyond
    F_J_ξ_36 = Array{S}(undef, N_ext, N_ext)   # ξ-axis component
    F_J_ζ_36 = Array{S}(undef, N_ext, N_ext)   # ζ-axis component

    # Trigonometric functions of integer multiples the longitude λ in the body-fixed coordinate system
    # See equations (165)-(168) in page 33 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract
    sin_mλ = Array{S}(undef, N_ext, N_ext, n1SEM[mo])
    cos_mλ = Array{S}(undef, N_ext, N_ext, n1SEM[mo])
    # Lunar teseral harmonics C_{nm}/S_{nm} * trigonometric function of integer times the longitude λ
    Cnm_cosmλ = Array{S}(undef, N_ext, N_ext, n1SEM[mo]+1, n1SEM[mo]+1)
    Cnm_sinmλ = Array{S}(undef, N_ext, N_ext, n1SEM[mo]+1, n1SEM[mo]+1)
    Snm_cosmλ = Array{S}(undef, N_ext, N_ext, n1SEM[mo]+1, n1SEM[mo]+1)
    Snm_sinmλ = Array{S}(undef, N_ext, N_ext, n1SEM[mo]+1, n1SEM[mo]+1)

    # Associated Legendre functions
    secϕ_P_nm = Array{S}(undef, N_ext, N_ext, n1SEM[mo]+1, n1SEM[mo]+1)   # secϕ P_n^m
    P_nm = Array{S}(undef, N_ext, N_ext, n1SEM[mo]+1, n1SEM[mo]+1)        # Vector of associated Legendre functions
    cosϕ_dP_nm = Array{S}(undef, N_ext, N_ext, n1SEM[mo]+1, n1SEM[mo]+1)  # cosϕ d/d(sin ϕ)P_n^m
    # Accelerations due to second zonal harmonic
    F_J_ξ = Array{S}(undef, N_ext, N_ext)   # ξ-axis component
    F_J_η = Array{S}(undef, N_ext, N_ext)   # η-axis component
    F_J_ζ = Array{S}(undef, N_ext, N_ext)   # ζ-axis component
    # Accelerations due to lunar tesseral harmonics C_{21} and S_{21}
    F_CS_ξ = Array{S}(undef, N_ext, N_ext)  # ξ-axis component
    F_CS_η = Array{S}(undef, N_ext, N_ext)  # η-axis component
    F_CS_ζ = Array{S}(undef, N_ext, N_ext)  # ζ-axis component
    # Sum of the zonal and tesseral (only for the moon) accelerations without mass parameter
    # in body-fixed frame
    F_JCS_ξ = Array{S}(undef, N_ext, N_ext) # ξ-axis component
    F_JCS_η = Array{S}(undef, N_ext, N_ext) # η-axis component
    F_JCS_ζ = Array{S}(undef, N_ext, N_ext) # ζ-axis component

    # Rotation matrices

    # R matrix body-fixed -> "primed" (ξ, η, ζ) frame
    # See equation (161) in page 32 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract
    Rb2p = Array{S}(undef, N_ext, N_ext, 3, 3)
    # G matrix "space-fixed" -> "primed" (ξ, η, ζ) frame
    # See equation (163) in page 32 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract
    Gc2p = Array{S}(undef, N_ext, N_ext, 3, 3)

    # Full extended-body accelerations
    accX = Array{S}(undef, N_ext)
    accY = Array{S}(undef, N_ext)
    accZ = Array{S}(undef, N_ext)

    # Lunar torques
    # See equation (43) in page 18 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract
    # Vector of lunar torques
    N_MfigM_pmA_x = Array{S}(undef, N_ext)   # x-axis component
    N_MfigM_pmA_y = Array{S}(undef, N_ext)   # y-axis component
    N_MfigM_pmA_z = Array{S}(undef, N_ext)   # z-axis component
    # Temporary array for the sum of lunar torques
    temp_N_M_x = Array{S}(undef, N_ext)       # x-axis component
    temp_N_M_y = Array{S}(undef, N_ext)       # y-axis component
    temp_N_M_z = Array{S}(undef, N_ext)       # z-axis component
    # Total lunar torque
    N_MfigM = Array{S}(undef, 3)
    N_MfigM[1] = zero_q_1                     # x-axis component
    N_MfigM[2] = zero_q_1                     # y-axis component
    N_MfigM[3] = zero_q_1                     # z-axis component

    # Rotations to and from Earth, Sun and Moon pole-oriented frames
    local αs = deg2rad(α_p_sun*one_t)           # Sun's rotation pole right ascension (radians)
    local δs = deg2rad(δ_p_sun*one_t)           # Sun's rotation pole right ascension (radians)
    # Space-fixed -> body-fixed coordinate transformations
    RotM = Array{S}(undef, 3, 3, 5)
    local RotM[:,:,ea] = c2t_jpl_de430(dsj2k)   # Earth
    local RotM[:,:,su] = pole_rotation(αs, δs)  # Sun
    # Lunar mantle Euler angles
    ϕ_m = q[6N+1]
    θ_m = q[6N+2]
    ψ_m = q[6N+3]
    # Lunar mantle space-fixed -> body-fixed coodinate transformations
    # See equations (10)-(13) in page 9 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract
    RotM[1,1,mo] = (cos(ϕ_m)*cos(ψ_m)) - (cos(θ_m)*(sin(ϕ_m)*sin(ψ_m)))
    RotM[2,1,mo] = (-cos(θ_m)*(cos(ψ_m)*sin(ϕ_m))) - (cos(ϕ_m)*sin(ψ_m))
    RotM[3,1,mo] = sin(θ_m)*sin(ϕ_m)
    RotM[1,2,mo] = (cos(ψ_m)*sin(ϕ_m)) + (cos(θ_m)*(cos(ϕ_m)*sin(ψ_m)))
    RotM[2,2,mo] = (cos(θ_m)*(cos(ϕ_m)*cos(ψ_m))) - (sin(ϕ_m)*sin(ψ_m))
    RotM[3,2,mo] = (-cos(ϕ_m))*sin(θ_m)
    RotM[1,3,mo] = sin(θ_m)*sin(ψ_m)
    RotM[2,3,mo] = cos(ψ_m)*sin(θ_m)
    RotM[3,3,mo] = cos(θ_m)
    # Lunar mantle frame -> inertial frame -> Lunar core-equatorial frame coord transformation
    # See equation (16) in page 10 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract
    mantlef2coref = Array{S}(undef, 3, 3)
    # Lunar core Euler angle
    ϕ_c = q[6N+7]
    # mantlef2coref = R_z(ϕ_c)*[ R_z(ψ_m)*R_x(θ_m)*R_z(ϕ_m) ]^T
    mantlef2coref[1,1] = (( RotM[1,1,mo])*cos(ϕ_c)) + (RotM[1,2,mo]*sin(ϕ_c))
    mantlef2coref[2,1] = ((-RotM[1,1,mo])*sin(ϕ_c)) + (RotM[1,2,mo]*cos(ϕ_c))
    mantlef2coref[3,1] =  RotM[1,3,mo]
    mantlef2coref[1,2] = (( RotM[2,1,mo])*cos(ϕ_c)) + (RotM[2,2,mo]*sin(ϕ_c))
    mantlef2coref[2,2] = ((-RotM[2,1,mo])*sin(ϕ_c)) + (RotM[2,2,mo]*cos(ϕ_c))
    mantlef2coref[3,2] =  RotM[2,3,mo]
    mantlef2coref[1,3] = (( RotM[3,1,mo])*cos(ϕ_c)) + (RotM[3,2,mo]*sin(ϕ_c))
    mantlef2coref[2,3] = ((-RotM[3,1,mo])*sin(ϕ_c)) + (RotM[3,2,mo]*cos(ϕ_c))
    mantlef2coref[3,3] =  RotM[3,3,mo]
    # Core angular velocity in core-equatorial frame
    ω_c_CE_1 = (mantlef2coref[1,1]*q[6N+10]) + ((mantlef2coref[1,2]*q[6N+11]) + (mantlef2coref[1,3]*q[6N+12]))
    ω_c_CE_2 = (mantlef2coref[2,1]*q[6N+10]) + ((mantlef2coref[2,2]*q[6N+11]) + (mantlef2coref[2,3]*q[6N+12]))
    ω_c_CE_3 = (mantlef2coref[3,1]*q[6N+10]) + ((mantlef2coref[3,2]*q[6N+11]) + (mantlef2coref[3,3]*q[6N+12]))

    # Second zonal harmonic coefficient
    # See Table 10 in page 50 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract
    local J2E_t = (J2E + J2EDOT*(dsj2k/yr))*(RE_au^2)  # Earth (considering a linear change in time with rate J2EDOT)
    local J2S_t = JSEM[su,2]*one_t                     # Sun (static)
    # Vector of second zonal harmonic coefficients
    J2_t = Array{S}(undef, 5)
    J2_t[su] = J2S_t             # Earth
    J2_t[ea] = J2E_t             # Sun
    # Lunar torques: overall numerical factor in equation (44) in page 18 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract
    local N_MfigM_figE_factor = 7.5*μ[ea]*J2E_t

    #=
    Compute point-mass Newtonian accelerations, all bodies
    See equation (35) in page 7 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract
    =#
    for j in 1:N
        # Fill point-mass Newton accelerations  with zeros
        newtonX[j] = zero_q_1
        newtonY[j] = zero_q_1
        newtonZ[j] = zero_q_1
        newtonianNb_Potential[j] = zero_q_1
        # Fill first 3N elements of dq with velocities
        dq[3j-2] = q[3(N+j)-2]
        dq[3j-1] = q[3(N+j)-1]
        dq[3j  ] = q[3(N+j)  ]
    end
    # Fill extended-body accelerations with zeros
    for j in 1:N_ext
        accX[j] = zero_q_1
        accY[j] = zero_q_1
        accZ[j] = zero_q_1
    end

    for j in 1:N
        for i in 1:N
            # i == j && continue
            if i == j
                continue
            else
                # Difference in position \mathbf{r_i} - \mathbf{r_j}
                X[i,j] = q[3i-2]-q[3j-2]      # X-axis component
                Y[i,j] = q[3i-1]-q[3j-1]      # Y-axis component
                Z[i,j] = q[3i]-q[3j]          # Z-axis component

                # Difference in velocity \mathbf{v_i} - \mathbf{v_j}
                U[i,j] = dq[3i-2]-dq[3j-2]    # X-axis component
                V[i,j] = dq[3i-1]-dq[3j-1]    # Y-axis component
                W[i,j] = dq[3i  ]-dq[3j  ]    # Z-axis component

                # Weighted difference in velocity 4\mathbf{v_i} - 3\mathbf{v_j}
                _4U_m_3X[i,j] = (4dq[3j-2])-(3dq[3i-2]) # X-axis component
                _4V_m_3Y[i,j] = (4dq[3j-1])-(3dq[3i-1]) # Y-axis component
                _4W_m_3Z[i,j] = (4dq[3j  ])-(3dq[3i  ]) # Z-axis component

                # Dot product inside [] in the second term
                pn2x = X[i,j]*_4U_m_3X[i,j]
                pn2y = Y[i,j]*_4V_m_3Y[i,j]
                pn2z = Z[i,j]*_4W_m_3Z[i,j]

                # Product of velocity components
                UU[i,j] = dq[3i-2]*dq[3j-2]   # v_{ix}v_{jx}
                VV[i,j] = dq[3i-1]*dq[3j-1]   # v_{iy}v_{jy}
                WW[i,j] = dq[3i  ]*dq[3j  ]   # v_{iz}v_{jz}

                # Dot product of velocities \mathbf{v_i}\cdot\mathbf{v_j}
                vi_dot_vj[i,j] = ( UU[i,j]+VV[i,j] ) + WW[i,j]

                # Distances r_{ij} = ||\mathbf{r_i} - \mathbf{r_j}||
                r_p2[i,j] = ( (X[i,j]^2)+(Y[i,j]^2) ) + (Z[i,j]^2) # r_{ij}^2
                r_p1d2[i,j] = sqrt(r_p2[i,j])                      # r_{ij}
                r_p3d2[i,j] = r_p2[i,j]^1.5                        # r_{ij}^3
                r_p7d2[i,j] = r_p2[i,j]^3.5                        # r_{ij}^7

                # Newtonian coefficient, i.e., mass parameter / distance^3 -> \mu_i / r_{ij}^3
                newtonianCoeff[i,j] =  μ[i]/r_p3d2[i,j]

                # Second term without (\mathbf{v}_i - \mathbf{v}_j)
                pn2[i,j] = newtonianCoeff[i,j]*(( pn2x+pn2y ) + pn2z)

                # Newtonian coefficient * difference between two positions, i.e.,
                # \mu_i * (\mathbf{r_i} - \mathbf{r_j}) / r_{ij}^3
                newton_acc_X[i,j] = X[i,j]*newtonianCoeff[i,j]    # X-axis component
                newton_acc_Y[i,j] = Y[i,j]*newtonianCoeff[i,j]    # Y-axis component
                newton_acc_Z[i,j] = Z[i,j]*newtonianCoeff[i,j]    # Z-axis component

                # Newtonian potential of 1 body \mu_i / r_{ij}
                newtonian1b_Potential[i,j] = μ[i]/r_p1d2[i, j]
                # Third term without newtonian accelerations \mathbf{a}_i
                pn3[i,j] = 3.5newtonian1b_Potential[i,j]
                # Full second term
                U_t_pn2[i,j] = pn2[i,j]*U[i,j]   # X-axis component
                V_t_pn2[i,j] = pn2[i,j]*V[i,j]   # Y-axis component
                W_t_pn2[i,j] = pn2[i,j]*W[i,j]   # Z-axis component

                # Newtonian accelerations \mathbf{a}_{i} = \sum_{i\neq j} mu_i * (\mathbf{r_i} - \mathbf{r_j}) / r_{ij}^3
                temp_001 = newtonX[j] + (X[i,j]*newtonianCoeff[i,j])  # X-axis component
                newtonX[j] = temp_001
                temp_002 = newtonY[j] + (Y[i,j]*newtonianCoeff[i,j])  # Y-axis component
                newtonY[j] = temp_002
                temp_003 = newtonZ[j] + (Z[i,j]*newtonianCoeff[i,j])  # Z-axis component
                newtonZ[j] = temp_003
                # Newtonian potential of N bodies
                # \sum_{i\neq l} \frac{\mu_i}{r_{il}}
                temp_004 = newtonianNb_Potential[j] + newtonian1b_Potential[i, j]
                newtonianNb_Potential[j] = temp_004
            end # else (i != j)
        end #for, i
        # Velocity magnitude squared ||\mathbf{v}_i||^2
        v2[j] = ( (dq[3j-2]^2)+(dq[3j-1]^2) ) + (dq[3j]^2)
    end #for, j

    #=
    Extended body accelerations
    See equation (28) in page 13 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract
    and equations (173) and (174) in page 33 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract
    =#

    # 2nd-order lunar zonal (J_2) and tesseral (C_2, S_2) harmonics coefficients
    # times the equatorial radius of the moon squared R_M^2
    # See equation (30) in page 13 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract
    J2M_t = ( I_M_t[3,3] - ((I_M_t[1,1]+I_M_t[2,2])/2) )/(μ[mo]) # J_{2,M}*R_M^2
    C22M_t = ((I_M_t[2,2] - I_M_t[1,1])/(μ[mo]))/4               # C_{22,M}*R_M^2
    C21M_t = (-I_M_t[1,3])/(μ[mo])                               # C_{21,M}*R_M^2
    S21M_t = (-I_M_t[3,2])/(μ[mo])                               # S_{21,M}*R_M^2
    S22M_t = ((-I_M_t[2,1])/(μ[mo]))/2                           # S_{22,M}*R_M^2
    J2_t[mo] = J2M_t

    for j in 1:N_ext
        for i in 1:N_ext
            # i == j && continue
            if i == j
                continue
            else
                # J_n, C_{nm}, S_{nm} accelerations, if j-th body is flattened
                if UJ_interaction[i,j]
                    # Rotate from (X, Y, Z) inertial frame to (X_bf, Y_bf, Z_by) extended-body frame
                    X_bf_1[i,j] = X[i,j]*RotM[1,1,j]
                    X_bf_2[i,j] = Y[i,j]*RotM[1,2,j]
                    X_bf_3[i,j] = Z[i,j]*RotM[1,3,j]
                    Y_bf_1[i,j] = X[i,j]*RotM[2,1,j]
                    Y_bf_2[i,j] = Y[i,j]*RotM[2,2,j]
                    Y_bf_3[i,j] = Z[i,j]*RotM[2,3,j]
                    Z_bf_1[i,j] = X[i,j]*RotM[3,1,j]
                    Z_bf_2[i,j] = Y[i,j]*RotM[3,2,j]
                    Z_bf_3[i,j] = Z[i,j]*RotM[3,3,j]
                    X_bf[i,j] = (X_bf_1[i,j] + X_bf_2[i,j]) + (X_bf_3[i,j]) # x-coordinate in body-fixed frame
                    Y_bf[i,j] = (Y_bf_1[i,j] + Y_bf_2[i,j]) + (Y_bf_3[i,j]) # y-coordinate in body-fixed frame
                    Z_bf[i,j] = (Z_bf_1[i,j] + Z_bf_2[i,j]) + (Z_bf_3[i,j]) # z-coordinate in body-fixed frame

                    # Trigonometric functions of latitude ϕ and longitude λ in the body-fixed coordinate system
                    # See equations (165)-(168) in page 33 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract
                    sin_ϕ[i,j] = Z_bf[i,j]/r_p1d2[i,j]               # eq. (165)
                    r_xy[i,j] = sqrt( (X_bf[i,j]^2)+(Y_bf[i,j]^2) )  # X-Y projection magnitude in body-fixed frame sqrt(x_b^2 + y_b^2)
                    cos_ϕ[i,j] = r_xy[i,j]/r_p1d2[i,j]               # eq. (166)
                    sin_λ[i,j] = Y_bf[i,j]/r_xy[i,j]                 # eq. (167)
                    cos_λ[i,j] = X_bf[i,j]/r_xy[i,j]                 # eq. (168)

                    # Legendre polynomials

                    # See equations (176) and (177) in page 33 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract
                    P_n[i,j,1] = one_t       # Zeroth Legendre polynomial
                    P_n[i,j,2] = sin_ϕ[i,j]  # First Legendre polynomial
                    dP_n[i,j,1] = zero_q_1   # d/d(sin_ϕ) of zeroth Legendre polynomial
                    dP_n[i,j,2] = one_t      # d/d(sin_ϕ) of first Legendre polynomial

                    for n in 2:n1SEM[j] # min(3,n1SEM[j])
                        # Recursion relation for the n-th Legenre polynomial
                        # See equation (175) in page 33 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract
                        P_n[i,j,n+1] = ((P_n[i,j,n]*sin_ϕ[i,j])*fact1_jsem[n]) - (P_n[i,j,n-1]*fact2_jsem[n])
                        # Recursion relation for d/d(sin_ϕ) of the n-th Legendre polynomial
                        # See equation (178) in page 33 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract
                        dP_n[i,j,n+1] = (dP_n[i,j,n]*sin_ϕ[i,j]) + (P_n[i,j,n]*fact3_jsem[n])
                        # r_{ij}^{n+2}
                        temp_rn[i,j,n] = r_p1d2[i,j]^fact5_jsem[n]
                    end
                    r_p4[i,j] = r_p2[i,j]^2  # r_{ij}^4

                    # Compute accelerations due to zonal harmonics J_n

                    # Second zonal harmonic J_2
                    # See equation (28) in page 13 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract
                    # and equation (173) in page 33 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract
                    F_J_ξ[i,j] = ((P_n[i,j,3]*fact4_jsem[2])*J2_t[j])/r_p4[i,j]   # ξ-axis
                    F_J_ζ[i,j] = (((-dP_n[i,j,3])*cos_ϕ[i,j])*J2_t[j])/r_p4[i,j]  # ζ-axis
                    # Beyond third zonal harmonic J_3,...
                    F_J_ξ_36[i,j] = zero_q_1
                    F_J_ζ_36[i,j] = zero_q_1
                    for n in 3:n1SEM[j] # min(3,n1SEM[j])
                        # ξ-axis
                        temp_fjξ[i,j,n] = (((P_n[i,j,n+1]*fact4_jsem[n])*JSEM[j,n])/temp_rn[i,j,n]) + F_J_ξ_36[i,j]
                        # ζ-axis
                        temp_fjζ[i,j,n] = ((((-dP_n[i,j,n+1])*cos_ϕ[i,j])*JSEM[j,n])/temp_rn[i,j,n]) + F_J_ζ_36[i,j]
                        F_J_ξ_36[i,j] = temp_fjξ[i,j,n]
                        F_J_ζ_36[i,j] = temp_fjζ[i,j,n]
                    end

                    # Associate Legendre functions (only for the moon)
                    if j == mo
                        for m in 1:n1SEM[mo]
                            if m == 1
                                # In this case associate Legendre functions reduce to Legendre polynomials
                                # See equations (167) and (168) in page 33 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract
                                sin_mλ[i,j,1] = sin_λ[i,j] # Moyer (1971), eq. (167)
                                cos_mλ[i,j,1] = cos_λ[i,j] # Moyer (1971), eq. (168)
                                # sec( Associate Legendre polynomial with m = n = 1 )
                                # See equation (181) in page 33 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract
                                secϕ_P_nm[i,j,1,1] = one_t
                                # Associate Legendre polynomial with m = n = 1
                                P_nm[i,j,1,1] = cos_ϕ[i,j]
                                # cosϕP_1^1'
                                # See equation (183) in page 34 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract
                                # Note: the second term equation (183) vanishes when n = m
                                cosϕ_dP_nm[i,j,1,1] = sin_ϕ[i,j]*lnm3[1]
                            else
                                # Trigonometric identity sin(λ + (m - 1)λ) and cos(λ + (m - 1)λ)
                                sin_mλ[i,j,m] = (cos_mλ[i,j,m-1]*sin_mλ[i,j,1]) + (sin_mλ[i,j,m-1]*cos_mλ[i,j,1])
                                cos_mλ[i,j,m] = (cos_mλ[i,j,m-1]*cos_mλ[i,j,1]) - (sin_mλ[i,j,m-1]*sin_mλ[i,j,1])
                                # secϕ P_n^n
                                # See equation (180) in page 33 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract
                                secϕ_P_nm[i,j,m,m] = (secϕ_P_nm[i,j,m-1,m-1]*cos_ϕ[i,j])*lnm5[m]
                                # Associate Legendre polynomial with n = m
                                P_nm[i,j,m,m] = secϕ_P_nm[i,j,m,m]*cos_ϕ[i,j]
                                # cosϕP_m^m'
                                # See equation (183) in page 34 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract
                                # Note: the second term equation (183) vanishes when n = m
                                cosϕ_dP_nm[i,j,m,m] = (secϕ_P_nm[i,j,m,m]*sin_ϕ[i,j])*lnm3[m]
                            end
                            for n in m+1:n1SEM[mo]
                                # secϕ P_n^m
                                # See equation (182) in page 34 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract
                                if n == m+1
                                    secϕ_P_nm[i,j,n,m] = (secϕ_P_nm[i,j,n-1,m]*sin_ϕ[i,j])*lnm1[n,m]
                                else
                                    secϕ_P_nm[i,j,n,m] = ((secϕ_P_nm[i,j,n-1,m]*sin_ϕ[i,j])*lnm1[n,m]) + (secϕ_P_nm[i,j,n-2,m]*lnm2[n,m])
                                end
                                # Associate Legendre polynomial of degree n and order m
                                P_nm[i,j,n,m] = secϕ_P_nm[i,j,n,m]*cos_ϕ[i,j]
                                # secϕ P_n^m
                                # See equation (183) in page 34 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract
                                cosϕ_dP_nm[i,j,n,m] = ((secϕ_P_nm[i,j,n,m]*sin_ϕ[i,j])*lnm3[n]) + (secϕ_P_nm[i,j,n-1,m]*lnm4[n,m])
                            end
                        end

                        # Moon: Compute accelerations due to tesseral harmonics C_{nm}, S_{nm}
                        # See equation (28) in page 13 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract
                        # and equation (174) in page 33 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract

                        # Accelerations due to lunar tesseral harmonics C_{21} and S_{21}
                        F_CS_ξ[i,j] = (   (  (P_nm[i,j,2,1]*lnm6[2]     )*( (C21M_t*cos_mλ[i,j,1]) + (S21M_t*sin_mλ[i,j,1]) )  ) + (  (P_nm[i,j,2,2]*lnm6[2]     )*( (C22M_t*cos_mλ[i,j,2]) + (S22M_t*sin_mλ[i,j,2]) )  )   )/r_p4[i,j]
                        F_CS_η[i,j] = (   (  (secϕ_P_nm[i,j,2,1]*lnm7[1])*( (S21M_t*cos_mλ[i,j,1]) - (C21M_t*sin_mλ[i,j,1]) )  ) + (  (secϕ_P_nm[i,j,2,2]*lnm7[2])*( (S22M_t*cos_mλ[i,j,2]) - (C22M_t*sin_mλ[i,j,2]) )  )   )/r_p4[i,j]
                        F_CS_ζ[i,j] = (   (  (cosϕ_dP_nm[i,j,2,1]       )*( (C21M_t*cos_mλ[i,j,1]) + (S21M_t*sin_mλ[i,j,1]) )  ) + (  (cosϕ_dP_nm[i,j,2,2]       )*( (C22M_t*cos_mλ[i,j,2]) + (S22M_t*sin_mλ[i,j,2]) )  )   )/r_p4[i,j]
                        # Accelerations due to lunar tesseral harmonics beyond C_{21} and S_{21}
                        F_CS_ξ_36[i,j] = zero_q_1
                        F_CS_η_36[i,j] = zero_q_1
                        F_CS_ζ_36[i,j] = zero_q_1
                        for n in 3:n2M
                            for m in 1:n
                                # Lunar teseral harmonics C_{nm}/S_{nm} * trigonometric function of integer times the longitude λ
                                Cnm_cosmλ[i,j,n,m] = CM[n,m]*cos_mλ[i,j,m]
                                Cnm_sinmλ[i,j,n,m] = CM[n,m]*sin_mλ[i,j,m]
                                Snm_cosmλ[i,j,n,m] = SM[n,m]*cos_mλ[i,j,m]
                                Snm_sinmλ[i,j,n,m] = SM[n,m]*sin_mλ[i,j,m]
                                # Vector sum in equation (173)
                                temp_CS_ξ[i,j,n,m] = (   (  (P_nm[i,j,n,m]*lnm6[n]     )*( Cnm_cosmλ[i,j,n,m] + Snm_sinmλ[i,j,n,m] )  )/temp_rn[i,j,n]   ) + F_CS_ξ_36[i,j]
                                temp_CS_η[i,j,n,m] = (   (  (secϕ_P_nm[i,j,n,m]*lnm7[m])*( Snm_cosmλ[i,j,n,m] - Cnm_sinmλ[i,j,n,m] )  )/temp_rn[i,j,n]   ) + F_CS_η_36[i,j]
                                temp_CS_ζ[i,j,n,m] = (   (  (cosϕ_dP_nm[i,j,n,m]       )*( Cnm_cosmλ[i,j,n,m] + Snm_sinmλ[i,j,n,m] )  )/temp_rn[i,j,n]   ) + F_CS_ζ_36[i,j]
                                F_CS_ξ_36[i,j] = temp_CS_ξ[i,j,n,m]
                                F_CS_η_36[i,j] = temp_CS_η[i,j,n,m]
                                F_CS_ζ_36[i,j] = temp_CS_ζ[i,j,n,m]
                            end
                        end
                        # Sum the zonal and tesseral (only for the moon) accelerations without mass parameter
                        F_JCS_ξ[i,j] = (F_J_ξ[i,j] + F_J_ξ_36[i,j]) + (F_CS_ξ[i,j]+F_CS_ξ_36[i,j])
                        F_JCS_η[i,j] = (F_CS_η[i,j]+F_CS_η_36[i,j])
                        F_JCS_ζ[i,j] = (F_J_ζ[i,j] + F_J_ζ_36[i,j]) + (F_CS_ζ[i,j]+F_CS_ζ_36[i,j])
                    else
                        # Sum the zonal accelerations without mass parameter
                        F_JCS_ξ[i,j] = (F_J_ξ[i,j] + F_J_ξ_36[i,j])
                        F_JCS_η[i,j] = zero_q_1
                        F_JCS_ζ[i,j] = (F_J_ζ[i,j] + F_J_ζ_36[i,j])
                    end

                    # R matrix: body-fixed -> "primed" (ξ, η, ζ) system
                    # See equation (161) in page 32 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract
                    Rb2p[i,j,1,1] = cos_ϕ[i,j]*cos_λ[i,j]
                    Rb2p[i,j,2,1] = -sin_λ[i,j]
                    Rb2p[i,j,3,1] = -sin_ϕ[i,j]*cos_λ[i,j]
                    Rb2p[i,j,1,2] = cos_ϕ[i,j]*sin_λ[i,j]
                    Rb2p[i,j,2,2] = cos_λ[i,j]
                    Rb2p[i,j,3,2] = -sin_ϕ[i,j]*sin_λ[i,j]
                    Rb2p[i,j,1,3] = sin_ϕ[i,j]
                    Rb2p[i,j,2,3] = zero_q_1
                    Rb2p[i,j,3,3] = cos_ϕ[i,j]
                    # G matrix: space-fixed -> body-fixed -> "primed" (ξ, η, ζ) system
                    # G_{i,j} = \sum_k R_{i,k} RotM{k,j} or G = RotM * R
                    # See equation (163) in page 32 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract
                    Gc2p[i,j,1,1] = ((Rb2p[i,j,1,1]*RotM[1,1,j]) + (Rb2p[i,j,1,2]*RotM[2,1,j])) + (Rb2p[i,j,1,3]*RotM[3,1,j])
                    Gc2p[i,j,2,1] = ((Rb2p[i,j,2,1]*RotM[1,1,j]) + (Rb2p[i,j,2,2]*RotM[2,1,j])) + (Rb2p[i,j,2,3]*RotM[3,1,j])
                    Gc2p[i,j,3,1] = ((Rb2p[i,j,3,1]*RotM[1,1,j]) + (Rb2p[i,j,3,2]*RotM[2,1,j])) + (Rb2p[i,j,3,3]*RotM[3,1,j])
                    Gc2p[i,j,1,2] = ((Rb2p[i,j,1,1]*RotM[1,2,j]) + (Rb2p[i,j,1,2]*RotM[2,2,j])) + (Rb2p[i,j,1,3]*RotM[3,2,j])
                    Gc2p[i,j,2,2] = ((Rb2p[i,j,2,1]*RotM[1,2,j]) + (Rb2p[i,j,2,2]*RotM[2,2,j])) + (Rb2p[i,j,2,3]*RotM[3,2,j])
                    Gc2p[i,j,3,2] = ((Rb2p[i,j,3,1]*RotM[1,2,j]) + (Rb2p[i,j,3,2]*RotM[2,2,j])) + (Rb2p[i,j,3,3]*RotM[3,2,j])
                    Gc2p[i,j,1,3] = ((Rb2p[i,j,1,1]*RotM[1,3,j]) + (Rb2p[i,j,1,2]*RotM[2,3,j])) + (Rb2p[i,j,1,3]*RotM[3,3,j])
                    Gc2p[i,j,2,3] = ((Rb2p[i,j,2,1]*RotM[1,3,j]) + (Rb2p[i,j,2,2]*RotM[2,3,j])) + (Rb2p[i,j,2,3]*RotM[3,3,j])
                    Gc2p[i,j,3,3] = ((Rb2p[i,j,3,1]*RotM[1,3,j]) + (Rb2p[i,j,3,2]*RotM[2,3,j])) + (Rb2p[i,j,3,3]*RotM[3,3,j])
                    # Compute cartesian coordinates of acceleration due to body figure in inertial frame (without mass parameter)
                    # See equation (169) in page 33 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract
                    F_JCS_x[i,j] = ((F_JCS_ξ[i,j]*Gc2p[i,j,1,1]) + (F_JCS_η[i,j]*Gc2p[i,j,2,1])) + (F_JCS_ζ[i,j]*Gc2p[i,j,3,1])
                    F_JCS_y[i,j] = ((F_JCS_ξ[i,j]*Gc2p[i,j,1,2]) + (F_JCS_η[i,j]*Gc2p[i,j,2,2])) + (F_JCS_ζ[i,j]*Gc2p[i,j,3,2])
                    F_JCS_z[i,j] = ((F_JCS_ξ[i,j]*Gc2p[i,j,1,3]) + (F_JCS_η[i,j]*Gc2p[i,j,2,3])) + (F_JCS_ζ[i,j]*Gc2p[i,j,3,3])
                end #if UJ_interaction[i,j]
            end # else (i != j)
        end #for i in 1:N_ext
    end #for j in 1:N_ext

    for j in 1:N_ext
        for i in 1:N_ext
            # i == j && continue
            if i == j
                continue
            else
                if UJ_interaction[i,j]
                    # Extended body accelerations
                    # J_n, C_{nm}, S_{nm} accelerations, if j-th body is flattened

                    # Add result to total acceleration upon j-th body figure due to i-th point mass
                    temp_accX_j[i,j] = accX[j] - (μ[i]*F_JCS_x[i,j])
                    accX[j] = temp_accX_j[i,j]
                    temp_accY_j[i,j] = accY[j] - (μ[i]*F_JCS_y[i,j])
                    accY[j] = temp_accY_j[i,j]
                    temp_accZ_j[i,j] = accZ[j] - (μ[i]*F_JCS_z[i,j])
                    accZ[j] = temp_accZ_j[i,j]

                    # Reaction force on i-th body
                    temp_accX_i[i,j] = accX[i] + (μ[j]*F_JCS_x[i,j])
                    accX[i] = temp_accX_i[i,j]
                    temp_accY_i[i,j] = accY[i] + (μ[j]*F_JCS_y[i,j])
                    accY[i] = temp_accY_i[i,j]
                    temp_accZ_i[i,j] = accZ[i] + (μ[j]*F_JCS_z[i,j])
                    accZ[i] = temp_accZ_i[i,j]

                    # Lunar torques
                    if j == mo
                        # Compute torques acting upon the body-figure of the Moon due to external point masses
                        # See equation (43) in page 13 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract
                        N_MfigM_pmA_x[i] = μ[i]*( (Y[i,j]*F_JCS_z[i,j]) - (Z[i,j]*F_JCS_y[i,j]) )
                        N_MfigM_pmA_y[i] = μ[i]*( (Z[i,j]*F_JCS_x[i,j]) - (X[i,j]*F_JCS_z[i,j]) )
                        N_MfigM_pmA_z[i] = μ[i]*( (X[i,j]*F_JCS_y[i,j]) - (Y[i,j]*F_JCS_x[i,j]) )
                        # Expressions below have minus sign since N_MfigM_pmA_{x,y,z} have inverted signs in cross product
                        temp_N_M_x[i] = N_MfigM[1] - (N_MfigM_pmA_x[i]*μ[j])
                        N_MfigM[1] = temp_N_M_x[i]
                        temp_N_M_y[i] = N_MfigM[2] - (N_MfigM_pmA_y[i]*μ[j])
                        N_MfigM[2] = temp_N_M_y[i]
                        temp_N_M_z[i] = N_MfigM[3] - (N_MfigM_pmA_z[i]*μ[j])
                        N_MfigM[3] = temp_N_M_z[i]
                    end
                end
            end # else (i != j)
        end
    end

    #=
    Post-Newtonian corrections to gravitational acceleration
    Post-Newtonian iterative procedure setup and initialization
    See equation (35) in page 7 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract
    =#
    for j in 1:N
        for i in 1:N
            # i == j && continue
            if i == j
                continue
            else
                # 4*\sum term inside {}
                _4ϕj[i,j] = 4newtonianNb_Potential[j]
                # 4*\sum + \sum terms inside {}
                ϕi_plus_4ϕj[i,j] = newtonianNb_Potential[i] + _4ϕj[i,j]
                # 2 * ||\mathbf{v_i}||^2
                _2v2[i,j] = 2v2[i]
                # \dot{s}_j^2 + 2\dot{s}_i^2 inside {}
                sj2_plus_2si2[i,j] = v2[j] + _2v2[i,j]
                # \dot{s}_j^2 + 2\dot{s}_i^2 - 4<, > terms inside {}
                sj2_plus_2si2_minus_4vivj[i,j] = sj2_plus_2si2[i,j] - (4vi_dot_vj[i,j])
                # -4\sum - \sum + \dot{s}_j^2 + 2\dot{s}_i^2  - 4<, > terms inside {}
                ϕs_and_vs[i,j] = sj2_plus_2si2_minus_4vivj[i,j] - ϕi_plus_4ϕj[i,j]
                # (\mathbf{r}_i - \mathbf{r}_j)\cdot\mathbf{v_i}
                Xij_t_Ui = X[i,j]*dq[3i-2]
                Yij_t_Vi = Y[i,j]*dq[3i-1]
                Zij_t_Wi = Z[i,j]*dq[3i]
                Rij_dot_Vi = ( Xij_t_Ui+Yij_t_Vi ) + Zij_t_Wi
                # The expression below inside the (...)^2 should have a minus sign in front of the numerator,
                # but upon squaring it is eliminated, so at the end of the day, it is irrelevant ;)
                # (\mathbf{r}_i - \mathbf{r}_j)\cdot\mathbf{v_i} / r_{ij}
                pn1t7 = (Rij_dot_Vi^2)/r_p2[i,j]
                # Everything inside the {} except for the first and last terms
                pn1t2_7 = ϕs_and_vs[i,j] - (1.5pn1t7)
                # Everything inside the {} except for the last term
                pn1t1_7[i,j] = c_p2+pn1t2_7
            end # else (i != j)
        end
        # Temporary post-Newtonian accelerations
        pntempX[j] = zero_q_1   # X-axis component
        pntempY[j] = zero_q_1   # Y-axis component
        pntempZ[j] = zero_q_1   # Z-axis component
    end

    for j in 1:N
        for i in 1:N
            # i == j && continue
            if i == j
                continue
            else

                # First term of equation (35) in page 7 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract

                # Last term inside the {}
                pNX_t_X[i,j] = newtonX[i]*X[i,j]   # X-axis component
                pNY_t_Y[i,j] = newtonY[i]*Y[i,j]   # Y-axis component
                pNZ_t_Z[i,j] = newtonZ[i]*Z[i,j]   # Z-axis component
                # Everything inside the {} in the first term
                pn1[i,j] = (  pn1t1_7[i,j]  +  0.5*( (pNX_t_X[i,j]+pNY_t_Y[i,j]) + pNZ_t_Z[i,j] )  )
                # Full first term
                X_t_pn1[i,j] = newton_acc_X[i,j]*pn1[i,j]   # X-axis component
                Y_t_pn1[i,j] = newton_acc_Y[i,j]*pn1[i,j]   # Y-axis component
                Z_t_pn1[i,j] = newton_acc_Z[i,j]*pn1[i,j]   # Z-axis component

                # Full third term of equation (35) in page 7 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract
                pNX_t_pn3[i,j] = newtonX[i]*pn3[i,j]   # X-axis component
                pNY_t_pn3[i,j] = newtonY[i]*pn3[i,j]   # Y-axis component
                pNZ_t_pn3[i,j] = newtonZ[i]*pn3[i,j]   # Z-axis component

                # Temporary post-Newtonian accelerations
                termpnx = ( X_t_pn1[i,j] + (U_t_pn2[i,j]+pNX_t_pn3[i,j]) )   # X-axis component
                sumpnx = pntempX[j] + termpnx
                pntempX[j] = sumpnx
                termpny = ( Y_t_pn1[i,j] + (V_t_pn2[i,j]+pNY_t_pn3[i,j]) )   # Y-axis component
                sumpny = pntempY[j] + termpny
                pntempY[j] = sumpny
                termpnz = ( Z_t_pn1[i,j] + (W_t_pn2[i,j]+pNZ_t_pn3[i,j]) )   # Z-axis component
                sumpnz = pntempZ[j] + termpnz
                pntempZ[j] = sumpnz
            end # else (i != j)
        end
        # Post-Newtonian acelerations
        postNewtonX[j] = pntempX[j]*c_m2
        postNewtonY[j] = pntempY[j]*c_m2
        postNewtonZ[j] = pntempZ[j]*c_m2
    end

    # Fill accelerations (post-Newtonian and extended body accelerations)
    # postNewton -> post-Newtonian accelerations (all bodies)
    # accX/Y/Z -> extended body accelerations (only first N_ext bodies)
    for i in 1:N_ext
        dq[3(N+i)-2] = postNewtonX[i] + accX[i]
        dq[3(N+i)-1] = postNewtonY[i] + accY[i]
        dq[3(N+i)  ] = postNewtonZ[i] + accZ[i]
    end
    for i in N_ext+1:N
        dq[3(N+i)-2] = postNewtonX[i]
        dq[3(N+i)-1] = postNewtonY[i]
        dq[3(N+i)  ] = postNewtonZ[i]
    end

    #=
    Lunar physical librations
    See equations (33)-(35) in pages 15-16 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract
    =#

    # Lunar moment of intertia I times angular velocity ω: Iω
    Iω_x = (I_m_t[1,1]*q[6N+4]) + ((I_m_t[1,2]*q[6N+5]) + (I_m_t[1,3]*q[6N+6])) # x-axis component
    Iω_y = (I_m_t[2,1]*q[6N+4]) + ((I_m_t[2,2]*q[6N+5]) + (I_m_t[2,3]*q[6N+6])) # y-axis component
    Iω_z = (I_m_t[3,1]*q[6N+4]) + ((I_m_t[3,2]*q[6N+5]) + (I_m_t[3,3]*q[6N+6])) # z-axis component

    # Cross product of angular velocity and Iω: ω × (I*ω)
    ωxIω_x = (q[6N+5]*Iω_z) - (q[6N+6]*Iω_y)   # x-axis component
    ωxIω_y = (q[6N+6]*Iω_x) - (q[6N+4]*Iω_z)   # y-axis component
    ωxIω_z = (q[6N+4]*Iω_y) - (q[6N+5]*Iω_x)   # z-axis component

    # Time derivative of moment of inertia times angular velocity: (dI/dt)*ω
    dIω_x = (dI_m_t[1,1]*q[6N+4]) + ((dI_m_t[1,2]*q[6N+5]) + (dI_m_t[1,3]*q[6N+6])) # x-axis component
    dIω_y = (dI_m_t[2,1]*q[6N+4]) + ((dI_m_t[2,2]*q[6N+5]) + (dI_m_t[2,3]*q[6N+6])) # y-axis component
    dIω_z = (dI_m_t[3,1]*q[6N+4]) + ((dI_m_t[3,2]*q[6N+5]) + (dI_m_t[3,3]*q[6N+6])) # z-axis component

    # Moon -> Earth radial unit vector (inertial coordinates)
    er_EM_I_1 = X[ea,mo]/r_p1d2[ea,mo]
    er_EM_I_2 = Y[ea,mo]/r_p1d2[ea,mo]
    er_EM_I_3 = Z[ea,mo]/r_p1d2[ea,mo]

    # Earth pole unit vector (inertial coordinates)
    p_E_I_1 = RotM[3,1,ea]
    p_E_I_2 = RotM[3,2,ea]
    p_E_I_3 = RotM[3,3,ea]

    # Transform Moon -> Earth radial unit vector (inertial coordinates) er_EM_I_i and
    # Earth pole unit vector p_E_I_i to lunar mantle frame coordinates
    er_EM_1 = (RotM[1,1,mo]*er_EM_I_1) + ((RotM[1,2,mo]*er_EM_I_2) + (RotM[1,3,mo]*er_EM_I_3))
    er_EM_2 = (RotM[2,1,mo]*er_EM_I_1) + ((RotM[2,2,mo]*er_EM_I_2) + (RotM[2,3,mo]*er_EM_I_3))
    er_EM_3 = (RotM[3,1,mo]*er_EM_I_1) + ((RotM[3,2,mo]*er_EM_I_2) + (RotM[3,3,mo]*er_EM_I_3))
    p_E_1 = (RotM[1,1,mo]*p_E_I_1) + ((RotM[1,2,mo]*p_E_I_2) + (RotM[1,3,mo]*p_E_I_3))
    p_E_2 = (RotM[2,1,mo]*p_E_I_1) + ((RotM[2,2,mo]*p_E_I_2) + (RotM[2,3,mo]*p_E_I_3))
    p_E_3 = (RotM[3,1,mo]*p_E_I_1) + ((RotM[3,2,mo]*p_E_I_2) + (RotM[3,3,mo]*p_E_I_3))

    # Evaluate equation (44) https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract
    # in lunar mantle frame coords

    # I*e_r
    I_er_EM_1 = (I_m_t[1,1]*er_EM_1) + ((I_m_t[1,2]*er_EM_2) + (I_m_t[1,3]*er_EM_3))
    I_er_EM_2 = (I_m_t[2,1]*er_EM_1) + ((I_m_t[2,2]*er_EM_2) + (I_m_t[2,3]*er_EM_3))
    I_er_EM_3 = (I_m_t[3,1]*er_EM_1) + ((I_m_t[3,2]*er_EM_2) + (I_m_t[3,3]*er_EM_3))

    # I*p_E
    I_p_E_1 = (I_m_t[1,1]*p_E_1) + ((I_m_t[1,2]*p_E_2) + (I_m_t[1,3]*p_E_3))
    I_p_E_2 = (I_m_t[2,1]*p_E_1) + ((I_m_t[2,2]*p_E_2) + (I_m_t[2,3]*p_E_3))
    I_p_E_3 = (I_m_t[3,1]*p_E_1) + ((I_m_t[3,2]*p_E_2) + (I_m_t[3,3]*p_E_3))

    # e_r × (I*e_r)
    er_EM_cross_I_er_EM_1 = (er_EM_2*I_er_EM_3) - (er_EM_3*I_er_EM_2)
    er_EM_cross_I_er_EM_2 = (er_EM_3*I_er_EM_1) - (er_EM_1*I_er_EM_3)
    er_EM_cross_I_er_EM_3 = (er_EM_1*I_er_EM_2) - (er_EM_2*I_er_EM_1)

    # e_r × (I*p_E)
    er_EM_cross_I_p_E_1 = (er_EM_2*I_p_E_3) - (er_EM_3*I_p_E_2)
    er_EM_cross_I_p_E_2 = (er_EM_3*I_p_E_1) - (er_EM_1*I_p_E_3)
    er_EM_cross_I_p_E_3 = (er_EM_1*I_p_E_2) - (er_EM_2*I_p_E_1)

    # p_E × (I*e_r)
    p_E_cross_I_er_EM_1 = (p_E_2*I_er_EM_3) - (p_E_3*I_er_EM_2)
    p_E_cross_I_er_EM_2 = (p_E_3*I_er_EM_1) - (p_E_1*I_er_EM_3)
    p_E_cross_I_er_EM_3 = (p_E_1*I_er_EM_2) - (p_E_2*I_er_EM_1)

    # p_E × (I*p_E)
    p_E_cross_I_p_E_1 = (p_E_2*I_p_E_3) - (p_E_3*I_p_E_2)
    p_E_cross_I_p_E_2 = (p_E_3*I_p_E_1) - (p_E_1*I_p_E_3)
    p_E_cross_I_p_E_3 = (p_E_1*I_p_E_2) - (p_E_2*I_p_E_1)

    # Coefficients of first and second terms inside {} in eq. (44)
    one_minus_7sin2ϕEM = one_t - (7((sin_ϕ[ea,mo])^2))
    two_sinϕEM = 2sin_ϕ[ea,mo]

    # Overall numerical factor in eq. (44) / r_{EM}^5
    N_MfigM_figE_factor_div_rEMp5 = (N_MfigM_figE_factor/(r_p1d2[mo,ea]^5))
    # Evaluation of eq. (44)
    N_MfigM_figE_1 = N_MfigM_figE_factor_div_rEMp5*( (one_minus_7sin2ϕEM*er_EM_cross_I_er_EM_1) + (two_sinϕEM*(er_EM_cross_I_p_E_1+p_E_cross_I_er_EM_1)) - (0.4p_E_cross_I_p_E_1))
    N_MfigM_figE_2 = N_MfigM_figE_factor_div_rEMp5*( (one_minus_7sin2ϕEM*er_EM_cross_I_er_EM_2) + (two_sinϕEM*(er_EM_cross_I_p_E_2+p_E_cross_I_er_EM_2)) - (0.4p_E_cross_I_p_E_2))
    N_MfigM_figE_3 = N_MfigM_figE_factor_div_rEMp5*( (one_minus_7sin2ϕEM*er_EM_cross_I_er_EM_3) + (two_sinϕEM*(er_EM_cross_I_p_E_3+p_E_cross_I_er_EM_3)) - (0.4p_E_cross_I_p_E_3))

    # RHS of equation (34) in page 16 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract

    # Torques acting upon lunar body-figure due to external point masses: transform coordinates from inertial frame to lunar mantle frame
    N_1_LMF = (RotM[1,1,mo]*N_MfigM[1]) + ((RotM[1,2,mo]*N_MfigM[2]) + (RotM[1,3,mo]*N_MfigM[3]))
    N_2_LMF = (RotM[2,1,mo]*N_MfigM[1]) + ((RotM[2,2,mo]*N_MfigM[2]) + (RotM[2,3,mo]*N_MfigM[3]))
    N_3_LMF = (RotM[3,1,mo]*N_MfigM[1]) + ((RotM[3,2,mo]*N_MfigM[2]) + (RotM[3,3,mo]*N_MfigM[3]))

    # Torque on the mantle due to the interaction between core and mantle (evaluated in mantle frame)
    # See equation (45) in page 18 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract
    N_cmb_1 = (k_ν*(q[6N+10]-q[6N+4])) - (C_c_m_A_c*(q[6N+12]*q[6N+11]))
    N_cmb_2 = (k_ν*(q[6N+11]-q[6N+5])) + (C_c_m_A_c*(q[6N+12]*q[6N+10]))
    N_cmb_3 = (k_ν*(q[6N+12]-q[6N+6]))

    # I*(dω/dt); i.e., I times RHS of equation (34) in page 16 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract
    I_dω_1 = ((N_1_LMF + N_MfigM_figE_1) + N_cmb_1) - (dIω_x + ωxIω_x)
    I_dω_2 = ((N_2_LMF + N_MfigM_figE_2) + N_cmb_2) - (dIω_y + ωxIω_y)
    I_dω_3 = ((N_3_LMF + N_MfigM_figE_3) + N_cmb_3) - (dIω_z + ωxIω_z)

    # RHS of equation (35) in page 16 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract

    # I_c * ω_c
    Ic_ωc_1 = I_c_t[1,1]*q[6N+10] # + ((I_c_t[1,2]*q[6N+11]) + (I_c_t[1,3]*q[6N+12]))
    Ic_ωc_2 = I_c_t[2,2]*q[6N+11] # + ((I_c_t[2,1]*q[6N+10]) + (I_c_t[2,3]*q[6N+12]))
    Ic_ωc_3 = I_c_t[3,3]*q[6N+12] # + ((I_c_t[3,1]*q[6N+10]) + (I_c_t[3,2]*q[6N+11]))

    # - ω_m × (I_c * ω_c)
    m_ωm_x_Icωc_1 = (q[6N+6]*Ic_ωc_2) - (q[6N+5]*Ic_ωc_3)
    m_ωm_x_Icωc_2 = (q[6N+4]*Ic_ωc_3) - (q[6N+6]*Ic_ωc_1)
    m_ωm_x_Icωc_3 = (q[6N+5]*Ic_ωc_1) - (q[6N+4]*Ic_ωc_2)

    # I_c*(dω_c/dt); i.e., I_c times RHS of of equation (35) in page 16 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract
    Ic_dωc_1 = m_ωm_x_Icωc_1 - N_cmb_1
    Ic_dωc_2 = m_ωm_x_Icωc_2 - N_cmb_2
    Ic_dωc_3 = m_ωm_x_Icωc_3 - N_cmb_3

    # Lunar mantle physical librations

    # Euler angles
    # See equation (14) in page 9 of  https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract
    dq[6N+1] = ((q[6N+4]*sin(q[6N+3])) + (q[6N+5]*cos(q[6N+3])) )/sin(q[6N+2])
    dq[6N+2] = (q[6N+4]*cos(q[6N+3])) - (q[6N+5]*sin(q[6N+3]))
    dq[6N+3] = q[6N+6] - (dq[6N+1]*cos(q[6N+2]))
    # Angular velocitiy
    # See equation (34) in page 16 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract
    dq[6N+4] = (inv_I_m_t[1,1]*I_dω_1) + ( (inv_I_m_t[1,2]*I_dω_2) + (inv_I_m_t[1,3]*I_dω_3) )
    dq[6N+5] = (inv_I_m_t[2,1]*I_dω_1) + ( (inv_I_m_t[2,2]*I_dω_2) + (inv_I_m_t[2,3]*I_dω_3) )
    dq[6N+6] = (inv_I_m_t[3,1]*I_dω_1) + ( (inv_I_m_t[3,2]*I_dω_2) + (inv_I_m_t[3,3]*I_dω_3) )

    # Lunar core physical librations

    # Euler angles
    # See equation (15) in page 9 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract
    # (core angular velocity components ω_c_CE_i represent lunar core-equator frame coordinates)
    dq[6N+9] = -(ω_c_CE_2/sin(q[6N+8])) ### evaluated first, since it's used below
    dq[6N+7] = ω_c_CE_3-(dq[6N+9]*cos(q[6N+8]))
    dq[6N+8] = ω_c_CE_1
    # Angular velocity of the core expressed in the mantle frame
    # See equation (35) in page 16 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract
    dq[6N+10] = inv_I_c_t[1,1]*Ic_dωc_1 # + ( (inv_I_c_t[1,2]*Ic_dωc_2) + (inv_I_c_t[1,3]*Ic_dωc_3) )
    dq[6N+11] = inv_I_c_t[2,2]*Ic_dωc_2 # + ( (inv_I_c_t[2,1]*Ic_dωc_1) + (inv_I_c_t[2,3]*Ic_dωc_3) )
    dq[6N+12] = inv_I_c_t[3,3]*Ic_dωc_3 # + ( (inv_I_c_t[3,1]*Ic_dωc_1) + (inv_I_c_t[3,2]*Ic_dωc_2) )

    # TT-TDB
    # TODO: implement TT-TDB integration
    dq[6N+13] = zero_q_1

    nothing
end

@doc raw"""
    NBP_pN_A_J23E_J23M_J2S_threads!(dq, q, params, t)

Threaded version of `NBP_pN_A_J23E_J23M_J2S!`.

See also [`NBP_pN_A_J23E_J23M_J2S!`](@ref).
""" NBP_pN_A_J23E_J23M_J2S_threads!
function NBP_pN_A_J23E_J23M_J2S_threads!(dq, q, params, t)
    # N: number of bodies
    # jd0: initial Julian date
    local N, jd0 = params
    local S = eltype(q)   # Type of positions/velocities components

    local zero_q_1 = zero(q[1])                  # Zero of type S
    local one_t = one(t)                         # One of the same type as time t
    local dsj2k = t+(jd0-J2000)                  # Days since J2000.0 (TDB)
    # Matrix elements of lunar moment of inertia at time t-τ_M (without tidal distortion)
    # See equations (36) to (41) in pages 16-17 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract
    # ITM(q_del_τ_M, eulang_del_τ_M)
    local I_m_t = (ITM_und-I_c).*one_t           # Undistorted moment of inertia of the mantle, see equation (40)
    local dI_m_t = ordpres_differentiate.(I_m_t) # Time-derivative of lunar mantle I at time t-τ_M
    local inv_I_m_t = inv(I_m_t)                 # Inverse of lunar mantle I matrix at time t-τ_M
    local I_c_t = I_c.*one_t                     # Lunar core I matrix, see equation (39)
    local inv_I_c_t = inv(I_c_t)                 # Inverse of lunar core I matrix
    local I_M_t = I_m_t+I_c_t                    # Total I matrix (mantle + core)

    #=
    Point-mass accelerations
    See equation (35) in page 7 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract
    =#

    # Note: All the following arrays are declared here in order to help @taylorize work

    # Difference between two positions (\mathbf{r}_i - \mathbf{r}_j)
    X = Array{S}(undef, N, N)         # X-axis component
    Y = Array{S}(undef, N, N)         # Y-axis component
    Z = Array{S}(undef, N, N)         # Z-axis component

    # Distance between two positions r_{ij} = ||\mathbf{r}_i - \mathbf{r}_j||
    r_p2 = Array{S}(undef, N, N)      # r_{ij}^2
    r_p1d2 = Array{S}(undef, N, N)    # sqrt(r_p2) <-> r_{ij}
    r_p3d2 = Array{S}(undef, N, N)    # r_p2^1.5 <-> r_{ij}^3
    r_p7d2 = Array{S}(undef, N, N)    # r_p2^3.5 <-> r_{ij}^7

    # Newtonian accelerations \mathbf{a}_{i} = \sum_{i\neq j} mu_i * (\mathbf{r_i} - \mathbf{r_j}) / r_{ij}^3
    newtonX = Array{S}(undef, N)      # X-axis component
    newtonY = Array{S}(undef, N)      # Y-axis component
    newtonZ = Array{S}(undef, N)      # Z-axis component
    # Newtonian coefficient, i.e., mass parameter / distance^3 -> \mu_i / r_{ij}^3
    newtonianCoeff = Array{S}(undef, N, N)

    # Post-Newtonian stuff

    # Difference between two velocities (\mathbf{v}_i - \mathbf{v}_j)
    U = Array{S}(undef, N, N)         # X-axis component
    V = Array{S}(undef, N, N)         # Y-axis component
    W = Array{S}(undef, N, N)         # Z-axis component

    # Weighted difference between two velocities (4\mathbf{v}_i - 3\mathbf{v}_j)
    _4U_m_3X = Array{S}(undef, N, N)  # X-axis component
    _4V_m_3Y = Array{S}(undef, N, N)  # Y-axis component
    _4W_m_3Z = Array{S}(undef, N, N)  # Z-axis component

    # Product of velocity components
    UU = Array{S}(undef, N, N)        # v_{ix}v_{jx}
    VV = Array{S}(undef, N, N)        # v_{iy}v_{jy}
    WW = Array{S}(undef, N, N)        # v_{iz}v_{jz}

    # Newtonian potential of 1 body \mu_i / r_{ij}
    newtonian1b_Potential = Array{S}(undef, N, N)
    # Newtonian potential of N bodies
    # \sum_{i\neq l} \frac{\mu_i}{r_{il}}
    newtonianNb_Potential = Array{S}(undef, N)

    # Newtonian coefficient * difference between two positions, i.e.,
    # \mu_i * (\mathbf{r_i} - \mathbf{r_j}) / r_{ij}^3
    newton_acc_X = Array{S}(undef, N, N)
    newton_acc_Y = Array{S}(undef, N, N)
    newton_acc_Z = Array{S}(undef, N, N)

    # Combinations of velocities
    v2 = Array{S}(undef, N)                # Velocity magnitude squared ||\mathbf{v}_i||^2
    _2v2 = Array{S}(undef, N, N)           # 2 * ||\mathbf{v_i}||^2
    vi_dot_vj = Array{S}(undef, N, N)      # Dot product of two velocities \mathbf{v}_i\cdot\mathbf{v}_j

    # Second term of equation (35) in page 7 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract

    # Second term without (\mathbf{v}_i - \mathbf{v}_j)
    pn2 = Array{S}(undef, N, N)            # \mu_i * [(\mathbf{r_i} - \mathbf{r_j})\cdot(4\mathbf{v_i} - 3\mathbf{v_j})]
    # Full second term
    U_t_pn2 = Array{S}(undef, N, N)        # X-axis component
    V_t_pn2 = Array{S}(undef, N, N)        # Y-axis component
    W_t_pn2 = Array{S}(undef, N, N)        # Z-axis component

    # Third term of equation (35) in page 7 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract

    # Third term without newtonian accelerations \mathbf{a}_i
    pn3 = Array{S}(undef, N, N)
    # Full third term of equation (35)
    pNX_t_pn3 = Array{S}(undef, N, N)      # X-axis component
    pNY_t_pn3 = Array{S}(undef, N, N)      # Y-axis component
    pNZ_t_pn3 = Array{S}(undef, N, N)      # Z-axis component

    # First term of equation (35) in page 7 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract

    _4ϕj = Array{S}(undef, N, N)            # 4*\sum term inside {}
    ϕi_plus_4ϕj = Array{S}(undef, N, N)     # 4*\sum + \sum terms inside {}
    sj2_plus_2si2 = Array{S}(undef, N, N)   # \dot{s}_j^2 + 2\dot{s}_i^2 inside {}
    sj2_plus_2si2_minus_4vivj = Array{S}(undef, N, N)  # \dot{s}_j^2 + 2\dot{s}_i^2 - 4<, > terms inside {}
    ϕs_and_vs = Array{S}(undef, N, N)       # -4\sum - \sum + \dot{s}_j^2 + 2\dot{s}_i^2  - 4<, > terms inside {}
    pn1t1_7 = Array{S}(undef, N, N)         # Everything inside the {} in the first term except for the term with accelerations (last)
    # Last term inside the {}
    pNX_t_X = Array{S}(undef, N, N)     # X-axis component
    pNY_t_Y = Array{S}(undef, N, N)     # Y-axis component
    pNZ_t_Z = Array{S}(undef, N, N)     # Z-axis component
    # Everything inside the {} in the first term
    pn1 = Array{S}(undef, N, N)
    # Full first term
    X_t_pn1 = Array{S}(undef, N, N)    # X-axis component
    Y_t_pn1 = Array{S}(undef, N, N)    # Y-axis component
    Z_t_pn1 = Array{S}(undef, N, N)    # Z-axis component

    # Temporary post-Newtonian accelerations
    pntempX = Array{S}(undef, N)        # X-axis component
    pntempY = Array{S}(undef, N)        # Y-axis component
    pntempZ = Array{S}(undef, N)        # Z-axis component
    # Full post-Newtonian accelerations
    postNewtonX = Array{S}(undef, N)    # X-axis component
    postNewtonY = Array{S}(undef, N)    # Y-axis component
    postNewtonZ = Array{S}(undef, N)    # Z-axis component

    #=
    Extended body accelerations
    See equation (28) in page 13 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract
    and equations (173) and (174) in page 33 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract
    =#

    # (J_n, C_{nm}, S_{nm}) acceleration auxiliaries

    # Auxiliaries to compute body-fixed frame coordinates
    X_bf_1 = Array{S}(undef, N_ext, N_ext)
    Y_bf_1 = Array{S}(undef, N_ext, N_ext)
    Z_bf_1 = Array{S}(undef, N_ext, N_ext)
    X_bf_2 = Array{S}(undef, N_ext, N_ext)
    Y_bf_2 = Array{S}(undef, N_ext, N_ext)
    Z_bf_2 = Array{S}(undef, N_ext, N_ext)
    X_bf_3 = Array{S}(undef, N_ext, N_ext)
    Y_bf_3 = Array{S}(undef, N_ext, N_ext)
    Z_bf_3 = Array{S}(undef, N_ext, N_ext)
    # Body-fixed frame coordinates
    X_bf = Array{S}(undef, N_ext, N_ext)
    Y_bf = Array{S}(undef, N_ext, N_ext)
    Z_bf = Array{S}(undef, N_ext, N_ext)

    # Extended body accelerations (without mass parameter) in the inertial frame
    F_JCS_x = Array{S}(undef, N_ext, N_ext)
    F_JCS_y = Array{S}(undef, N_ext, N_ext)
    F_JCS_z = Array{S}(undef, N_ext, N_ext)
    # Temporary arrays for the sum of full extended body accelerations
    temp_accX_j = Array{S}(undef, N_ext, N_ext)
    temp_accY_j = Array{S}(undef, N_ext, N_ext)
    temp_accZ_j = Array{S}(undef, N_ext, N_ext)
    temp_accX_i = Array{S}(undef, N_ext, N_ext)
    temp_accY_i = Array{S}(undef, N_ext, N_ext)
    temp_accZ_i = Array{S}(undef, N_ext, N_ext)

    # Trigonometric functions of latitude ϕ and longitude λ in the body-fixed coordinate system
    # See equations (165)-(168) in page 33 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract
    sin_ϕ = Array{S}(undef, N_ext, N_ext)
    cos_ϕ = Array{S}(undef, N_ext, N_ext)
    sin_λ = Array{S}(undef, N_ext, N_ext)
    cos_λ = Array{S}(undef, N_ext, N_ext)

    # Distances
    r_xy = Array{S}(undef, N_ext, N_ext)  # X-Y projection magnitude in body-fixed frame sqrt(x_b^2 + y_b^2)
    r_p4 = Array{S}(undef, N_ext, N_ext)  # r_{ij}^4
    # Legendre polynomials
    P_n = Array{S}(undef, N_ext, N_ext, maximum(n1SEM)+1)   # Vector of Legendre polynomials
    dP_n = Array{S}(undef, N_ext, N_ext, maximum(n1SEM)+1)  # Vector of d/d(sin ϕ) of Legendre polynomials

    # Temporary arrays for the sum of accelerations due to zonal harmonics J_n
    temp_fjξ = Array{S}(undef, N_ext, N_ext, maximum(n1SEM)+1)       # ξ-axis component
    temp_fjζ = Array{S}(undef, N_ext, N_ext, maximum(n1SEM)+1)       # ζ-axis component
    temp_rn = Array{S}(undef, N_ext, N_ext, maximum(n1SEM)+1)        # r_{ij}^{n+2}
    # Temporary arrays for the vector sum in equation (173) in page 33 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract
    temp_CS_ξ = Array{S}(undef, N_ext, N_ext, n1SEM[mo], n1SEM[mo])  # ξ-axis component
    temp_CS_η = Array{S}(undef, N_ext, N_ext, n1SEM[mo], n1SEM[mo])  # η-axis component
    temp_CS_ζ = Array{S}(undef, N_ext, N_ext, n1SEM[mo], n1SEM[mo])  # ζ-axis component
    # Accelerations due to lunar tesseral harmonics beyond C_{21} and S_{21}
    F_CS_ξ_36 = Array{S}(undef, N_ext, N_ext)  # ξ-axis component
    F_CS_η_36 = Array{S}(undef, N_ext, N_ext)  # η-axis component
    F_CS_ζ_36 = Array{S}(undef, N_ext, N_ext)  # ζ-axis component
    # Accelerations due to third zonal harmonic and beyond
    F_J_ξ_36 = Array{S}(undef, N_ext, N_ext)   # ξ-axis component
    F_J_ζ_36 = Array{S}(undef, N_ext, N_ext)   # ζ-axis component

    # Trigonometric functions of integer multiples the longitude λ in the body-fixed coordinate system
    # See equations (165)-(168) in page 33 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract
    sin_mλ = Array{S}(undef, N_ext, N_ext, n1SEM[mo])
    cos_mλ = Array{S}(undef, N_ext, N_ext, n1SEM[mo])
    # Lunar teseral harmonics C_{nm}/S_{nm} * trigonometric function of integer times the longitude λ
    Cnm_cosmλ = Array{S}(undef, N_ext, N_ext, n1SEM[mo]+1, n1SEM[mo]+1)
    Cnm_sinmλ = Array{S}(undef, N_ext, N_ext, n1SEM[mo]+1, n1SEM[mo]+1)
    Snm_cosmλ = Array{S}(undef, N_ext, N_ext, n1SEM[mo]+1, n1SEM[mo]+1)
    Snm_sinmλ = Array{S}(undef, N_ext, N_ext, n1SEM[mo]+1, n1SEM[mo]+1)

    # Associated Legendre functions
    secϕ_P_nm = Array{S}(undef, N_ext, N_ext, n1SEM[mo]+1, n1SEM[mo]+1)   # secϕ P_n^m
    P_nm = Array{S}(undef, N_ext, N_ext, n1SEM[mo]+1, n1SEM[mo]+1)        # Vector of associated Legendre functions
    cosϕ_dP_nm = Array{S}(undef, N_ext, N_ext, n1SEM[mo]+1, n1SEM[mo]+1)  # cosϕ d/d(sin ϕ)P_n^m# Accelerations due to second zonal harmonic
    # Accelerations due to second zonal harmonic
    F_J_ξ = Array{S}(undef, N_ext, N_ext)   # ξ-axis component
    F_J_η = Array{S}(undef, N_ext, N_ext)   # η-axis component
    F_J_ζ = Array{S}(undef, N_ext, N_ext)   # ζ-axis component
    # Accelerations due to lunar tesseral harmonics C_{21} and S_{21}
    F_CS_ξ = Array{S}(undef, N_ext, N_ext)  # ξ-axis component
    F_CS_η = Array{S}(undef, N_ext, N_ext)  # η-axis component
    F_CS_ζ = Array{S}(undef, N_ext, N_ext)  # ζ-axis component
    # Sum of the zonal and tesseral (only for the moon) accelerations without mass parameter
    # in body-fixed frame
    F_JCS_ξ = Array{S}(undef, N_ext, N_ext) # ξ-axis component
    F_JCS_η = Array{S}(undef, N_ext, N_ext) # η-axis component
    F_JCS_ζ = Array{S}(undef, N_ext, N_ext) # ζ-axis component

    # Rotation matrices

    # R matrix body-fixed -> "primed" (ξ, η, ζ) frame
    # See equation (161) in page 32 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract
    Rb2p = Array{S}(undef, N_ext, N_ext, 3, 3)
    # G matrix "space-fixed" -> "primed" (ξ, η, ζ) frame
    # See equation (163) in page 32 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract
    Gc2p = Array{S}(undef, N_ext, N_ext, 3, 3)

    # Full extended-body accelerations
    accX = Array{S}(undef, N_ext)
    accY = Array{S}(undef, N_ext)
    accZ = Array{S}(undef, N_ext)

    # Lunar torques
    # See equation (43) in page 18 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract
    # Vector of lunar torques
    N_MfigM_pmA_x = Array{S}(undef, N_ext)   # x-axis component
    N_MfigM_pmA_y = Array{S}(undef, N_ext)   # y-axis component
    N_MfigM_pmA_z = Array{S}(undef, N_ext)   # z-axis component
    # Temporary array for the sum of lunar torques
    temp_N_M_x = Array{S}(undef, N_ext)       # x-axis component
    temp_N_M_y = Array{S}(undef, N_ext)       # y-axis component
    temp_N_M_z = Array{S}(undef, N_ext)       # z-axis component
    # Total lunar torque
    N_MfigM = Array{S}(undef, 3)
    N_MfigM[1] = zero_q_1                     # x-axis component
    N_MfigM[2] = zero_q_1                     # y-axis component
    N_MfigM[3] = zero_q_1                     # z-axis component

    # Rotations to and from Earth, Sun and Moon pole-oriented frames
    local αs = deg2rad(α_p_sun*one_t)           # Sun's rotation pole right ascension (radians)
    local δs = deg2rad(δ_p_sun*one_t)           # Sun's rotation pole right ascension (radians)
    # Space-fixed -> body-fixed coordinate transformations
    RotM = Array{S}(undef, 3, 3, 5)
    local RotM[:,:,ea] = c2t_jpl_de430(dsj2k)   # Earth
    local RotM[:,:,su] = pole_rotation(αs, δs)  # Sun
    # Lunar mantle Euler angles
    ϕ_m = q[6N+1]
    θ_m = q[6N+2]
    ψ_m = q[6N+3]
    # Lunar mantle space-fixed -> body-fixed coodinate transformations
    # See equations (10)-(13) in page 9 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract
    RotM[1,1,mo] = (cos(ϕ_m)*cos(ψ_m)) - (cos(θ_m)*(sin(ϕ_m)*sin(ψ_m)))
    RotM[2,1,mo] = (-cos(θ_m)*(cos(ψ_m)*sin(ϕ_m))) - (cos(ϕ_m)*sin(ψ_m))
    RotM[3,1,mo] = sin(θ_m)*sin(ϕ_m)
    RotM[1,2,mo] = (cos(ψ_m)*sin(ϕ_m)) + (cos(θ_m)*(cos(ϕ_m)*sin(ψ_m)))
    RotM[2,2,mo] = (cos(θ_m)*(cos(ϕ_m)*cos(ψ_m))) - (sin(ϕ_m)*sin(ψ_m))
    RotM[3,2,mo] = (-cos(ϕ_m))*sin(θ_m)
    RotM[1,3,mo] = sin(θ_m)*sin(ψ_m)
    RotM[2,3,mo] = cos(ψ_m)*sin(θ_m)
    RotM[3,3,mo] = cos(θ_m)
    # Lunar mantle frame -> inertial frame -> Lunar core-equatorial frame coord transformation
    # See equation (16) in page 10 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract
    mantlef2coref = Array{S}(undef, 3, 3) # lunar mantle frame -> inertial frame -> lunar core-equatorial frame coord transformation
    # Lunar core Euler angle
    ϕ_c = q[6N+7]
    # mantlef2coref = R_z(ϕ_c)*[ R_z(ψ_m)*R_x(θ_m)*R_z(ϕ_m) ]^T
    mantlef2coref[1,1] = (( RotM[1,1,mo])*cos(ϕ_c)) + (RotM[1,2,mo]*sin(ϕ_c))
    mantlef2coref[2,1] = ((-RotM[1,1,mo])*sin(ϕ_c)) + (RotM[1,2,mo]*cos(ϕ_c))
    mantlef2coref[3,1] =  RotM[1,3,mo]
    mantlef2coref[1,2] = (( RotM[2,1,mo])*cos(ϕ_c)) + (RotM[2,2,mo]*sin(ϕ_c))
    mantlef2coref[2,2] = ((-RotM[2,1,mo])*sin(ϕ_c)) + (RotM[2,2,mo]*cos(ϕ_c))
    mantlef2coref[3,2] =  RotM[2,3,mo]
    mantlef2coref[1,3] = (( RotM[3,1,mo])*cos(ϕ_c)) + (RotM[3,2,mo]*sin(ϕ_c))
    mantlef2coref[2,3] = ((-RotM[3,1,mo])*sin(ϕ_c)) + (RotM[3,2,mo]*cos(ϕ_c))
    mantlef2coref[3,3] =  RotM[3,3,mo]
    # Core angular velocity in core-equatorial frame
    ω_c_CE_1 = (mantlef2coref[1,1]*q[6N+10]) + ((mantlef2coref[1,2]*q[6N+11]) + (mantlef2coref[1,3]*q[6N+12]))
    ω_c_CE_2 = (mantlef2coref[2,1]*q[6N+10]) + ((mantlef2coref[2,2]*q[6N+11]) + (mantlef2coref[2,3]*q[6N+12]))
    ω_c_CE_3 = (mantlef2coref[3,1]*q[6N+10]) + ((mantlef2coref[3,2]*q[6N+11]) + (mantlef2coref[3,3]*q[6N+12]))

    # Second zonal harmonic coefficient
    # See Table 10 in page 50 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract                          # Earth's radius in au
    local J2E_t = (J2E + J2EDOT*(dsj2k/yr))*(RE_au^2)  # Earth (considering a linear change in time with rate J2EDOT)
    local J2S_t = JSEM[su,2]*one_t                     # Sun (static)
    # Vector of second zonal harmonic coefficients
    J2_t = Array{S}(undef, 5)
    J2_t[su] = J2S_t             # Earth
    J2_t[ea] = J2E_t             # Sun
    # Lunar torques: overall numerical factor in equation (44) in page 18 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract
    local N_MfigM_figE_factor = 7.5*μ[ea]*J2E_t

    #=
    Compute point-mass Newtonian accelerations, all bodies
    See equation (35) in page 7 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract
    =#

    Threads.@threads for j in 1:N
        # Fill point-mass Newton accelerations  with zeros
        newtonX[j] = zero_q_1
        newtonY[j] = zero_q_1
        newtonZ[j] = zero_q_1
        newtonianNb_Potential[j] = zero_q_1
        # Fill first 3N elements of dq with velocities
        dq[3j-2] = q[3(N+j)-2]
        dq[3j-1] = q[3(N+j)-1]
        dq[3j  ] = q[3(N+j)  ]
    end
    # Fill extended-body accelerations with zeros
    Threads.@threads for j in 1:N_ext
        accX[j] = zero_q_1
        accY[j] = zero_q_1
        accZ[j] = zero_q_1
    end

    Threads.@threads for j in 1:N
        for i in 1:N
            # i == j && continue
            if i == j
                continue
            else
                # Difference in position \mathbf{r_i} - \mathbf{r_j}
                X[i,j] = q[3i-2]-q[3j-2]      # X-axis component
                Y[i,j] = q[3i-1]-q[3j-1]      # Y-axis component
                Z[i,j] = q[3i]-q[3j]          # Z-axis component

                # Difference in velocity \mathbf{v_i} - \mathbf{v_j}
                U[i,j] = dq[3i-2]-dq[3j-2]    # X-axis component
                V[i,j] = dq[3i-1]-dq[3j-1]    # Y-axis component
                W[i,j] = dq[3i  ]-dq[3j  ]    # Z-axis component

                # Weighted difference in velocity 4\mathbf{v_i} - 3\mathbf{v_j}
                _4U_m_3X[i,j] = (4dq[3j-2])-(3dq[3i-2]) # X-axis component
                _4V_m_3Y[i,j] = (4dq[3j-1])-(3dq[3i-1]) # Y-axis component
                _4W_m_3Z[i,j] = (4dq[3j  ])-(3dq[3i  ]) # Z-axis component

                # Dot product inside [] in the second term
                pn2x = X[i,j]*_4U_m_3X[i,j]
                pn2y = Y[i,j]*_4V_m_3Y[i,j]
                pn2z = Z[i,j]*_4W_m_3Z[i,j]

                # Product of velocity components
                UU[i,j] = dq[3i-2]*dq[3j-2]   # v_{ix}v_{jx}
                VV[i,j] = dq[3i-1]*dq[3j-1]   # v_{iy}v_{jy}
                WW[i,j] = dq[3i  ]*dq[3j  ]   # v_{iz}v_{jz}

                # Dot product of velocities \mathbf{v_i}\cdot\mathbf{v_j}
                vi_dot_vj[i,j] = ( UU[i,j]+VV[i,j] ) + WW[i,j]

                # Distances r_{ij} = ||\mathbf{r_i} - \mathbf{r_j}||
                r_p2[i,j] = ( (X[i,j]^2)+(Y[i,j]^2) ) + (Z[i,j]^2) # r_{ij}^2
                r_p1d2[i,j] = sqrt(r_p2[i,j])                      # r_{ij}
                r_p3d2[i,j] = r_p2[i,j]^1.5                        # r_{ij}^3
                r_p7d2[i,j] = r_p2[i,j]^3.5                        # r_{ij}^7

                # Newtonian coefficient, i.e., mass parameter / distance^3 -> \mu_i / r_{ij}^3
                newtonianCoeff[i,j] =  μ[i]/r_p3d2[i,j]

                # Second term without (\mathbf{v}_i - \mathbf{v}_j)
                pn2[i,j] = newtonianCoeff[i,j]*(( pn2x+pn2y ) + pn2z)

                # Newtonian coefficient * difference between two positions, i.e.,
                # \mu_i * (\mathbf{r_i} - \mathbf{r_j}) / r_{ij}^3
                newton_acc_X[i,j] = X[i,j]*newtonianCoeff[i,j]    # X-axis component
                newton_acc_Y[i,j] = Y[i,j]*newtonianCoeff[i,j]    # Y-axis component
                newton_acc_Z[i,j] = Z[i,j]*newtonianCoeff[i,j]    # Z-axis component

                # Newtonian potential of 1 body \mu_i / r_{ij}
                newtonian1b_Potential[i,j] = μ[i]/r_p1d2[i, j]
                pn3[i,j] = 3.5newtonian1b_Potential[i,j]
                U_t_pn2[i,j] = pn2[i,j]*U[i,j]
                V_t_pn2[i,j] = pn2[i,j]*V[i,j]
                W_t_pn2[i,j] = pn2[i,j]*W[i,j]

                temp_001 = newtonX[j] + (X[i,j]*newtonianCoeff[i,j])  # X-axis component
                newtonX[j] = temp_001
                temp_002 = newtonY[j] + (Y[i,j]*newtonianCoeff[i,j])  # Y-axis component
                newtonY[j] = temp_002
                temp_003 = newtonZ[j] + (Z[i,j]*newtonianCoeff[i,j])  # Z-axis component
                newtonZ[j] = temp_003
                # Newtonian potential of N bodies
                # \sum_{i\neq l} \frac{\mu_i}{r_{il}}
                temp_004 = newtonianNb_Potential[j] + newtonian1b_Potential[i, j]
                newtonianNb_Potential[j] = temp_004
            end # else (i != j)
        end #for, i
        # Velocity magnitude squared ||\mathbf{v}_i||^2
        v2[j] = ( (dq[3j-2]^2)+(dq[3j-1]^2) ) + (dq[3j]^2)
    end #for, j

    #=
    Extended body accelerations
    See equation (28) in page 13 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract
    and equations (173) and (174) in page 33 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract
    =#

    # 2nd-order lunar zonal (J_2) and tesseral (C_2, S_2) harmonics coefficients
    # times the equatorial radius of the moon squared R_M^2
    # See equation (30) in page 13 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract
    J2M_t = ( I_M_t[3,3] - ((I_M_t[1,1]+I_M_t[2,2])/2) )/(μ[mo]) # J_{2,M}*R_M^2
    C22M_t = ((I_M_t[2,2] - I_M_t[1,1])/(μ[mo]))/4               # C_{22,M}*R_M^2
    C21M_t = (-I_M_t[1,3])/(μ[mo])                               # C_{21,M}*R_M^2
    S21M_t = (-I_M_t[3,2])/(μ[mo])                               # S_{21,M}*R_M^2
    S22M_t = ((-I_M_t[2,1])/(μ[mo]))/2                           # S_{22,M}*R_M^2
    J2_t[mo] = J2M_t

    Threads.@threads for j in 1:N_ext
        for i in 1:N_ext
            # i == j && continue
            if i == j
                continue
            else
                # J_n, C_{nm}, S_{nm} accelerations, if j-th body is flattened
                if UJ_interaction[i,j]
                    # Rotate from (X, Y, Z) inertial frame to (X_bf, Y_bf, Z_by) extended-body frame
                    X_bf_1[i,j] = X[i,j]*RotM[1,1,j]
                    X_bf_2[i,j] = Y[i,j]*RotM[1,2,j]
                    X_bf_3[i,j] = Z[i,j]*RotM[1,3,j]
                    Y_bf_1[i,j] = X[i,j]*RotM[2,1,j]
                    Y_bf_2[i,j] = Y[i,j]*RotM[2,2,j]
                    Y_bf_3[i,j] = Z[i,j]*RotM[2,3,j]
                    Z_bf_1[i,j] = X[i,j]*RotM[3,1,j]
                    Z_bf_2[i,j] = Y[i,j]*RotM[3,2,j]
                    Z_bf_3[i,j] = Z[i,j]*RotM[3,3,j]
                    X_bf[i,j] = (X_bf_1[i,j] + X_bf_2[i,j]) + (X_bf_3[i,j]) # x-coordinate in body-fixed frame
                    Y_bf[i,j] = (Y_bf_1[i,j] + Y_bf_2[i,j]) + (Y_bf_3[i,j]) # y-coordinate in body-fixed frame
                    Z_bf[i,j] = (Z_bf_1[i,j] + Z_bf_2[i,j]) + (Z_bf_3[i,j]) # z-coordinate in body-fixed frame

                    # Trigonometric functions of latitude ϕ and longitude λ in the body-fixed coordinate system
                    # See equations (165)-(168) in page 33 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract
                    sin_ϕ[i,j] = Z_bf[i,j]/r_p1d2[i,j]               # eq. (165)
                    r_xy[i,j] = sqrt( (X_bf[i,j]^2)+(Y_bf[i,j]^2) )  # X-Y projection magnitude in body-fixed frame sqrt(x_b^2 + y_b^2)
                    cos_ϕ[i,j] = r_xy[i,j]/r_p1d2[i,j]               # eq. (166)
                    sin_λ[i,j] = Y_bf[i,j]/r_xy[i,j]                 # eq. (167)
                    cos_λ[i,j] = X_bf[i,j]/r_xy[i,j]                 # eq. (168)

                    # Legendre polynomials

                    # See equations (176) and (177) in page 33 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract
                    P_n[i,j,1] = one_t       # Zeroth Legendre polynomial
                    P_n[i,j,2] = sin_ϕ[i,j]  # First Legendre polynomial
                    dP_n[i,j,1] = zero_q_1   # d/d(sin_ϕ) of zeroth Legendre polynomial
                    dP_n[i,j,2] = one_t      # d/d(sin_ϕ) of first Legendre polynomial

                    for n in 2:n1SEM[j] # min(3,n1SEM[j])
                        # Recursion relation for the n-th Legenre polynomial
                        # See equation (175) in page 33 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract
                        P_n[i,j,n+1] = ((P_n[i,j,n]*sin_ϕ[i,j])*fact1_jsem[n]) - (P_n[i,j,n-1]*fact2_jsem[n])
                        # Recursion relation for d/d(sin_ϕ) of the n-th Legendre polynomial
                        # See equation (178) in page 33 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract
                        dP_n[i,j,n+1] = (dP_n[i,j,n]*sin_ϕ[i,j]) + (P_n[i,j,n]*fact3_jsem[n])
                        # r_{ij}^{n+2}
                        temp_rn[i,j,n] = r_p1d2[i,j]^fact5_jsem[n]
                    end
                    r_p4[i,j] = r_p2[i,j]^2  # r_{ij}^4


                    # Compute accelerations due to zonal harmonics J_n

                    # Second zonal harmonic J_2
                    # See equation (28) in page 13 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract
                    # and equation (173) in page 33 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract
                    F_J_ξ[i,j] = ((P_n[i,j,3]*fact4_jsem[2])*J2_t[j])/r_p4[i,j]   # ξ-axis
                    F_J_ζ[i,j] = (((-dP_n[i,j,3])*cos_ϕ[i,j])*J2_t[j])/r_p4[i,j]  # ζ-axis
                    # Beyond third zonal harmonic J_3,...
                    F_J_ξ_36[i,j] = zero_q_1
                    F_J_ζ_36[i,j] = zero_q_1
                    for n in 3:n1SEM[j] # min(3,n1SEM[j])
                        # ξ-axis
                        temp_fjξ[i,j,n] = (((P_n[i,j,n+1]*fact4_jsem[n])*JSEM[j,n])/temp_rn[i,j,n]) + F_J_ξ_36[i,j]
                        # ζ-axis
                        temp_fjζ[i,j,n] = ((((-dP_n[i,j,n+1])*cos_ϕ[i,j])*JSEM[j,n])/temp_rn[i,j,n]) + F_J_ζ_36[i,j]
                        F_J_ξ_36[i,j] = temp_fjξ[i,j,n]
                        F_J_ζ_36[i,j] = temp_fjζ[i,j,n]
                    end

                    # Associate Legendre functions (only for the moon)
                    if j == mo
                        for m in 1:n1SEM[mo]
                            if m == 1
                                # In this case associate Legendre functions reduce to Legendre polynomials
                                # See equations (167) and (168) in page 33 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract
                                sin_mλ[i,j,1] = sin_λ[i,j] # Moyer (1971), eq. (167)
                                cos_mλ[i,j,1] = cos_λ[i,j] # Moyer (1971), eq. (168)
                                # sec( Associate Legendre polynomial with m = n = 1 )
                                # See equation (181) in page 33 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract
                                secϕ_P_nm[i,j,1,1] = one_t
                                # Associate Legendre polynomial with m = n = 1
                                P_nm[i,j,1,1] = cos_ϕ[i,j]
                                # cosϕP_1^1'
                                # See equation (183) in page 34 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract
                                # Note: the second term equation (183) vanishes when n = m
                                cosϕ_dP_nm[i,j,1,1] = sin_ϕ[i,j]*lnm3[1]
                            else
                                # Trigonometric identity sin(λ + (m - 1)λ) and cos(λ + (m - 1)λ)
                                sin_mλ[i,j,m] = (cos_mλ[i,j,m-1]*sin_mλ[i,j,1]) + (sin_mλ[i,j,m-1]*cos_mλ[i,j,1])
                                cos_mλ[i,j,m] = (cos_mλ[i,j,m-1]*cos_mλ[i,j,1]) - (sin_mλ[i,j,m-1]*sin_mλ[i,j,1])
                                # secϕ P_n^n
                                # See equation (180) in page 33 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract
                                secϕ_P_nm[i,j,m,m] = (secϕ_P_nm[i,j,m-1,m-1]*cos_ϕ[i,j])*lnm5[m]
                                # Associate Legendre polynomial with n = m
                                P_nm[i,j,m,m] = secϕ_P_nm[i,j,m,m]*cos_ϕ[i,j]
                                # cosϕP_m^m'
                                # See equation (183) in page 34 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract
                                # Note: the second term equation (183) vanishes when n = m
                                cosϕ_dP_nm[i,j,m,m] = (secϕ_P_nm[i,j,m,m]*sin_ϕ[i,j])*lnm3[m]
                            end
                            for n in m+1:n1SEM[mo]
                                # secϕ P_n^m
                                # See equation (182) in page 34 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract
                                if n == m+1
                                    secϕ_P_nm[i,j,n,m] = (secϕ_P_nm[i,j,n-1,m]*sin_ϕ[i,j])*lnm1[n,m]
                                else
                                    secϕ_P_nm[i,j,n,m] = ((secϕ_P_nm[i,j,n-1,m]*sin_ϕ[i,j])*lnm1[n,m]) + (secϕ_P_nm[i,j,n-2,m]*lnm2[n,m])
                                end
                                # Associate Legendre polynomial of degree n and order m
                                P_nm[i,j,n,m] = secϕ_P_nm[i,j,n,m]*cos_ϕ[i,j]
                                # secϕ P_n^m
                                # See equation (183) in page 34 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract
                                cosϕ_dP_nm[i,j,n,m] = ((secϕ_P_nm[i,j,n,m]*sin_ϕ[i,j])*lnm3[n]) + (secϕ_P_nm[i,j,n-1,m]*lnm4[n,m])
                            end
                        end

                        # Moon: Compute accelerations due to tesseral harmonics C_{nm}, S_{nm}
                        # See equation (28) in page 13 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract
                        # and equation (174) in page 33 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract

                        # Accelerations due to lunar tesseral harmonics C_{21} and S_{21}
                        F_CS_ξ[i,j] = (   (  (P_nm[i,j,2,1]*lnm6[2]     )*( (C21M_t*cos_mλ[i,j,1]) + (S21M_t*sin_mλ[i,j,1]) )  ) + (  (P_nm[i,j,2,2]*lnm6[2]     )*( (C22M_t*cos_mλ[i,j,2]) + (S22M_t*sin_mλ[i,j,2]) )  )   )/r_p4[i,j]
                        F_CS_η[i,j] = (   (  (secϕ_P_nm[i,j,2,1]*lnm7[1])*( (S21M_t*cos_mλ[i,j,1]) - (C21M_t*sin_mλ[i,j,1]) )  ) + (  (secϕ_P_nm[i,j,2,2]*lnm7[2])*( (S22M_t*cos_mλ[i,j,2]) - (C22M_t*sin_mλ[i,j,2]) )  )   )/r_p4[i,j]
                        F_CS_ζ[i,j] = (   (  (cosϕ_dP_nm[i,j,2,1]       )*( (C21M_t*cos_mλ[i,j,1]) + (S21M_t*sin_mλ[i,j,1]) )  ) + (  (cosϕ_dP_nm[i,j,2,2]       )*( (C22M_t*cos_mλ[i,j,2]) + (S22M_t*sin_mλ[i,j,2]) )  )   )/r_p4[i,j]
                        # Accelerations due to lunar tesseral harmonics beyond C_{21} and S_{21}
                        F_CS_ξ_36[i,j] = zero_q_1
                        F_CS_η_36[i,j] = zero_q_1
                        F_CS_ζ_36[i,j] = zero_q_1
                        for n in 3:n2M
                            for m in 1:n
                                # Lunar teseral harmonics C_{nm}/S_{nm} * trigonometric function of integer times the longitude λ
                                Cnm_cosmλ[i,j,n,m] = CM[n,m]*cos_mλ[i,j,m]
                                Cnm_sinmλ[i,j,n,m] = CM[n,m]*sin_mλ[i,j,m]
                                Snm_cosmλ[i,j,n,m] = SM[n,m]*cos_mλ[i,j,m]
                                Snm_sinmλ[i,j,n,m] = SM[n,m]*sin_mλ[i,j,m]
                                # Vector sum in equation (173)
                                temp_CS_ξ[i,j,n,m] = (   (  (P_nm[i,j,n,m]*lnm6[n]     )*( Cnm_cosmλ[i,j,n,m] + Snm_sinmλ[i,j,n,m] )  )/temp_rn[i,j,n]   ) + F_CS_ξ_36[i,j]
                                temp_CS_η[i,j,n,m] = (   (  (secϕ_P_nm[i,j,n,m]*lnm7[m])*( Snm_cosmλ[i,j,n,m] - Cnm_sinmλ[i,j,n,m] )  )/temp_rn[i,j,n]   ) + F_CS_η_36[i,j]
                                temp_CS_ζ[i,j,n,m] = (   (  (cosϕ_dP_nm[i,j,n,m]       )*( Cnm_cosmλ[i,j,n,m] + Snm_sinmλ[i,j,n,m] )  )/temp_rn[i,j,n]   ) + F_CS_ζ_36[i,j]
                                F_CS_ξ_36[i,j] = temp_CS_ξ[i,j,n,m]
                                F_CS_η_36[i,j] = temp_CS_η[i,j,n,m]
                                F_CS_ζ_36[i,j] = temp_CS_ζ[i,j,n,m]
                            end
                        end
                        # Sum the zonal and tesseral (only for the moon) accelerations without mass parameter
                        F_JCS_ξ[i,j] = (F_J_ξ[i,j] + F_J_ξ_36[i,j]) + (F_CS_ξ[i,j]+F_CS_ξ_36[i,j])
                        F_JCS_η[i,j] = (F_CS_η[i,j]+F_CS_η_36[i,j])
                        F_JCS_ζ[i,j] = (F_J_ζ[i,j] + F_J_ζ_36[i,j]) + (F_CS_ζ[i,j]+F_CS_ζ_36[i,j])
                    else
                        # Sum the zonal accelerations without mass parameter
                        F_JCS_ξ[i,j] = (F_J_ξ[i,j] + F_J_ξ_36[i,j])
                        F_JCS_η[i,j] = zero_q_1
                        F_JCS_ζ[i,j] = (F_J_ζ[i,j] + F_J_ζ_36[i,j])
                    end

                    # R matrix: body-fixed -> "primed" (ξ, η, ζ) system
                    # See equation (161) in page 32 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract
                    Rb2p[i,j,1,1] = cos_ϕ[i,j]*cos_λ[i,j]
                    Rb2p[i,j,2,1] = -sin_λ[i,j]
                    Rb2p[i,j,3,1] = -sin_ϕ[i,j]*cos_λ[i,j]
                    Rb2p[i,j,1,2] = cos_ϕ[i,j]*sin_λ[i,j]
                    Rb2p[i,j,2,2] = cos_λ[i,j]
                    Rb2p[i,j,3,2] = -sin_ϕ[i,j]*sin_λ[i,j]
                    Rb2p[i,j,1,3] = sin_ϕ[i,j]
                    Rb2p[i,j,2,3] = zero_q_1
                    Rb2p[i,j,3,3] = cos_ϕ[i,j]
                    # G matrix: space-fixed -> body-fixed -> "primed" (ξ, η, ζ) system
                    # G_{i,j} = \sum_k R_{i,k} RotM{k,j} or G = RotM * R
                    # See equation (163) in page 32 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract
                    Gc2p[i,j,1,1] = ((Rb2p[i,j,1,1]*RotM[1,1,j]) + (Rb2p[i,j,1,2]*RotM[2,1,j])) + (Rb2p[i,j,1,3]*RotM[3,1,j])
                    Gc2p[i,j,2,1] = ((Rb2p[i,j,2,1]*RotM[1,1,j]) + (Rb2p[i,j,2,2]*RotM[2,1,j])) + (Rb2p[i,j,2,3]*RotM[3,1,j])
                    Gc2p[i,j,3,1] = ((Rb2p[i,j,3,1]*RotM[1,1,j]) + (Rb2p[i,j,3,2]*RotM[2,1,j])) + (Rb2p[i,j,3,3]*RotM[3,1,j])
                    Gc2p[i,j,1,2] = ((Rb2p[i,j,1,1]*RotM[1,2,j]) + (Rb2p[i,j,1,2]*RotM[2,2,j])) + (Rb2p[i,j,1,3]*RotM[3,2,j])
                    Gc2p[i,j,2,2] = ((Rb2p[i,j,2,1]*RotM[1,2,j]) + (Rb2p[i,j,2,2]*RotM[2,2,j])) + (Rb2p[i,j,2,3]*RotM[3,2,j])
                    Gc2p[i,j,3,2] = ((Rb2p[i,j,3,1]*RotM[1,2,j]) + (Rb2p[i,j,3,2]*RotM[2,2,j])) + (Rb2p[i,j,3,3]*RotM[3,2,j])
                    Gc2p[i,j,1,3] = ((Rb2p[i,j,1,1]*RotM[1,3,j]) + (Rb2p[i,j,1,2]*RotM[2,3,j])) + (Rb2p[i,j,1,3]*RotM[3,3,j])
                    Gc2p[i,j,2,3] = ((Rb2p[i,j,2,1]*RotM[1,3,j]) + (Rb2p[i,j,2,2]*RotM[2,3,j])) + (Rb2p[i,j,2,3]*RotM[3,3,j])
                    Gc2p[i,j,3,3] = ((Rb2p[i,j,3,1]*RotM[1,3,j]) + (Rb2p[i,j,3,2]*RotM[2,3,j])) + (Rb2p[i,j,3,3]*RotM[3,3,j])
                    # Compute cartesian coordinates of acceleration due to body figure in inertial frame (without mass parameter)
                    # See equation (169) in page 33 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract
                    F_JCS_x[i,j] = ((F_JCS_ξ[i,j]*Gc2p[i,j,1,1]) + (F_JCS_η[i,j]*Gc2p[i,j,2,1])) + (F_JCS_ζ[i,j]*Gc2p[i,j,3,1])
                    F_JCS_y[i,j] = ((F_JCS_ξ[i,j]*Gc2p[i,j,1,2]) + (F_JCS_η[i,j]*Gc2p[i,j,2,2])) + (F_JCS_ζ[i,j]*Gc2p[i,j,3,2])
                    F_JCS_z[i,j] = ((F_JCS_ξ[i,j]*Gc2p[i,j,1,3]) + (F_JCS_η[i,j]*Gc2p[i,j,2,3])) + (F_JCS_ζ[i,j]*Gc2p[i,j,3,3])
                end #if UJ_interaction[i,j]
            end # else (i != j)
        end #for i in 1:N_ext
    end #for j in 1:N_ext

    for j in 1:N_ext
        for i in 1:N_ext
            # i == j && continue
            if i == j
                continue
            else
                if UJ_interaction[i,j]
                    # Extended body accelerations
                    # J_n, C_{nm}, S_{nm} accelerations, if j-th body is flattened

                    # Add result to total acceleration upon j-th body figure due to i-th point mass
                    temp_accX_j[i,j] = accX[j] - (μ[i]*F_JCS_x[i,j])
                    accX[j] = temp_accX_j[i,j]
                    temp_accY_j[i,j] = accY[j] - (μ[i]*F_JCS_y[i,j])
                    accY[j] = temp_accY_j[i,j]
                    temp_accZ_j[i,j] = accZ[j] - (μ[i]*F_JCS_z[i,j])
                    accZ[j] = temp_accZ_j[i,j]

                    # Reaction force on i-th body
                    temp_accX_i[i,j] = accX[i] + (μ[j]*F_JCS_x[i,j])
                    accX[i] = temp_accX_i[i,j]
                    temp_accY_i[i,j] = accY[i] + (μ[j]*F_JCS_y[i,j])
                    accY[i] = temp_accY_i[i,j]
                    temp_accZ_i[i,j] = accZ[i] + (μ[j]*F_JCS_z[i,j])
                    accZ[i] = temp_accZ_i[i,j]

                    # Lunar torques
                    if j == mo
                        # Compute torques acting upon the body-figure of the Moon due to external point masses
                        # See equation (43) in page 13 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract
                        N_MfigM_pmA_x[i] = μ[i]*( (Y[i,j]*F_JCS_z[i,j]) - (Z[i,j]*F_JCS_y[i,j]) )
                        N_MfigM_pmA_y[i] = μ[i]*( (Z[i,j]*F_JCS_x[i,j]) - (X[i,j]*F_JCS_z[i,j]) )
                        N_MfigM_pmA_z[i] = μ[i]*( (X[i,j]*F_JCS_y[i,j]) - (Y[i,j]*F_JCS_x[i,j]) )
                        # Expressions below have minus sign since N_MfigM_pmA_{x,y,z} have inverted signs in cross product
                        temp_N_M_x[i] = N_MfigM[1] - N_MfigM_pmA_x[i]
                        N_MfigM[1] = temp_N_M_x[i]
                        temp_N_M_y[i] = N_MfigM[2] - N_MfigM_pmA_y[i]
                        N_MfigM[2] = temp_N_M_y[i]
                        temp_N_M_z[i] = N_MfigM[3] - N_MfigM_pmA_z[i]
                        N_MfigM[3] = temp_N_M_z[i]
                    end
                end
            end # else (i != j)
        end
    end

    #=
    Post-Newtonian corrections to gravitational acceleration
    Post-Newtonian iterative procedure setup and initialization
    See equation (35) in page 7 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract
    =#
    Threads.@threads for j in 1:N
        for i in 1:N
            # i == j && continue
            if i == j
                continue
            else
                # 4*\sum term inside {}
                _4ϕj[i,j] = 4newtonianNb_Potential[j]
                # 4*\sum + \sum terms inside {}
                ϕi_plus_4ϕj[i,j] = newtonianNb_Potential[i] + _4ϕj[i,j]
                # 2 * ||\mathbf{v_i}||^2
                _2v2[i,j] = 2v2[i]
                # \dot{s}_j^2 + 2\dot{s}_i^2 inside {}
                sj2_plus_2si2[i,j] = v2[j] + _2v2[i,j]
                # \dot{s}_j^2 + 2\dot{s}_i^2 - 4<, > terms inside {}
                sj2_plus_2si2_minus_4vivj[i,j] = sj2_plus_2si2[i,j] - (4vi_dot_vj[i,j])
                # -4\sum - \sum + \dot{s}_j^2 + 2\dot{s}_i^2  - 4<, > terms inside {}
                ϕs_and_vs[i,j] = sj2_plus_2si2_minus_4vivj[i,j] - ϕi_plus_4ϕj[i,j]
                # (\mathbf{r}_i - \mathbf{r}_j)\cdot\mathbf{v_i}
                Xij_t_Ui = X[i,j]*dq[3i-2]
                Yij_t_Vi = Y[i,j]*dq[3i-1]
                Zij_t_Wi = Z[i,j]*dq[3i]
                Rij_dot_Vi = ( Xij_t_Ui+Yij_t_Vi ) + Zij_t_Wi
                # The expression below inside the (...)^2 should have a minus sign in front of the numerator,
                # but upon squaring it is eliminated, so at the end of the day, it is irrelevant ;)
                # (\mathbf{r}_i - \mathbf{r}_j)\cdot\mathbf{v_i} / r_{ij}
                pn1t7 = (Rij_dot_Vi^2)/r_p2[i,j]
                # Everything inside the {} except for the first and last terms
                pn1t2_7 = ϕs_and_vs[i,j] - (1.5pn1t7)
                # Everything inside the {} except for the last term
                pn1t1_7[i,j] = c_p2+pn1t2_7
            end # else (i != j)
        end
        # Temporary post-Newtonian accelerations
        pntempX[j] = zero_q_1   # X-axis component
        pntempY[j] = zero_q_1   # Y-axis component
        pntempZ[j] = zero_q_1   # Z-axis component
    end

    Threads.@threads for j in 1:N
        for i in 1:N
            # i == j && continue
            if i == j
                continue
            else

                # First term of equation (35) in page 7 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract

                # Last term inside the {}
                pNX_t_X[i,j] = newtonX[i]*X[i,j]   # X-axis component
                pNY_t_Y[i,j] = newtonY[i]*Y[i,j]   # Y-axis component
                pNZ_t_Z[i,j] = newtonZ[i]*Z[i,j]   # Z-axis component
                # Everything inside the {} in the first term
                pn1[i,j] = (  pn1t1_7[i,j]  +  0.5*( (pNX_t_X[i,j]+pNY_t_Y[i,j]) + pNZ_t_Z[i,j] )  )
                # Full first term
                X_t_pn1[i,j] = newton_acc_X[i,j]*pn1[i,j]   # X-axis component
                Y_t_pn1[i,j] = newton_acc_Y[i,j]*pn1[i,j]   # Y-axis component
                Z_t_pn1[i,j] = newton_acc_Z[i,j]*pn1[i,j]   # Z-axis component

                # Full third term of equation (35) in page 7 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract
                pNX_t_pn3[i,j] = newtonX[i]*pn3[i,j]   # X-axis component
                pNY_t_pn3[i,j] = newtonY[i]*pn3[i,j]   # Y-axis component
                pNZ_t_pn3[i,j] = newtonZ[i]*pn3[i,j]   # Z-axis component

                # Temporary post-Newtonian accelerations
                termpnx = ( X_t_pn1[i,j] + (U_t_pn2[i,j]+pNX_t_pn3[i,j]) )   # X-axis component
                sumpnx = pntempX[j] + termpnx
                pntempX[j] = sumpnx
                termpny = ( Y_t_pn1[i,j] + (V_t_pn2[i,j]+pNY_t_pn3[i,j]) )   # Y-axis component
                sumpny = pntempY[j] + termpny
                pntempY[j] = sumpny
                termpnz = ( Z_t_pn1[i,j] + (W_t_pn2[i,j]+pNZ_t_pn3[i,j]) )   # Z-axis component
                sumpnz = pntempZ[j] + termpnz
                pntempZ[j] = sumpnz
            end # else (i != j)
        end
        # Post-Newtonian acelerations
        postNewtonX[j] = pntempX[j]*c_m2
        postNewtonY[j] = pntempY[j]*c_m2
        postNewtonZ[j] = pntempZ[j]*c_m2
    end

    # Fill accelerations (post-Newtonian and extended body accelerations)
    # postNewton -> post-Newtonian accelerations (all bodies)
    # accX/Y/Z -> extended body accelerations (only first N_ext bodies)
    Threads.@threads for i in 1:N_ext
        dq[3(N+i)-2] = postNewtonX[i] + accX[i]
        dq[3(N+i)-1] = postNewtonY[i] + accY[i]
        dq[3(N+i)  ] = postNewtonZ[i] + accZ[i]
    end
    Threads.@threads for i in N_ext+1:N
        dq[3(N+i)-2] = postNewtonX[i]
        dq[3(N+i)-1] = postNewtonY[i]
        dq[3(N+i)  ] = postNewtonZ[i]
    end

    #=
    Lunar physical librations
    See equations (33)-(35) in pages 15-16 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract
    =#

    # Lunar moment of intertia I times angular velocity ω: Iω
    Iω_x = (I_m_t[1,1]*q[6N+4]) + ((I_m_t[1,2]*q[6N+5]) + (I_m_t[1,3]*q[6N+6])) # x-axis component
    Iω_y = (I_m_t[2,1]*q[6N+4]) + ((I_m_t[2,2]*q[6N+5]) + (I_m_t[2,3]*q[6N+6])) # y-axis component
    Iω_z = (I_m_t[3,1]*q[6N+4]) + ((I_m_t[3,2]*q[6N+5]) + (I_m_t[3,3]*q[6N+6])) # z-axis component

    # Cross product of angular velocity and Iω: ω × (I*ω)
    ωxIω_x = (q[6N+5]*Iω_z) - (q[6N+6]*Iω_y)   # x-axis component
    ωxIω_y = (q[6N+6]*Iω_x) - (q[6N+4]*Iω_z)   # y-axis component
    ωxIω_z = (q[6N+4]*Iω_y) - (q[6N+5]*Iω_x)   # z-axis component

    # Time derivative of moment of inertia times angular velocity: (dI/dt)*ω
    dIω_x = (dI_m_t[1,1]*q[6N+4]) + ((dI_m_t[1,2]*q[6N+5]) + (dI_m_t[1,3]*q[6N+6])) # x-axis component
    dIω_y = (dI_m_t[2,1]*q[6N+4]) + ((dI_m_t[2,2]*q[6N+5]) + (dI_m_t[2,3]*q[6N+6])) # y-axis component
    dIω_z = (dI_m_t[3,1]*q[6N+4]) + ((dI_m_t[3,2]*q[6N+5]) + (dI_m_t[3,3]*q[6N+6])) # z-axis component

    # Moon -> Earth radial unit vector (inertial coordinates)
    er_EM_I_1 = X[ea,mo]/r_p1d2[ea,mo]
    er_EM_I_2 = Y[ea,mo]/r_p1d2[ea,mo]
    er_EM_I_3 = Z[ea,mo]/r_p1d2[ea,mo]

    # Earth pole unit vector (inertial coordinates)
    p_E_I_1 = RotM[3,1,ea]
    p_E_I_2 = RotM[3,2,ea]
    p_E_I_3 = RotM[3,3,ea]

    # Transform Moon -> Earth radial unit vector (inertial coordinates) er_EM_I_i and
    # Earth pole unit vector p_E_I_i to lunar mantle frame coordinates
    er_EM_1 = (RotM[1,1,mo]*er_EM_I_1) + ((RotM[1,2,mo]*er_EM_I_2) + (RotM[1,3,mo]*er_EM_I_3))
    er_EM_2 = (RotM[2,1,mo]*er_EM_I_1) + ((RotM[2,2,mo]*er_EM_I_2) + (RotM[2,3,mo]*er_EM_I_3))
    er_EM_3 = (RotM[3,1,mo]*er_EM_I_1) + ((RotM[3,2,mo]*er_EM_I_2) + (RotM[3,3,mo]*er_EM_I_3))
    p_E_1 = (RotM[1,1,mo]*p_E_I_1) + ((RotM[1,2,mo]*p_E_I_2) + (RotM[1,3,mo]*p_E_I_3))
    p_E_2 = (RotM[2,1,mo]*p_E_I_1) + ((RotM[2,2,mo]*p_E_I_2) + (RotM[2,3,mo]*p_E_I_3))
    p_E_3 = (RotM[3,1,mo]*p_E_I_1) + ((RotM[3,2,mo]*p_E_I_2) + (RotM[3,3,mo]*p_E_I_3))

    # Evaluate equation (44) https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract
    # in lunar mantle frame coords

    # I*e_r
    I_er_EM_1 = (I_m_t[1,1]*er_EM_1) + ((I_m_t[1,2]*er_EM_2) + (I_m_t[1,3]*er_EM_3))
    I_er_EM_2 = (I_m_t[2,1]*er_EM_1) + ((I_m_t[2,2]*er_EM_2) + (I_m_t[2,3]*er_EM_3))
    I_er_EM_3 = (I_m_t[3,1]*er_EM_1) + ((I_m_t[3,2]*er_EM_2) + (I_m_t[3,3]*er_EM_3))

    # I*p_E
    I_p_E_1 = (I_m_t[1,1]*p_E_1) + ((I_m_t[1,2]*p_E_2) + (I_m_t[1,3]*p_E_3))
    I_p_E_2 = (I_m_t[2,1]*p_E_1) + ((I_m_t[2,2]*p_E_2) + (I_m_t[2,3]*p_E_3))
    I_p_E_3 = (I_m_t[3,1]*p_E_1) + ((I_m_t[3,2]*p_E_2) + (I_m_t[3,3]*p_E_3))

    # e_r × (I*e_r)
    er_EM_cross_I_er_EM_1 = (er_EM_2*I_er_EM_3) - (er_EM_3*I_er_EM_2)
    er_EM_cross_I_er_EM_2 = (er_EM_3*I_er_EM_1) - (er_EM_1*I_er_EM_3)
    er_EM_cross_I_er_EM_3 = (er_EM_1*I_er_EM_2) - (er_EM_2*I_er_EM_1)

    # e_r × (I*p_E)
    er_EM_cross_I_p_E_1 = (er_EM_2*I_p_E_3) - (er_EM_3*I_p_E_2)
    er_EM_cross_I_p_E_2 = (er_EM_3*I_p_E_1) - (er_EM_1*I_p_E_3)
    er_EM_cross_I_p_E_3 = (er_EM_1*I_p_E_2) - (er_EM_2*I_p_E_1)

    # p_E × (I*e_r)
    p_E_cross_I_er_EM_1 = (p_E_2*I_er_EM_3) - (p_E_3*I_er_EM_2)
    p_E_cross_I_er_EM_2 = (p_E_3*I_er_EM_1) - (p_E_1*I_er_EM_3)
    p_E_cross_I_er_EM_3 = (p_E_1*I_er_EM_2) - (p_E_2*I_er_EM_1)

    # p_E × (I*p_E)
    p_E_cross_I_p_E_1 = (p_E_2*I_p_E_3) - (p_E_3*I_p_E_2)
    p_E_cross_I_p_E_2 = (p_E_3*I_p_E_1) - (p_E_1*I_p_E_3)
    p_E_cross_I_p_E_3 = (p_E_1*I_p_E_2) - (p_E_2*I_p_E_1)

    # Coefficients of first and second terms inside {} in eq. (44)
    one_minus_7sin2ϕEM = one_t - (7((sin_ϕ[ea,mo])^2))
    two_sinϕEM = 2sin_ϕ[ea,mo]

    # Overall numerical factor in eq. (44) / r_{EM}^5
    N_MfigM_figE_factor_div_rEMp5 = (N_MfigM_figE_factor/(r_p1d2[mo,ea]^5))
    # Evaluation of eq. (44)
    N_MfigM_figE_1 = N_MfigM_figE_factor_div_rEMp5*( (one_minus_7sin2ϕEM*er_EM_cross_I_er_EM_1) + (two_sinϕEM*(er_EM_cross_I_p_E_1+p_E_cross_I_er_EM_1)) - (0.4p_E_cross_I_p_E_1))
    N_MfigM_figE_2 = N_MfigM_figE_factor_div_rEMp5*( (one_minus_7sin2ϕEM*er_EM_cross_I_er_EM_2) + (two_sinϕEM*(er_EM_cross_I_p_E_2+p_E_cross_I_er_EM_2)) - (0.4p_E_cross_I_p_E_2))
    N_MfigM_figE_3 = N_MfigM_figE_factor_div_rEMp5*( (one_minus_7sin2ϕEM*er_EM_cross_I_er_EM_3) + (two_sinϕEM*(er_EM_cross_I_p_E_3+p_E_cross_I_er_EM_3)) - (0.4p_E_cross_I_p_E_3))

    # RHS of equation (34) in page 16 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract

    # Torques acting upon lunar body-figure due to external point masses: transform coordinates from inertial frame to lunar mantle frame
    N_1_LMF = (RotM[1,1,mo]*N_MfigM[1]) + ((RotM[1,2,mo]*N_MfigM[2]) + (RotM[1,3,mo]*N_MfigM[3]))
    N_2_LMF = (RotM[2,1,mo]*N_MfigM[1]) + ((RotM[2,2,mo]*N_MfigM[2]) + (RotM[2,3,mo]*N_MfigM[3]))
    N_3_LMF = (RotM[3,1,mo]*N_MfigM[1]) + ((RotM[3,2,mo]*N_MfigM[2]) + (RotM[3,3,mo]*N_MfigM[3]))

    # Torque on the mantle due to the interaction between core and mantle (evaluated in mantle frame)
    # See equation (45) in page 18 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract
    N_cmb_1 = (k_ν*(q[6N+10]-q[6N+4])) - (C_c_m_A_c*(q[6N+12]*q[6N+11]))
    N_cmb_2 = (k_ν*(q[6N+11]-q[6N+5])) + (C_c_m_A_c*(q[6N+12]*q[6N+10]))
    N_cmb_3 = (k_ν*(q[6N+12]-q[6N+6]))

    # I*(dω/dt); i.e., I times RHS of equation (34) in page 16 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract
    I_dω_1 = ((N_MfigM_figE_1 + (μ[mo]*N_1_LMF)) + N_cmb_1) - (dIω_x + ωxIω_x)
    I_dω_2 = ((N_MfigM_figE_2 + (μ[mo]*N_2_LMF)) + N_cmb_2) - (dIω_y + ωxIω_y)
    I_dω_3 = ((N_MfigM_figE_3 + (μ[mo]*N_3_LMF)) + N_cmb_3) - (dIω_z + ωxIω_z)

    # RHS of equation (35) in page 16 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract

    # I_c * ω_c
    Ic_ωc_1 = I_c_t[1,1]*q[6N+10] # + ((I_c_t[1,2]*q[6N+11]) + (I_c_t[1,3]*q[6N+12]))
    Ic_ωc_2 = I_c_t[2,2]*q[6N+11] # + ((I_c_t[2,1]*q[6N+10]) + (I_c_t[2,3]*q[6N+12]))
    Ic_ωc_3 = I_c_t[3,3]*q[6N+12] # + ((I_c_t[3,1]*q[6N+10]) + (I_c_t[3,2]*q[6N+11]))

    # - ω_m × (I_c * ω_c)
    m_ωm_x_Icωc_1 = (q[6N+6]*Ic_ωc_2) - (q[6N+5]*Ic_ωc_3)
    m_ωm_x_Icωc_2 = (q[6N+4]*Ic_ωc_3) - (q[6N+6]*Ic_ωc_1)
    m_ωm_x_Icωc_3 = (q[6N+5]*Ic_ωc_1) - (q[6N+4]*Ic_ωc_2)

    # I_c*(dω_c/dt); i.e., I_c times RHS of of equation (35) in page 16 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract
    Ic_dωc_1 = m_ωm_x_Icωc_1 - N_cmb_1
    Ic_dωc_2 = m_ωm_x_Icωc_2 - N_cmb_2
    Ic_dωc_3 = m_ωm_x_Icωc_3 - N_cmb_3

    # Lunar mantle physical librations

    # Euler angles
    # See equation (14) in page 9 of  https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract
    dq[6N+1] = ((q[6N+4]*sin(q[6N+3])) + (q[6N+5]*cos(q[6N+3])) )/sin(q[6N+2])
    dq[6N+2] = (q[6N+4]*cos(q[6N+3])) - (q[6N+5]*sin(q[6N+3]))
    dq[6N+3] = q[6N+6] - (dq[6N+1]*cos(q[6N+2]))
    # Angular velocitiy
    # See equation (34) in page 16 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract
    dq[6N+4] = (inv_I_m_t[1,1]*I_dω_1) + ( (inv_I_m_t[1,2]*I_dω_2) + (inv_I_m_t[1,3]*I_dω_3) )
    dq[6N+5] = (inv_I_m_t[2,1]*I_dω_1) + ( (inv_I_m_t[2,2]*I_dω_2) + (inv_I_m_t[2,3]*I_dω_3) )
    dq[6N+6] = (inv_I_m_t[3,1]*I_dω_1) + ( (inv_I_m_t[3,2]*I_dω_2) + (inv_I_m_t[3,3]*I_dω_3) )

    # Lunar core physical librations

    # Euler angles
    # See equation (15) in page 9 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract
    # (core angular velocity components ω_c_CE_i represent lunar core-equator frame coordinates)
    dq[6N+9] = -(ω_c_CE_2/sin(q[6N+8])) ### evaluated first, since it's used below
    dq[6N+7] = ω_c_CE_3-(dq[6N+9]*cos(q[6N+8]))
    dq[6N+8] = ω_c_CE_1
    # Angular velocity of the core expressed in the mantle frame
    # See equation (35) in page 16 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract
    dq[6N+10] = inv_I_c_t[1,1]*Ic_dωc_1 # + ( (inv_I_c_t[1,2]*Ic_dωc_2) + (inv_I_c_t[1,3]*Ic_dωc_3) )
    dq[6N+11] = inv_I_c_t[2,2]*Ic_dωc_2 # + ( (inv_I_c_t[2,1]*Ic_dωc_1) + (inv_I_c_t[2,3]*Ic_dωc_3) )
    dq[6N+12] = inv_I_c_t[3,3]*Ic_dωc_3 # + ( (inv_I_c_t[3,1]*Ic_dωc_1) + (inv_I_c_t[3,2]*Ic_dωc_2) )

    # TT-TDB
    # TODO: implement TT-TDB integration
    dq[6N+13] = zero_q_1

    nothing
end

@doc raw"""
    DE430!(dq, q, params, t)

Solar System (JPL DE430/431) dynamical model. This function uses threads and includes all the effects in
`NBP_pN_A_J23E_J23M_J2S!` plus

- Tidal secular acceleration of Moon due to rides raised on Earth by both the Moon and the Sun: see equation (32) in page 14 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract
```math
\begin{align*}
    \mathbf{a}_{M,tide} & = \frac{3}{2}\left(\frac{m_E + m_M}{m_E}\right)\frac{Gm_T R_E^5}{r^5}
    \left\lbrace
        \frac{k_{20,E}}{r_0^*\,^5}\left(\left[2z_0^*\,^2\mathbf{z} + \rho_0^*\,^2\mathbf{\rho}\right]
        - \frac{5\left[\left(zz_0^*\right)^2 + \frac{1}{2}\left(\rho\rho_0^*\right)^2\right]\mathbf{r}}{r^2}
        + r_0^*\,^2\mathbf{r} \right) \right. \\
        & \hspace{0.5cm} + \frac{k_{21,E}}{r_1^*\,^5}\left(2\left[\left(\mathbf{\rho}\cdot\mathbf{\rho}_1^*\right)\mathbf{z}_1^*
        + zz_1^*\mathbf{\rho}_1^*\right]
        - \frac{10zz_1^*\left(\mathbf{\rho}\cdot\mathbf{\rho}_1^*\right)\mathbf{r}}{r^2}\right) \\
        & \hspace{0.5cm} + \left. \frac{k_{22,E}}{r_2^*\,^5}\left(
        \left[2\left(\mathbf{\rho}\cdot\mathbf{\rho}_2^*\right)\mathbf{\rho}_2^* - \rho_2^*\,^2\mathbf{\rho}\right]
        - \frac{5\left[\left(\mathbf{\rho}\cdot\mathbf{\rho}_2^*\right)^2 - \frac{1}{2}\left(\rho\rho_2^*\right)^2\right]\mathbf{r}}{r^2}
        \right)
    \right\rbrace,
\end{align*}
```
where ``k_{2j,E}`` with ``j = 0, 1, 2`` are the degree-2 Love numbers corresponding to tides
with long-period, diurnal, and semi-diurnal periods, respectively; ``m_E``, ``m_M`` and ``m_T`` are the masses
of the Earth, the Moon and the tide-raising body, respectively; ``R_E`` is the Earth's radius;
the position vectors ``\mathbf{r}`` and ``\mathbf{r}_j^*`` with ``j = 0, 1, 2`` are
expressed in cylindrical coordinates with the ``Z`` axis perpendicular to the Earth’s equator,
so that ``\mathbf{r} = \mathbf{\rho} + \mathbf{z}`` and the time-delayed position of the tide-raising
body is given by ``\mathbf{r}_j^* = \mathbf{\rho}_j^* + \mathbf{z}_j^*``; and ``\mathbf{a}_{M,tide}`` is
the acceleration of the Moon with respect to Earth, for each tide-raising body.

See also [`NBP_pN_A_J23E_J23M_J2S!`](@ref) and [`NBP_pN_A_J23E_J23M_J2S_threads!`](@ref).
""" DE430!
@taylorize function DE430!(dq, q, params, t)
    # N: number of bodies
    # jd0: initial Julian date
    local N, jd0 = params
    # Time Taylor variable
    local __t = Taylor1(numtype(t), t.order)
    # Type of positions/velocities components
    local S = eltype(q)
    # Zero of type S
    local zero_q_1 = zero(q[1])
    # One of the same type as time t
    local one_t = one(t)
    # Days since J2000.0 (TDB)
    local dsj2k = t+(jd0-J2000)

    # Short backward integration needed to evaluate time-delayed tidal interactions

    # Parameters
    local params_bwd = (N_bwd, jd0)
    # Positions
    local qq_bwd = Taylor1.(constant_term.(  q[ union(nbodyind(N,1:N_bwd),6N+1:6N+13) ]), t.order )::Vector{S}
    # Velocities
    local dqq_bwd = similar(qq_bwd)
    # Vector of auxiliaries
    local xaux_bwd = similar(qq_bwd)
    # Backward integration
    # TO DO: Used taylorized method instead of default jetcoeffs!
    local jc = TaylorIntegration.jetcoeffs!(NBP_pN_A_J23E_J23M_J2S_threads!, t, qq_bwd, dqq_bwd, xaux_bwd, params_bwd)

    # Evaluation of time-delayed positions
    local q_del_τ_M = special_eval(qq_bwd, __t-τ_M)   # τ_M
    local q_del_τ_0 = special_eval(qq_bwd, __t-τ_0p)  # τ_0p
    local q_del_τ_1 = special_eval(qq_bwd, __t-τ_1p)  # τ_1p
    local q_del_τ_2 = special_eval(qq_bwd, __t-τ_2p)  # τ_2p
    # Lunar mantle euler angles delayed τ_M
    local eulang_del_τ_M = q_del_τ_M[6N_bwd+1:6N_bwd+3]::Vector{S}
    # Lunar mantle angular velocity delayed τ_M
    local ω_m_del_τ_M = q_del_τ_M[6N_bwd+4:6N_bwd+6]::Vector{S}

    # Matrix elements of lunar mantle moment of inertia at time t-τ_M (including tidal distortion)
    # See equations (36) to (41) in pages 16-17 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract
    local I_m_t = ITM(q_del_τ_M, eulang_del_τ_M, ω_m_del_τ_M)::Matrix{S}
    local dI_m_t = ordpres_differentiate.(I_m_t) # Time-derivative of lunar mantle I at time t-τ_M
    local inv_I_m_t = inv(I_m_t)                 # Inverse of lunar mantle I matrix at time t-τ_M
    local I_c_t = I_c.*one_t                     # Lunar core I matrix, see equation (39)
    local inv_I_c_t = inv(I_c_t)                 # Inverse of lunar core I matrix
    local I_M_t = I_m_t+I_c_t                    # Total I matrix (mantle + core)

    #=
    Point-mass accelerations
    See equation (35) in page 7 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract
    =#

    # Note: All the following arrays are declared here in order to help @taylorize work

    # Difference between two positions (\mathbf{r}_i - \mathbf{r}_j)
    X = Array{S}(undef, N, N)         # X-axis component
    Y = Array{S}(undef, N, N)         # Y-axis component
    Z = Array{S}(undef, N, N)         # Z-axis component

    # Distance between two positions r_{ij} = ||\mathbf{r}_i - \mathbf{r}_j||
    r_p2 = Array{S}(undef, N, N)      # r_{ij}^2
    r_p1d2 = Array{S}(undef, N, N)    # sqrt(r_p2) <-> r_{ij}
    r_p3d2 = Array{S}(undef, N, N)    # r_p2^1.5 <-> r_{ij}^3
    r_p7d2 = Array{S}(undef, N, N)    # r_p2^3.5 <-> r_{ij}^7

    # Newtonian accelerations \mathbf{a}_{i} = \sum_{i\neq j} mu_i * (\mathbf{r_i} - \mathbf{r_j}) / r_{ij}^3
    newtonX = Array{S}(undef, N)      # X-axis component
    newtonY = Array{S}(undef, N)      # Y-axis component
    newtonZ = Array{S}(undef, N)      # Z-axis component
    # Newtonian coefficient, i.e., mass parameter / distance^3 -> \mu_i / r_{ij}^3
    newtonianCoeff = Array{S}(undef, N, N)

    # Post-Newtonian stuff

    # Difference between two velocities (\mathbf{v}_i - \mathbf{v}_j)
    U = Array{S}(undef, N, N)         # X-axis component
    V = Array{S}(undef, N, N)         # Y-axis component
    W = Array{S}(undef, N, N)         # Z-axis component

    # Weighted difference between two velocities (4\mathbf{v}_i - 3\mathbf{v}_j)
    _4U_m_3X = Array{S}(undef, N, N)  # X-axis component
    _4V_m_3Y = Array{S}(undef, N, N)  # Y-axis component
    _4W_m_3Z = Array{S}(undef, N, N)  # Z-axis component

    # Product of velocity components
    UU = Array{S}(undef, N, N)        # v_{ix}v_{jx}
    VV = Array{S}(undef, N, N)        # v_{iy}v_{jy}
    WW = Array{S}(undef, N, N)        # v_{iz}v_{jz}

    # Newtonian potential of 1 body \mu_i / r_{ij}
    newtonian1b_Potential = Array{S}(undef, N, N)
    # Newtonian potential of N bodies
    # \sum_{i\neq l} \frac{\mu_i}{r_{il}}
    newtonianNb_Potential = Array{S}(undef, N)

    # Newtonian coefficient * difference between two positions, i.e.,
    # \mu_i * (\mathbf{r_i} - \mathbf{r_j}) / r_{ij}^3
    newton_acc_X = Array{S}(undef, N, N)
    newton_acc_Y = Array{S}(undef, N, N)
    newton_acc_Z = Array{S}(undef, N, N)

    # Combinations of velocities
    v2 = Array{S}(undef, N)                       # Velocity magnitude squared ||\mathbf{v}_i||^2
    _2v2 = Array{S}(undef, N, N)                  # 2 * ||\mathbf{v_i}||^2
    vi_dot_vj = Array{S}(undef, N, N)             # Dot product of two velocities \mathbf{v}_i\cdot\mathbf{v}_j
    rij_dot_vi_div_rij_sq = Array{S}(undef, N, N) # || \mathbf{r}_{i,j} \cdot \mathbf{v}_j ||^2  /  ||\mathbf{r}_{i,j}||^2

    # Second term of equation (35) in page 7 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract

    # Second term without (\mathbf{v}_i - \mathbf{v}_j)
    pn2 = Array{S}(undef, N, N)            # \mu_i * [(\mathbf{r_i} - \mathbf{r_j})\cdot(4\mathbf{v_i} - 3\mathbf{v_j})]
    # Full second term
    U_t_pn2 = Array{S}(undef, N, N)        # X-axis component
    V_t_pn2 = Array{S}(undef, N, N)        # Y-axis component
    W_t_pn2 = Array{S}(undef, N, N)        # Z-axis component

    # Third term of equation (35) in page 7 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract

    # Third term without Newtonian accelerations \mathbf{a}_i
    pn3 = Array{S}(undef, N, N)
    # Full third term of equation (35)
    pNX_t_pn3 = Array{S}(undef, N, N)      # X-axis component
    pNY_t_pn3 = Array{S}(undef, N, N)      # Y-axis component
    pNZ_t_pn3 = Array{S}(undef, N, N)      # Z-axis component

    # First term of equation (35) in page 7 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract

    _4ϕj = Array{S}(undef, N, N)            # 4*\sum term inside {}
    ϕi_plus_4ϕj = Array{S}(undef, N, N)     # 4*\sum + \sum terms inside {}
    sj2_plus_2si2 = Array{S}(undef, N, N)   # \dot{s}_j^2 + 2\dot{s}_i^2 inside {}
    sj2_plus_2si2_minus_4vivj = Array{S}(undef, N, N)  # \dot{s}_j^2 + 2\dot{s}_i^2 - 4<, > terms inside {}
    ϕs_and_vs = Array{S}(undef, N, N)       # -4\sum - \sum + \dot{s}_j^2 + 2\dot{s}_i^2 - 4<, > terms inside {}
    pn1t1_7 = Array{S}(undef, N, N)         # Everything inside the {} in the first term except for the term with accelerations (last)
    # Last term inside the {}
    pNX_t_X = Array{S}(undef, N, N)     # X-axis component
    pNY_t_Y = Array{S}(undef, N, N)     # Y-axis component
    pNZ_t_Z = Array{S}(undef, N, N)     # Z-axis component
    # Everything inside the {} in the first term
    pn1 = Array{S}(undef, N, N)
    # Full first term
    X_t_pn1 = Array{S}(undef, N, N)    # X-axis component
    Y_t_pn1 = Array{S}(undef, N, N)    # Y-axis component
    Z_t_pn1 = Array{S}(undef, N, N)    # Z-axis component

    # Temporary post-Newtonian accelerations
    pntempX = Array{S}(undef, N)        # X-axis component
    pntempY = Array{S}(undef, N)        # Y-axis component
    pntempZ = Array{S}(undef, N)        # Z-axis component
    # Full post-Newtonian accelerations
    postNewtonX = Array{S}(undef, N)    # X-axis component
    postNewtonY = Array{S}(undef, N)    # Y-axis component
    postNewtonZ = Array{S}(undef, N)    # Z-axis component

    #=
    Extended body accelerations
    See equation (28) in page 13 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract
    and equations (173) and (174) in page 33 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract
    =#

    # (J_n, C_{nm}, S_{nm}) acceleration auxiliaries

    # Auxiliaries to compute body-fixed frame coordinates
    X_bf_1 = Array{S}(undef, N_ext, N_ext)
    Y_bf_1 = Array{S}(undef, N_ext, N_ext)
    Z_bf_1 = Array{S}(undef, N_ext, N_ext)
    X_bf_2 = Array{S}(undef, N_ext, N_ext)
    Y_bf_2 = Array{S}(undef, N_ext, N_ext)
    Z_bf_2 = Array{S}(undef, N_ext, N_ext)
    X_bf_3 = Array{S}(undef, N_ext, N_ext)
    Y_bf_3 = Array{S}(undef, N_ext, N_ext)
    Z_bf_3 = Array{S}(undef, N_ext, N_ext)
    # Body-fixed frame coordinates
    X_bf = Array{S}(undef, N_ext, N_ext)
    Y_bf = Array{S}(undef, N_ext, N_ext)
    Z_bf = Array{S}(undef, N_ext, N_ext)

    # Extended body accelerations (without mass parameter) in the inertial frame
    F_JCS_x = Array{S}(undef, N_ext, N_ext)
    F_JCS_y = Array{S}(undef, N_ext, N_ext)
    F_JCS_z = Array{S}(undef, N_ext, N_ext)
    # Temporary arrays for the sum of full extended body accelerations
    temp_accX_j = Array{S}(undef, N_ext, N_ext)
    temp_accY_j = Array{S}(undef, N_ext, N_ext)
    temp_accZ_j = Array{S}(undef, N_ext, N_ext)
    temp_accX_i = Array{S}(undef, N_ext, N_ext)
    temp_accY_i = Array{S}(undef, N_ext, N_ext)
    temp_accZ_i = Array{S}(undef, N_ext, N_ext)

    # Trigonometric functions of latitude ϕ and longitude λ in the body-fixed coordinate system
    # See equations (165)-(168) in page 33 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract
    sin_ϕ = Array{S}(undef, N_ext, N_ext)
    cos_ϕ = Array{S}(undef, N_ext, N_ext)
    sin_λ = Array{S}(undef, N_ext, N_ext)
    cos_λ = Array{S}(undef, N_ext, N_ext)

    # Distances
    r_xy = Array{S}(undef, N_ext, N_ext)  # X-Y projection magnitude in body-fixed frame sqrt(x_b^2 + y_b^2)
    r_p4 = Array{S}(undef, N_ext, N_ext)  # r_{ij}^4
    # Legendre polynomials
    P_n = Array{S}(undef, N_ext, N_ext, maximum(n1SEM)+1)   # Vector of Legendre polynomials
    dP_n = Array{S}(undef, N_ext, N_ext, maximum(n1SEM)+1)  # Vector of d/d(sin ϕ) of Legendre polynomials

    # Temporary arrays for the sum of accelerations due to zonal harmonics J_n
    temp_fjξ = Array{S}(undef, N_ext, N_ext, maximum(n1SEM)+1)       # ξ-axis component
    temp_fjζ = Array{S}(undef, N_ext, N_ext, maximum(n1SEM)+1)       # ζ-axis component
    temp_rn = Array{S}(undef, N_ext, N_ext, maximum(n1SEM)+1)        # r_{ij}^{n+2}
    # Temporary arrays for the vector sum in equation (173) in page 33 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract
    temp_CS_ξ = Array{S}(undef, N_ext, N_ext, n1SEM[mo], n1SEM[mo])  # ξ-axis component
    temp_CS_η = Array{S}(undef, N_ext, N_ext, n1SEM[mo], n1SEM[mo])  # η-axis component
    temp_CS_ζ = Array{S}(undef, N_ext, N_ext, n1SEM[mo], n1SEM[mo])  # ζ-axis component
    # Accelerations due to lunar tesseral harmonics beyond C_{21} and S_{21}
    F_CS_ξ_36 = Array{S}(undef, N_ext, N_ext)  # ξ-axis component
    F_CS_η_36 = Array{S}(undef, N_ext, N_ext)  # η-axis component
    F_CS_ζ_36 = Array{S}(undef, N_ext, N_ext)  # ζ-axis component
    # Accelerations due to third zonal harmonic and beyond
    F_J_ξ_36 = Array{S}(undef, N_ext, N_ext)   # ξ-axis component
    F_J_ζ_36 = Array{S}(undef, N_ext, N_ext)   # ζ-axis component

    # Trigonometric functions of integer multiples the longitude λ in the body-fixed coordinate system
    # See equations (165)-(168) in page 33 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract
    sin_mλ = Array{S}(undef, N_ext, N_ext, n1SEM[mo])
    cos_mλ = Array{S}(undef, N_ext, N_ext, n1SEM[mo])
    # Lunar teseral harmonics C_{nm}/S_{nm} * trigonometric function of integer times the longitude λ
    Cnm_cosmλ = Array{S}(undef, N_ext, N_ext, n1SEM[mo]+1, n1SEM[mo]+1)
    Cnm_sinmλ = Array{S}(undef, N_ext, N_ext, n1SEM[mo]+1, n1SEM[mo]+1)
    Snm_cosmλ = Array{S}(undef, N_ext, N_ext, n1SEM[mo]+1, n1SEM[mo]+1)
    Snm_sinmλ = Array{S}(undef, N_ext, N_ext, n1SEM[mo]+1, n1SEM[mo]+1)

    # Associated Legendre functions
    secϕ_P_nm = Array{S}(undef, N_ext, N_ext, n1SEM[mo]+1, n1SEM[mo]+1)   # secϕ P_n^m
    P_nm = Array{S}(undef, N_ext, N_ext, n1SEM[mo]+1, n1SEM[mo]+1)        # Vector of associated Legendre functions
    cosϕ_dP_nm = Array{S}(undef, N_ext, N_ext, n1SEM[mo]+1, n1SEM[mo]+1)  # cosϕ d/d(sin ϕ)P_n^m# Accelerations due to second zonal harmonic
    # Accelerations due to second zonal harmonic
    F_J_ξ = Array{S}(undef, N_ext, N_ext)   # ξ-axis component
    F_J_η = Array{S}(undef, N_ext, N_ext)   # η-axis component
    F_J_ζ = Array{S}(undef, N_ext, N_ext)   # ζ-axis component
    # Accelerations due to lunar tesseral harmonics C_{21} and S_{21}
    F_CS_ξ = Array{S}(undef, N_ext, N_ext)  # ξ-axis component
    F_CS_η = Array{S}(undef, N_ext, N_ext)  # η-axis component
    F_CS_ζ = Array{S}(undef, N_ext, N_ext)  # ζ-axis component
    # Sum of the zonal and tesseral (only for the moon) accelerations without mass parameter
    # in body-fixed frame
    F_JCS_ξ = Array{S}(undef, N_ext, N_ext) # ξ-axis component
    F_JCS_η = Array{S}(undef, N_ext, N_ext) # η-axis component
    F_JCS_ζ = Array{S}(undef, N_ext, N_ext) # ζ-axis component

    # Rotation matrices

    # R matrix body-fixed -> "primed" (ξ, η, ζ) frame
    # See equation (161) in page 32 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract
    Rb2p = Array{S}(undef, N_ext, N_ext, 3, 3)
    # G matrix "space-fixed" -> "primed" (ξ, η, ζ) frame
    # See equation (163) in page 32 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract
    Gc2p = Array{S}(undef, N_ext, N_ext, 3, 3)

    # Full extended-body accelerations
    accX = Array{S}(undef, N_ext)
    accY = Array{S}(undef, N_ext)
    accZ = Array{S}(undef, N_ext)

    # Lunar torques
    # See equation (43) in page 18 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract
    # Vector of lunar torques
    N_MfigM_pmA_x = Array{S}(undef, N_ext)   # x-axis component
    N_MfigM_pmA_y = Array{S}(undef, N_ext)   # y-axis component
    N_MfigM_pmA_z = Array{S}(undef, N_ext)   # z-axis component
    # Temporary array for the sum of lunar torques
    temp_N_M_x = Array{S}(undef, N_ext)       # x-axis component
    temp_N_M_y = Array{S}(undef, N_ext)       # y-axis component
    temp_N_M_z = Array{S}(undef, N_ext)       # z-axis component
    # Total lunar torque
    N_MfigM = Array{S}(undef, 3)
    N_MfigM[1] = zero_q_1                     # x-axis component
    N_MfigM[2] = zero_q_1                     # y-axis component
    N_MfigM[3] = zero_q_1                     # z-axis component

    # Rotations to and from Earth, Sun and Moon pole-oriented frames
    local αs = deg2rad(α_p_sun*one_t)           # Sun's rotation pole right ascension (radians)
    local δs = deg2rad(δ_p_sun*one_t)           # Sun's rotation pole right ascension (radians)
    # Space-fixed -> body-fixed coordinate transformations
    RotM = Array{S}(undef, 3, 3, 5)
    local RotM[:,:,ea] = c2t_jpl_de430(dsj2k)   # Earth
    local RotM[:,:,su] = pole_rotation(αs, δs)  # Sun
    # Lunar mantle Euler angles
    ϕ_m = q[6N+1]
    θ_m = q[6N+2]
    ψ_m = q[6N+3]
    # Lunar mantle space-fixed -> body-fixed coodinate transformations
    # See equations (10)-(13) in page 9 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract
    RotM[1,1,mo] = (cos(ϕ_m)*cos(ψ_m)) - (cos(θ_m)*(sin(ϕ_m)*sin(ψ_m)))
    RotM[2,1,mo] = (-cos(θ_m)*(cos(ψ_m)*sin(ϕ_m))) - (cos(ϕ_m)*sin(ψ_m))
    RotM[3,1,mo] = sin(θ_m)*sin(ϕ_m)
    RotM[1,2,mo] = (cos(ψ_m)*sin(ϕ_m)) + (cos(θ_m)*(cos(ϕ_m)*sin(ψ_m)))
    RotM[2,2,mo] = (cos(θ_m)*(cos(ϕ_m)*cos(ψ_m))) - (sin(ϕ_m)*sin(ψ_m))
    RotM[3,2,mo] = (-cos(ϕ_m))*sin(θ_m)
    RotM[1,3,mo] = sin(θ_m)*sin(ψ_m)
    RotM[2,3,mo] = cos(ψ_m)*sin(θ_m)
    RotM[3,3,mo] = cos(θ_m)
    # Lunar mantle frame -> inertial frame -> Lunar core-equatorial frame coord transformation
    # See equation (16) in page 10 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract
    mantlef2coref = Array{S}(undef, 3, 3) # lunar mantle frame -> inertial frame -> lunar core-equatorial frame coord transformation
    # Lunar core Euler angle
    ϕ_c = q[6N+7]
    # mantlef2coref = R_z(ϕ_c)*[ R_z(ψ_m)*R_x(θ_m)*R_z(ϕ_m) ]^T
    mantlef2coref[1,1] = (( RotM[1,1,mo])*cos(ϕ_c)) + (RotM[1,2,mo]*sin(ϕ_c))
    mantlef2coref[2,1] = ((-RotM[1,1,mo])*sin(ϕ_c)) + (RotM[1,2,mo]*cos(ϕ_c))
    mantlef2coref[3,1] =  RotM[1,3,mo]
    mantlef2coref[1,2] = (( RotM[2,1,mo])*cos(ϕ_c)) + (RotM[2,2,mo]*sin(ϕ_c))
    mantlef2coref[2,2] = ((-RotM[2,1,mo])*sin(ϕ_c)) + (RotM[2,2,mo]*cos(ϕ_c))
    mantlef2coref[3,2] =  RotM[2,3,mo]
    mantlef2coref[1,3] = (( RotM[3,1,mo])*cos(ϕ_c)) + (RotM[3,2,mo]*sin(ϕ_c))
    mantlef2coref[2,3] = ((-RotM[3,1,mo])*sin(ϕ_c)) + (RotM[3,2,mo]*cos(ϕ_c))
    mantlef2coref[3,3] =  RotM[3,3,mo]
    # Core angular velocity in core-equatorial frame
    ω_c_CE_1 = (mantlef2coref[1,1]*q[6N+10]) + ((mantlef2coref[1,2]*q[6N+11]) + (mantlef2coref[1,3]*q[6N+12]))
    ω_c_CE_2 = (mantlef2coref[2,1]*q[6N+10]) + ((mantlef2coref[2,2]*q[6N+11]) + (mantlef2coref[2,3]*q[6N+12]))
    ω_c_CE_3 = (mantlef2coref[3,1]*q[6N+10]) + ((mantlef2coref[3,2]*q[6N+11]) + (mantlef2coref[3,3]*q[6N+12]))

    # Second zonal harmonic coefficient
    # See Table 10 in page 50 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract                          # Earth's radius in au
    local J2E_t = (J2E + J2EDOT*(dsj2k/yr))*(RE_au^2)  # Earth (considering a linear change in time with rate J2EDOT)
    local J2S_t = JSEM[su,2]*one_t                     # Sun (static)
    # Vector of second zonal harmonic coefficients
    J2_t = Array{S}(undef, 5)
    J2_t[su] = J2S_t             # Earth
    J2_t[ea] = J2E_t             # Sun
    # Lunar torques: overall numerical factor in equation (44) in page 18 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract
    local N_MfigM_figE_factor = 7.5*μ[ea]*J2E_t

    #=
    Tidal accelerations
    See equation (32) in page 14 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract
    =#

    # Time-delayed geocentric Moon position
    # See equation (31) in page 14 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract
    local q_ME_τ_0 = q_del_τ_0[3mo-2:3mo] .- q_del_τ_0[3ea-2:3ea]
    local q_ME_τ_1 = q_del_τ_1[3mo-2:3mo] .- q_del_τ_1[3ea-2:3ea]
    local q_ME_τ_2 = q_del_τ_2[3mo-2:3mo] .- q_del_τ_2[3ea-2:3ea]
    # Time-delayed geocentric Sun position
    local q_SE_τ_0 = q_del_τ_0[3su-2:3su] .- q_del_τ_0[3ea-2:3ea]
    local q_SE_τ_1 = q_del_τ_1[3su-2:3su] .- q_del_τ_1[3ea-2:3ea]
    local q_SE_τ_2 = q_del_τ_2[3su-2:3su] .- q_del_τ_2[3ea-2:3ea]
    # R3X: geocentric space-fixed -> rotational time-delay -> geocentric Earth true-equator-of-date frame
    local R30 = RotM[:,:,ea]               # *Rz(-ω_E*τ_0) == Id(3x3), since τ_0=0
    local R31 = Rz(-ω_E*τ_1)*RotM[:,:,ea]
    local R32 = Rz(-ω_E*τ_2)*RotM[:,:,ea]

    # Position vectors of tide-raising body in geocentric Earth true-equator-of-date frame
    local r_star_M_0 = R30*q_ME_τ_0 # r-star 0, Moon
    local r_star_M_1 = R31*q_ME_τ_1 # r-star 1, Moon
    local r_star_M_2 = R32*q_ME_τ_2 # r-star 2, Moon
    local r_star_S_0 = R30*q_SE_τ_0 # r-star 0, Sun
    local r_star_S_1 = R31*q_SE_τ_1 # r-star 1, Sun
    local r_star_S_2 = R32*q_SE_τ_2 # r-star 2, Sun

    #=
    Compute point-mass Newtonian accelerations, all bodies
    See equation (35) in page 7 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract
    =#
    for j in 1:N
        # Fill point-mass Newton accelerations  with zeros
        newtonX[j] = zero_q_1
        newtonY[j] = zero_q_1
        newtonZ[j] = zero_q_1
        newtonianNb_Potential[j] = zero_q_1
        # Fill first 3N elements of dq with velocities
        dq[3j-2] = q[3(N+j)-2]
        dq[3j-1] = q[3(N+j)-1]
        dq[3j  ] = q[3(N+j)  ]
    end
    # Fill extended-body accelerations with zeros
    for j in 1:N_ext
        accX[j] = zero_q_1
        accY[j] = zero_q_1
        accZ[j] = zero_q_1
    end

    Threads.@threads for j in 1:N
        for i in 1:N
            # i == j && continue
            if i == j
                continue
            else
                # Difference in position \mathbf{r_i} - \mathbf{r_j}
                X[i,j] = q[3i-2]-q[3j-2]      # X-axis component
                Y[i,j] = q[3i-1]-q[3j-1]      # Y-axis component
                Z[i,j] = q[3i]-q[3j]          # Z-axis component

                # Difference in velocity \mathbf{v_i} - \mathbf{v_j}
                U[i,j] = dq[3i-2]-dq[3j-2]    # X-axis component
                V[i,j] = dq[3i-1]-dq[3j-1]    # Y-axis component
                W[i,j] = dq[3i  ]-dq[3j  ]    # Z-axis component

                # Weighted difference in velocity 4\mathbf{v_i} - 3\mathbf{v_j}
                _4U_m_3X[i,j] = (4dq[3j-2])-(3dq[3i-2]) # X-axis component
                _4V_m_3Y[i,j] = (4dq[3j-1])-(3dq[3i-1]) # Y-axis component
                _4W_m_3Z[i,j] = (4dq[3j  ])-(3dq[3i  ]) # Z-axis component

                # Dot product inside [] in the second term
                pn2x = X[i,j]*_4U_m_3X[i,j]
                pn2y = Y[i,j]*_4V_m_3Y[i,j]
                pn2z = Z[i,j]*_4W_m_3Z[i,j]

                # Product of velocity components
                UU[i,j] = dq[3i-2]*dq[3j-2]   # v_{ix}v_{jx}
                VV[i,j] = dq[3i-1]*dq[3j-1]   # v_{iy}v_{jy}
                WW[i,j] = dq[3i  ]*dq[3j  ]   # v_{iz}v_{jz}

                # Dot product of velocities \mathbf{v_i}\cdot\mathbf{v_j}
                vi_dot_vj[i,j] = ( UU[i,j]+VV[i,j] ) + WW[i,j]

                # Distances r_{ij} = ||\mathbf{r_i} - \mathbf{r_j}||
                r_p2[i,j] = ( (X[i,j]^2)+(Y[i,j]^2) ) + (Z[i,j]^2) # r_{ij}^2
                r_p1d2[i,j] = sqrt(r_p2[i,j])                      # r_{ij}
                r_p3d2[i,j] = r_p2[i,j]^1.5                        # r_{ij}^3
                r_p7d2[i,j] = r_p2[i,j]^3.5                        # r_{ij}^7

                # Newtonian coefficient, i.e., mass parameter / distance^3 -> \mu_i / r_{ij}^3
                newtonianCoeff[i,j] =  μ[i]/r_p3d2[i,j]

                # Second term without (\mathbf{v}_i - \mathbf{v}_j)
                pn2[i,j] = newtonianCoeff[i,j]*(( pn2x+pn2y ) + pn2z)

                # Newtonian coefficient * difference between two positions, i.e.,
                # \mu_i * (\mathbf{r_i} - \mathbf{r_j}) / r_{ij}^3
                newton_acc_X[i,j] = X[i,j]*newtonianCoeff[i,j]    # X-axis component
                newton_acc_Y[i,j] = Y[i,j]*newtonianCoeff[i,j]    # Y-axis component
                newton_acc_Z[i,j] = Z[i,j]*newtonianCoeff[i,j]    # Z-axis component

                # Newtonian potential of 1 body \mu_i / r_{ij}
                newtonian1b_Potential[i,j] = μ[i]/r_p1d2[i, j]
                # Third term without newtonian accelerations \mathbf{a}_i
                pn3[i,j] = 3.5newtonian1b_Potential[i,j]
                # Full second term
                U_t_pn2[i,j] = pn2[i,j]*U[i,j]   # X-axis component
                V_t_pn2[i,j] = pn2[i,j]*V[i,j]   # Y-axis component
                W_t_pn2[i,j] = pn2[i,j]*W[i,j]   # Z-axis component

                # Newtonian accelerations \mathbf{a}_{i} = \sum_{i\neq j} mu_i * (\mathbf{r_i} - \mathbf{r_j}) / r_{ij}^3
                temp_001 = newtonX[j] + (X[i,j]*newtonianCoeff[i,j])  # X-axis component
                newtonX[j] = temp_001
                temp_002 = newtonY[j] + (Y[i,j]*newtonianCoeff[i,j])  # Y-axis component
                newtonY[j] = temp_002
                temp_003 = newtonZ[j] + (Z[i,j]*newtonianCoeff[i,j])  # Z-axis component
                newtonZ[j] = temp_003
                # Newtonian potential of N bodies
                # \sum_{i\neq l} \frac{\mu_i}{r_{il}}
                temp_004 = newtonianNb_Potential[j] + newtonian1b_Potential[i, j]
                newtonianNb_Potential[j] = temp_004
            end # else (i != j)
        end #for, i
        # Velocity magnitude squared ||\mathbf{v}_i||^2
        v2[j] = ( (dq[3j-2]^2)+(dq[3j-1]^2) ) + (dq[3j]^2)
    end #for, j

    #=
    Extended body accelerations
    See equation (28) in page 13 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract
    and equations (173) and (174) in page 33 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract
    =#

    # 2nd-order lunar zonal (J_2) and tesseral (C_2, S_2) harmonics coefficients
    # times the equatorial radius of the moon squared R_M^2
    # See equation (30) in page 13 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract
    J2M_t = ( I_M_t[3,3] - ((I_M_t[1,1]+I_M_t[2,2])/2) )/(μ[mo]) # J_{2,M}*R_M^2
    C22M_t = ((I_M_t[2,2] - I_M_t[1,1])/(μ[mo]))/4               # C_{22,M}*R_M^2
    C21M_t = (-I_M_t[1,3])/(μ[mo])                               # C_{21,M}*R_M^2
    S21M_t = (-I_M_t[3,2])/(μ[mo])                               # S_{21,M}*R_M^2
    S22M_t = ((-I_M_t[2,1])/(μ[mo]))/2                           # S_{22,M}*R_M^2
    J2_t[mo] = J2M_t

    Threads.@threads for j in 1:N_ext
        for i in 1:N_ext
            # i == j && continue
            if i == j
                continue
            else
                # J_n, C_{nm}, S_{nm} accelerations, if j-th body is flattened
                if UJ_interaction[i,j]
                    # Rotate from (X, Y, Z) inertial frame to (X_bf, Y_bf, Z_by) extended-body frame
                    X_bf_1[i,j] = X[i,j]*RotM[1,1,j]
                    X_bf_2[i,j] = Y[i,j]*RotM[1,2,j]
                    X_bf_3[i,j] = Z[i,j]*RotM[1,3,j]
                    Y_bf_1[i,j] = X[i,j]*RotM[2,1,j]
                    Y_bf_2[i,j] = Y[i,j]*RotM[2,2,j]
                    Y_bf_3[i,j] = Z[i,j]*RotM[2,3,j]
                    Z_bf_1[i,j] = X[i,j]*RotM[3,1,j]
                    Z_bf_2[i,j] = Y[i,j]*RotM[3,2,j]
                    Z_bf_3[i,j] = Z[i,j]*RotM[3,3,j]
                    X_bf[i,j] = (X_bf_1[i,j] + X_bf_2[i,j]) + (X_bf_3[i,j]) # x-coordinate in body-fixed frame
                    Y_bf[i,j] = (Y_bf_1[i,j] + Y_bf_2[i,j]) + (Y_bf_3[i,j]) # y-coordinate in body-fixed frame
                    Z_bf[i,j] = (Z_bf_1[i,j] + Z_bf_2[i,j]) + (Z_bf_3[i,j]) # z-coordinate in body-fixed frame

                    # Trigonometric functions of latitude ϕ and longitude λ in the body-fixed coordinate system
                    # See equations (165)-(168) in page 33 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract
                    sin_ϕ[i,j] = Z_bf[i,j]/r_p1d2[i,j]               # eq. (165)
                    r_xy[i,j] = sqrt( (X_bf[i,j]^2)+(Y_bf[i,j]^2) )  # X-Y projection magnitude in body-fixed frame sqrt(x_b^2 + y_b^2)
                    cos_ϕ[i,j] = r_xy[i,j]/r_p1d2[i,j]               # eq. (166)
                    sin_λ[i,j] = Y_bf[i,j]/r_xy[i,j]                 # eq. (167)
                    cos_λ[i,j] = X_bf[i,j]/r_xy[i,j]                 # eq. (168)

                    # Legendre polynomials

                    # See equations (176) and (177) in page 33 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract
                    P_n[i,j,1] = one_t       # Zeroth Legendre polynomial
                    P_n[i,j,2] = sin_ϕ[i,j]  # First Legendre polynomial
                    dP_n[i,j,1] = zero_q_1   # d/d(sin_ϕ) of zeroth Legendre polynomial
                    dP_n[i,j,2] = one_t      # d/d(sin_ϕ) of first Legendre polynomial

                    for n in 2:n1SEM[j] # min(3,n1SEM[j])
                        # Recursion relation for the n-th Legenre polynomial
                        # See equation (175) in page 33 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract
                        P_n[i,j,n+1] = ((P_n[i,j,n]*sin_ϕ[i,j])*fact1_jsem[n]) - (P_n[i,j,n-1]*fact2_jsem[n])
                        # Recursion relation for d/d(sin_ϕ) of the n-th Legendre polynomial
                        # See equation (178) in page 33 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract
                        dP_n[i,j,n+1] = (dP_n[i,j,n]*sin_ϕ[i,j]) + (P_n[i,j,n]*fact3_jsem[n])
                        # r_{ij}^{n+2}
                        temp_rn[i,j,n] = r_p1d2[i,j]^fact5_jsem[n]
                    end
                    r_p4[i,j] = r_p2[i,j]^2  # r_{ij}^4

                    # Compute accelerations due to zonal harmonics J_n

                    # Second zonal harmonic J_2
                    # See equation (28) in page 13 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract
                    # and equation (173) in page 33 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract
                    F_J_ξ[i,j] = ((P_n[i,j,3]*fact4_jsem[2])*J2_t[j])/r_p4[i,j]   # ξ-axis
                    F_J_ζ[i,j] = (((-dP_n[i,j,3])*cos_ϕ[i,j])*J2_t[j])/r_p4[i,j]  # ζ-axis
                    # Beyond third zonal harmonic J_3,...
                    F_J_ξ_36[i,j] = zero_q_1
                    F_J_ζ_36[i,j] = zero_q_1
                    for n in 3:n1SEM[j] # min(3,n1SEM[j])
                        # ξ-axis
                        temp_fjξ[i,j,n] = (((P_n[i,j,n+1]*fact4_jsem[n])*JSEM[j,n])/temp_rn[i,j,n]) + F_J_ξ_36[i,j]
                        # ζ-axis
                        temp_fjζ[i,j,n] = ((((-dP_n[i,j,n+1])*cos_ϕ[i,j])*JSEM[j,n])/temp_rn[i,j,n]) + F_J_ζ_36[i,j]
                        F_J_ξ_36[i,j] = temp_fjξ[i,j,n]
                        F_J_ζ_36[i,j] = temp_fjζ[i,j,n]
                    end

                    # Associate Legendre functions (only for the moon)
                    if j == mo
                        for m in 1:n1SEM[mo]
                            if m == 1
                                # In this case associate Legendre functions reduce to Legendre polynomials
                                # See equations (167) and (168) in page 33 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract
                                sin_mλ[i,j,1] = sin_λ[i,j] # Moyer (1971), eq. (167)
                                cos_mλ[i,j,1] = cos_λ[i,j] # Moyer (1971), eq. (168)
                                # sec( Associate Legendre polynomial with m = n = 1 )
                                # See equation (181) in page 33 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract
                                secϕ_P_nm[i,j,1,1] = one_t
                                # Associate Legendre polynomial with m = n = 1
                                P_nm[i,j,1,1] = cos_ϕ[i,j]
                                # cosϕP_1^1'
                                # See equation (183) in page 34 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract
                                # Note: the second term equation (183) vanishes when n = m
                                cosϕ_dP_nm[i,j,1,1] = sin_ϕ[i,j]*lnm3[1]
                            else
                                # Trigonometric identity sin(λ + (m - 1)λ) and cos(λ + (m - 1)λ)
                                sin_mλ[i,j,m] = (cos_mλ[i,j,m-1]*sin_mλ[i,j,1]) + (sin_mλ[i,j,m-1]*cos_mλ[i,j,1])
                                cos_mλ[i,j,m] = (cos_mλ[i,j,m-1]*cos_mλ[i,j,1]) - (sin_mλ[i,j,m-1]*sin_mλ[i,j,1])
                                # secϕ P_n^n
                                # See equation (180) in page 33 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract
                                secϕ_P_nm[i,j,m,m] = (secϕ_P_nm[i,j,m-1,m-1]*cos_ϕ[i,j])*lnm5[m]
                                # Associate Legendre polynomial with n = m
                                P_nm[i,j,m,m] = secϕ_P_nm[i,j,m,m]*cos_ϕ[i,j]
                                # cosϕP_m^m'
                                # See equation (183) in page 34 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract
                                # Note: the second term equation (183) vanishes when n = m
                                cosϕ_dP_nm[i,j,m,m] = (secϕ_P_nm[i,j,m,m]*sin_ϕ[i,j])*lnm3[m]
                            end
                            for n in m+1:n1SEM[mo]
                                # secϕ P_n^m
                                # See equation (182) in page 34 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract
                                if n == m+1
                                    secϕ_P_nm[i,j,n,m] = (secϕ_P_nm[i,j,n-1,m]*sin_ϕ[i,j])*lnm1[n,m]
                                else
                                    secϕ_P_nm[i,j,n,m] = ((secϕ_P_nm[i,j,n-1,m]*sin_ϕ[i,j])*lnm1[n,m]) + (secϕ_P_nm[i,j,n-2,m]*lnm2[n,m])
                                end
                                # Associate Legendre polynomial of degree n and order m
                                P_nm[i,j,n,m] = secϕ_P_nm[i,j,n,m]*cos_ϕ[i,j]
                                # secϕ P_n^m
                                # See equation (183) in page 34 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract
                                cosϕ_dP_nm[i,j,n,m] = ((secϕ_P_nm[i,j,n,m]*sin_ϕ[i,j])*lnm3[n]) + (secϕ_P_nm[i,j,n-1,m]*lnm4[n,m])
                            end
                        end

                        # Moon: Compute accelerations due to tesseral harmonics C_{nm}, S_{nm}
                        # See equation (28) in page 13 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract
                        # and equation (174) in page 33 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract

                        # Accelerations due to lunar tesseral harmonics C_{21} and S_{21}
                        F_CS_ξ[i,j] = (   (  (P_nm[i,j,2,1]*lnm6[2]     )*( (C21M_t*cos_mλ[i,j,1]) + (S21M_t*sin_mλ[i,j,1]) )  ) + (  (P_nm[i,j,2,2]*lnm6[2]     )*( (C22M_t*cos_mλ[i,j,2]) + (S22M_t*sin_mλ[i,j,2]) )  )   )/r_p4[i,j]
                        F_CS_η[i,j] = (   (  (secϕ_P_nm[i,j,2,1]*lnm7[1])*( (S21M_t*cos_mλ[i,j,1]) - (C21M_t*sin_mλ[i,j,1]) )  ) + (  (secϕ_P_nm[i,j,2,2]*lnm7[2])*( (S22M_t*cos_mλ[i,j,2]) - (C22M_t*sin_mλ[i,j,2]) )  )   )/r_p4[i,j]
                        F_CS_ζ[i,j] = (   (  (cosϕ_dP_nm[i,j,2,1]       )*( (C21M_t*cos_mλ[i,j,1]) + (S21M_t*sin_mλ[i,j,1]) )  ) + (  (cosϕ_dP_nm[i,j,2,2]       )*( (C22M_t*cos_mλ[i,j,2]) + (S22M_t*sin_mλ[i,j,2]) )  )   )/r_p4[i,j]
                        # Accelerations due to lunar tesseral harmonics beyond C_{21} and S_{21}
                        F_CS_ξ_36[i,j] = zero_q_1
                        F_CS_η_36[i,j] = zero_q_1
                        F_CS_ζ_36[i,j] = zero_q_1
                        for n in 3:n2M
                            for m in 1:n
                                # Lunar teseral harmonics C_{nm}/S_{nm} * trigonometric function of integer times the longitude λ
                                Cnm_cosmλ[i,j,n,m] = CM[n,m]*cos_mλ[i,j,m]
                                Cnm_sinmλ[i,j,n,m] = CM[n,m]*sin_mλ[i,j,m]
                                Snm_cosmλ[i,j,n,m] = SM[n,m]*cos_mλ[i,j,m]
                                Snm_sinmλ[i,j,n,m] = SM[n,m]*sin_mλ[i,j,m]
                                # Vector sum in equation (173)
                                temp_CS_ξ[i,j,n,m] = (   (  (P_nm[i,j,n,m]*lnm6[n]     )*( Cnm_cosmλ[i,j,n,m] + Snm_sinmλ[i,j,n,m] )  )/temp_rn[i,j,n]   ) + F_CS_ξ_36[i,j]
                                temp_CS_η[i,j,n,m] = (   (  (secϕ_P_nm[i,j,n,m]*lnm7[m])*( Snm_cosmλ[i,j,n,m] - Cnm_sinmλ[i,j,n,m] )  )/temp_rn[i,j,n]   ) + F_CS_η_36[i,j]
                                temp_CS_ζ[i,j,n,m] = (   (  (cosϕ_dP_nm[i,j,n,m]       )*( Cnm_cosmλ[i,j,n,m] + Snm_sinmλ[i,j,n,m] )  )/temp_rn[i,j,n]   ) + F_CS_ζ_36[i,j]
                                F_CS_ξ_36[i,j] = temp_CS_ξ[i,j,n,m]
                                F_CS_η_36[i,j] = temp_CS_η[i,j,n,m]
                                F_CS_ζ_36[i,j] = temp_CS_ζ[i,j,n,m]
                            end
                        end
                        # Sum the zonal and tesseral (only for the moon) accelerations without mass parameter
                        F_JCS_ξ[i,j] = (F_J_ξ[i,j] + F_J_ξ_36[i,j]) + (F_CS_ξ[i,j]+F_CS_ξ_36[i,j])
                        F_JCS_η[i,j] = (F_CS_η[i,j]+F_CS_η_36[i,j])
                        F_JCS_ζ[i,j] = (F_J_ζ[i,j] + F_J_ζ_36[i,j]) + (F_CS_ζ[i,j]+F_CS_ζ_36[i,j])
                    else
                        # Sum the zonal accelerations without mass parameter
                        F_JCS_ξ[i,j] = (F_J_ξ[i,j] + F_J_ξ_36[i,j])
                        F_JCS_η[i,j] = zero_q_1
                        F_JCS_ζ[i,j] = (F_J_ζ[i,j] + F_J_ζ_36[i,j])
                    end

                    # R matrix: body-fixed -> "primed" (ξ, η, ζ) system
                    # See equation (161) in page 32 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract
                    Rb2p[i,j,1,1] = cos_ϕ[i,j]*cos_λ[i,j]
                    Rb2p[i,j,2,1] = -sin_λ[i,j]
                    Rb2p[i,j,3,1] = -sin_ϕ[i,j]*cos_λ[i,j]
                    Rb2p[i,j,1,2] = cos_ϕ[i,j]*sin_λ[i,j]
                    Rb2p[i,j,2,2] = cos_λ[i,j]
                    Rb2p[i,j,3,2] = -sin_ϕ[i,j]*sin_λ[i,j]
                    Rb2p[i,j,1,3] = sin_ϕ[i,j]
                    Rb2p[i,j,2,3] = zero_q_1
                    Rb2p[i,j,3,3] = cos_ϕ[i,j]
                    # G matrix: space-fixed -> body-fixed -> "primed" (ξ, η, ζ) system
                    # G_{i,j} = \sum_k R_{i,k} RotM{k,j} or G = RotM * R
                    # See equation (163) in page 32 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract
                    Gc2p[i,j,1,1] = ((Rb2p[i,j,1,1]*RotM[1,1,j]) + (Rb2p[i,j,1,2]*RotM[2,1,j])) + (Rb2p[i,j,1,3]*RotM[3,1,j])
                    Gc2p[i,j,2,1] = ((Rb2p[i,j,2,1]*RotM[1,1,j]) + (Rb2p[i,j,2,2]*RotM[2,1,j])) + (Rb2p[i,j,2,3]*RotM[3,1,j])
                    Gc2p[i,j,3,1] = ((Rb2p[i,j,3,1]*RotM[1,1,j]) + (Rb2p[i,j,3,2]*RotM[2,1,j])) + (Rb2p[i,j,3,3]*RotM[3,1,j])
                    Gc2p[i,j,1,2] = ((Rb2p[i,j,1,1]*RotM[1,2,j]) + (Rb2p[i,j,1,2]*RotM[2,2,j])) + (Rb2p[i,j,1,3]*RotM[3,2,j])
                    Gc2p[i,j,2,2] = ((Rb2p[i,j,2,1]*RotM[1,2,j]) + (Rb2p[i,j,2,2]*RotM[2,2,j])) + (Rb2p[i,j,2,3]*RotM[3,2,j])
                    Gc2p[i,j,3,2] = ((Rb2p[i,j,3,1]*RotM[1,2,j]) + (Rb2p[i,j,3,2]*RotM[2,2,j])) + (Rb2p[i,j,3,3]*RotM[3,2,j])
                    Gc2p[i,j,1,3] = ((Rb2p[i,j,1,1]*RotM[1,3,j]) + (Rb2p[i,j,1,2]*RotM[2,3,j])) + (Rb2p[i,j,1,3]*RotM[3,3,j])
                    Gc2p[i,j,2,3] = ((Rb2p[i,j,2,1]*RotM[1,3,j]) + (Rb2p[i,j,2,2]*RotM[2,3,j])) + (Rb2p[i,j,2,3]*RotM[3,3,j])
                    Gc2p[i,j,3,3] = ((Rb2p[i,j,3,1]*RotM[1,3,j]) + (Rb2p[i,j,3,2]*RotM[2,3,j])) + (Rb2p[i,j,3,3]*RotM[3,3,j])
                    # Compute cartesian coordinates of acceleration due to body figure in inertial frame (without mass parameter)
                    # See equation (169) in page 33 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract
                    F_JCS_x[i,j] = ((F_JCS_ξ[i,j]*Gc2p[i,j,1,1]) + (F_JCS_η[i,j]*Gc2p[i,j,2,1])) + (F_JCS_ζ[i,j]*Gc2p[i,j,3,1])
                    F_JCS_y[i,j] = ((F_JCS_ξ[i,j]*Gc2p[i,j,1,2]) + (F_JCS_η[i,j]*Gc2p[i,j,2,2])) + (F_JCS_ζ[i,j]*Gc2p[i,j,3,2])
                    F_JCS_z[i,j] = ((F_JCS_ξ[i,j]*Gc2p[i,j,1,3]) + (F_JCS_η[i,j]*Gc2p[i,j,2,3])) + (F_JCS_ζ[i,j]*Gc2p[i,j,3,3])
                end #if UJ_interaction[i,j]
            end # else (i != j)
        end #for i in 1:N_ext
    end #for j in 1:N_ext

    for j in 1:N_ext
        for i in 1:N_ext
            # i == j && continue
            if i == j
                continue
            else
                if UJ_interaction[i,j]
                    # Extended body accelerations
                    # J_n, C_{nm}, S_{nm} accelerations, if j-th body is flattened

                    # Add result to total acceleration upon j-th body figure due to i-th point mass
                    temp_accX_j[i,j] = accX[j] - (μ[i]*F_JCS_x[i,j])
                    accX[j] = temp_accX_j[i,j]
                    temp_accY_j[i,j] = accY[j] - (μ[i]*F_JCS_y[i,j])
                    accY[j] = temp_accY_j[i,j]
                    temp_accZ_j[i,j] = accZ[j] - (μ[i]*F_JCS_z[i,j])
                    accZ[j] = temp_accZ_j[i,j]

                    # Reaction force on i-th body
                    temp_accX_i[i,j] = accX[i] + (μ[j]*F_JCS_x[i,j])
                    accX[i] = temp_accX_i[i,j]
                    temp_accY_i[i,j] = accY[i] + (μ[j]*F_JCS_y[i,j])
                    accY[i] = temp_accY_i[i,j]
                    temp_accZ_i[i,j] = accZ[i] + (μ[j]*F_JCS_z[i,j])
                    accZ[i] = temp_accZ_i[i,j]

                    # Lunar torques
                    if j == mo
                        # Compute torques acting upon the body-figure of the Moon due to external point masses
                        # See equation (43) in page 13 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract
                        N_MfigM_pmA_x[i] = μ[i]*( (Y[i,j]*F_JCS_z[i,j]) - (Z[i,j]*F_JCS_y[i,j]) )
                        N_MfigM_pmA_y[i] = μ[i]*( (Z[i,j]*F_JCS_x[i,j]) - (X[i,j]*F_JCS_z[i,j]) )
                        N_MfigM_pmA_z[i] = μ[i]*( (X[i,j]*F_JCS_y[i,j]) - (Y[i,j]*F_JCS_x[i,j]) )
                        # Expressions below have minus sign since N_MfigM_pmA_{x,y,z} have inverted signs in cross product
                        temp_N_M_x[i] = N_MfigM[1] - N_MfigM_pmA_x[i]
                        N_MfigM[1] = temp_N_M_x[i]
                        temp_N_M_y[i] = N_MfigM[2] - N_MfigM_pmA_y[i]
                        N_MfigM[2] = temp_N_M_y[i]
                        temp_N_M_z[i] = N_MfigM[3] - N_MfigM_pmA_z[i]
                        N_MfigM[3] = temp_N_M_z[i]
                    end
                end
            end # else (i != j)
        end
    end

    #=
    Post-Newtonian corrections to gravitational acceleration
    Post-Newtonian iterative procedure setup and initialization
    See equation (35) in page 7 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract
    =#
    Threads.@threads for j in 1:N
        for i in 1:N
            # i == j && continue
            if i == j
                continue
            else
                # 4*\sum term inside {}
                _4ϕj[i,j] = 4newtonianNb_Potential[j]
                # 4*\sum + \sum terms inside {}
                ϕi_plus_4ϕj[i,j] = newtonianNb_Potential[i] + _4ϕj[i,j]
                # 2 * ||\mathbf{v_i}||^2
                _2v2[i,j] = 2v2[i]
                # \dot{s}_j^2 + 2\dot{s}_i^2 inside {}
                sj2_plus_2si2[i,j] = v2[j] + _2v2[i,j]
                # \dot{s}_j^2 + 2\dot{s}_i^2 - 4<, > terms inside {}
                sj2_plus_2si2_minus_4vivj[i,j] = sj2_plus_2si2[i,j] - (4vi_dot_vj[i,j])
                # -4\sum - \sum + \dot{s}_j^2 + 2\dot{s}_i^2  - 4<, > terms inside {}
                ϕs_and_vs[i,j] = sj2_plus_2si2_minus_4vivj[i,j] - ϕi_plus_4ϕj[i,j]
                # (\mathbf{r}_i - \mathbf{r}_j)\cdot\mathbf{v_i}
                Xij_t_Ui = X[i,j]*dq[3i-2]
                Yij_t_Vi = Y[i,j]*dq[3i-1]
                Zij_t_Wi = Z[i,j]*dq[3i]
                Rij_dot_Vi = ( Xij_t_Ui+Yij_t_Vi ) + Zij_t_Wi
                # The expression below inside the (...)^2 should have a minus sign in front of the numerator,
                # but upon squaring it is eliminated, so at the end of the day, it is irrelevant ;)
                # (\mathbf{r}_i - \mathbf{r}_j)\cdot\mathbf{v_i} / r_{ij}
                rij_dot_vi_div_rij_sq[i,j] = (Rij_dot_Vi^2)/r_p2[i,j]
                # Everything inside the {} except for the first and last terms
                pn1t2_7 = ϕs_and_vs[i,j] - (1.5rij_dot_vi_div_rij_sq[i,j])
                # Everything inside the {} except for the last term
                pn1t1_7[i,j] = c_p2+pn1t2_7
            end # else (i != j)
        end
        # Temporary post-Newtonian accelerations
        pntempX[j] = zero_q_1   # X-axis component
        pntempY[j] = zero_q_1   # Y-axis component
        pntempZ[j] = zero_q_1   # Z-axis component
    end

    Threads.@threads for j in 1:N
        for i in 1:N
            # i == j && continue
            if i == j
                continue
            else

                # First term of equation (35) in page 7 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract

                # Last term inside the {}
                pNX_t_X[i,j] = newtonX[i]*X[i,j]   # X-axis component
                pNY_t_Y[i,j] = newtonY[i]*Y[i,j]   # Y-axis component
                pNZ_t_Z[i,j] = newtonZ[i]*Z[i,j]   # Z-axis component
                # Everything inside the {} in the first term
                pn1[i,j] = (  pn1t1_7[i,j]  +  0.5*( (pNX_t_X[i,j]+pNY_t_Y[i,j]) + pNZ_t_Z[i,j] )  )
                # Full first term
                X_t_pn1[i,j] = newton_acc_X[i,j]*pn1[i,j]   # X-axis component
                Y_t_pn1[i,j] = newton_acc_Y[i,j]*pn1[i,j]   # Y-axis component
                Z_t_pn1[i,j] = newton_acc_Z[i,j]*pn1[i,j]   # Z-axis component

                # Full third term of equation (35) in page 7 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract
                pNX_t_pn3[i,j] = newtonX[i]*pn3[i,j]   # X-axis component
                pNY_t_pn3[i,j] = newtonY[i]*pn3[i,j]   # Y-axis component
                pNZ_t_pn3[i,j] = newtonZ[i]*pn3[i,j]   # Z-axis component

                # Temporary post-Newtonian accelerations
                termpnx = ( X_t_pn1[i,j] + (U_t_pn2[i,j]+pNX_t_pn3[i,j]) )   # X-axis component
                sumpnx = pntempX[j] + termpnx
                pntempX[j] = sumpnx
                termpny = ( Y_t_pn1[i,j] + (V_t_pn2[i,j]+pNY_t_pn3[i,j]) )   # Y-axis component
                sumpny = pntempY[j] + termpny
                pntempY[j] = sumpny
                termpnz = ( Z_t_pn1[i,j] + (W_t_pn2[i,j]+pNZ_t_pn3[i,j]) )   # Z-axis component
                sumpnz = pntempZ[j] + termpnz
                pntempZ[j] = sumpnz
            end # else (i != j)
        end
        # Post-Newtonian acelerations
        postNewtonX[j] = pntempX[j]*c_m2
        postNewtonY[j] = pntempY[j]*c_m2
        postNewtonZ[j] = pntempZ[j]*c_m2
    end

    #=
    Compute tidal acceleration of the Moon due to tides raised on Earth by the Sun and Moon
    See equation (32) in page 14 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract
    X_bf[mo, ea] are geocentric, Earth-fixed position coordinates in "unprimed" (inertial) frame of perturbed body (Moon)
    =#

    # Long-period tidal accelerations
    # See first term inside the {} in equation (32)

    # Position vector of the Moon in geocentric Earth true-equator-of-date frame
    x0s_M = r_star_M_0[1]
    y0s_M = r_star_M_0[2]
    z0s_M = r_star_M_0[3]

    # Cylindrical coordinates
    ρ0s2_M = (x0s_M^2) + (y0s_M^2)    # Radial (cylindrical) coordinate squared
    ρ0s_M = sqrt(ρ0s2_M)              # Radial (cylindrical) coordinate
    z0s2_M = z0s_M^2                  # z coordinate squared
    r0s2_M = ρ0s2_M + z0s2_M          # Radial (spherical) coordinate squared
    r0s_M = sqrt(r0s2_M)              # Radial (spherical) coordinate
    r0s5_M = r0s_M^5                  # Radial (spherical) coordinate ^5

    # Position vector of the Sun in geocentric Earth true-equator-of-date frame
    x0s_S = r_star_S_0[1]
    y0s_S = r_star_S_0[2]
    z0s_S = r_star_S_0[3]

    # Cylindrical coordinates
    ρ0s2_S = (x0s_S^2) + (y0s_S^2)    # Radial (cylindrical) coordinate squared
    ρ0s_S = sqrt(ρ0s2_S)              # Radial (cylindrical) coordinate
    z0s2_S = z0s_S^2                  # z coordinate squared
    r0s2_S = ρ0s2_S + z0s2_S          # Radial (spherical) coordinate squared
    r0s_S = sqrt(r0s2_S)              # Radial (spherical) coordinate
    r0s5_S = r0s_S^5                  # Radial (spherical) coordinate ^5

    # Moon and Sun coefficients
    # See second and third terms inside the first () inside the {} in equation (32)
    coeff0_M = r0s2_M - 5( ( ((Z_bf[mo,ea]*r_star_M_0[3])^2) + 0.5((r_xy[mo,ea]*ρ0s_M)^2) )/r_p2[mo,ea] )
    coeff0_S = r0s2_S - 5( ( ((Z_bf[mo,ea]*r_star_S_0[3])^2) + 0.5((r_xy[mo,ea]*ρ0s_S)^2) )/r_p2[mo,ea] )

    # Love number (2, 0) / distance^5
    # See first term inside the {} in equation (32)
    k_20E_div_r0s5_M = k_20E/r0s5_M
    k_20E_div_r0s5_S = k_20E/r0s5_S

    # Long-period tidal accelerations

    # Moon
    a_tid_0_M_x = k_20E_div_r0s5_M*((ρ0s2_M + coeff0_M)*X_bf[mo,ea])
    a_tid_0_M_y = k_20E_div_r0s5_M*((ρ0s2_M + coeff0_M)*Y_bf[mo,ea])
    a_tid_0_M_z = k_20E_div_r0s5_M*(((2z0s2_M) + coeff0_M)*Z_bf[mo,ea])
    # Sun
    a_tid_0_S_x = k_20E_div_r0s5_S*((ρ0s2_S + coeff0_S)*X_bf[mo,ea])
    a_tid_0_S_y = k_20E_div_r0s5_S*((ρ0s2_S + coeff0_S)*Y_bf[mo,ea])
    a_tid_0_S_z = k_20E_div_r0s5_S*(((2z0s2_S) + coeff0_S)*Z_bf[mo,ea])

    # Diurnal tidal accelerations
    # See first term inside the {} in equation (32)

    # Position vector of the Moon in geocentric Earth true-equator-of-date frame
    x1s_M = r_star_M_1[1]
    y1s_M = r_star_M_1[2]
    z1s_M = r_star_M_1[3]

    # Cylindrical coordinates
    ρ1s2_M = (x1s_M^2) + (y1s_M^2)    # Radial (cylindrical) coordinate squared
    ρ1s_M = sqrt(ρ1s2_M)              # Radial (cylindrical) coordinate
    z1s2_M = z1s_M^2                  # z coordinate squared
    r1s2_M = ρ1s2_M + z1s2_M          # Radial (spherical) coordinate squared
    r1s_M = sqrt(r1s2_M)              # Radial (spherical) coordinate
    r1s5_M = r1s_M^5                  # Radial (spherical) coordinate ^5

    # Position vector of the Sun in geocentric Earth true-equator-of-date frame
    x1s_S = r_star_S_1[1]
    y1s_S = r_star_S_1[2]
    z1s_S = r_star_S_1[3]

    # Cylindrical coordinates
    ρ1s2_S = (x1s_S^2) + (y1s_S^2)    # Radial (cylindrical) coordinate squared
    ρ1s_S = sqrt(ρ1s2_S)              # Radial (cylindrical) coordinate
    z1s2_S = z1s_S^2                  # z coordinate squared
    r1s2_S = ρ1s2_S + z1s2_S          # Radial (spherical) coordinate squared
    r1s_S = sqrt(r1s2_S)              # Radial (spherical) coordinate
    r1s5_S = r1s_S^5                  # Radial (spherical) coordinate ^5

    # Moon and Sun coefficients
    # See second term inside the {} in equation (32)
    coeff1_1_M = (X_bf[mo,ea]*r_star_M_1[1]) + (Y_bf[mo,ea]*r_star_M_1[2])
    coeff1_1_S = (X_bf[mo,ea]*r_star_S_1[1]) + (Y_bf[mo,ea]*r_star_S_1[2])

    coeff2_1_M = Z_bf[mo,ea]*r_star_M_1[3]
    coeff2_1_S = Z_bf[mo,ea]*r_star_S_1[3]

    coeff3_1_M = ((10coeff1_1_M)*coeff2_1_M)/r_p2[mo,ea]
    coeff3_1_S = ((10coeff1_1_S)*coeff2_1_S)/r_p2[mo,ea]

    # Love number (2, 1) / distance^5
    # See second term inside the {} in equation (32)
    k_21E_div_r1s5_M = k_21E/r1s5_M
    k_21E_div_r1s5_S = k_21E/r1s5_S

    # Diurnal tidal accelerations

    # Moon
    a_tid_1_M_x = k_21E_div_r1s5_M*((2coeff2_1_M*r_star_M_1[1]) - (coeff3_1_M*X_bf[mo,ea]))
    a_tid_1_M_y = k_21E_div_r1s5_M*((2coeff2_1_M*r_star_M_1[2]) - (coeff3_1_M*Y_bf[mo,ea]))
    a_tid_1_M_z = k_21E_div_r1s5_M*((2coeff1_1_M*r_star_M_1[3]) - (coeff3_1_M*Z_bf[mo,ea]))
    # Sun
    a_tid_1_S_x = k_21E_div_r1s5_S*((2coeff2_1_S*r_star_S_1[1]) - (coeff3_1_S*X_bf[mo,ea]))
    a_tid_1_S_y = k_21E_div_r1s5_S*((2coeff2_1_S*r_star_S_1[2]) - (coeff3_1_S*Y_bf[mo,ea]))
    a_tid_1_S_z = k_21E_div_r1s5_S*((2coeff1_1_S*r_star_S_1[3]) - (coeff3_1_S*Z_bf[mo,ea]))

    # Semi-diurnal tidal accelerations
    # See third term inside the {} in equation (32)

    # Position vector of the Moon in geocentric Earth true-equator-of-date frame
    x2s_M = r_star_M_2[1]
    y2s_M = r_star_M_2[2]
    z2s_M = r_star_M_2[3]


    ρ2s2_M = (x2s_M^2) + (y2s_M^2)    # Radial (cylindrical) coordinate squared
    ρ2s_M = sqrt(ρ2s2_M)              # Radial (cylindrical) coordinate
    z2s2_M = z2s_M^2                  # z coordinate squared
    r2s2_M = ρ2s2_M + z2s2_M          # Radial (spherical) coordinate squared
    r2s_M = sqrt(r2s2_M)              # Radial (spherical) coordinate
    r2s5_M = r2s_M^5                  # Radial (spherical) coordinate ^5

    # Position vector of the Sun in geocentric Earth true-equator-of-date frame
    x2s_S = r_star_S_2[1]
    y2s_S = r_star_S_2[2]
    z2s_S = r_star_S_2[3]

    # Cylindrical coordinates
    ρ2s2_S = (x2s_S^2) + (y2s_S^2)    # Radial (cylindrical) coordinate squared
    ρ2s_S = sqrt(ρ2s2_S)              # Radial (cylindrical) coordinate
    z2s2_S = z2s_S^2                  # z coordinate squared
    r2s2_S = ρ2s2_S + z2s2_S          # Radial (spherical) coordinate squared
    r2s_S = sqrt(r2s2_S)              # Radial (spherical) coordinate
    r2s5_S = r2s_S^5                  # Radial (spherical) coordinate ^5

    # Moon and Sun coefficients
    # See third term inside the {} in equation (32)
    coeff1_2_M = (X_bf[mo,ea]*r_star_M_2[1]) + (Y_bf[mo,ea]*r_star_M_2[2])
    coeff1_2_S = (X_bf[mo,ea]*r_star_S_2[1]) + (Y_bf[mo,ea]*r_star_S_2[2])

    coeff3_2_M = 5( (coeff1_2_M^2) - 0.5((r_xy[mo,ea]^2)*ρ2s2_M) )/r_p2[mo,ea]
    coeff3_2_S = 5( (coeff1_2_S^2) - 0.5((r_xy[mo,ea]^2)*ρ2s2_S) )/r_p2[mo,ea]

    # Love number (2, 2) / distance^5
    # See third term inside the {} in equation (32)
    k_22E_div_r2s5_M = k_22E/r2s5_M
    k_22E_div_r2s5_S = k_22E/r2s5_S

    # Semi-diurnal tidal accelerations

    # Moon
    a_tid_2_M_x = k_22E_div_r2s5_M*(  (2coeff1_2_M*r_star_M_2[1])-(ρ2s2_M+coeff3_2_M)*X_bf[mo,ea]  )
    a_tid_2_M_y = k_22E_div_r2s5_M*(  (2coeff1_2_M*r_star_M_2[2])-(ρ2s2_M+coeff3_2_M)*Y_bf[mo,ea]  )
    a_tid_2_M_z = k_22E_div_r2s5_M*(  -coeff3_2_M*Z_bf[mo,ea]  )
    # Sun
    a_tid_2_S_x = k_22E_div_r2s5_S*(  (2coeff1_2_S*r_star_S_2[1])-(ρ2s2_S+coeff3_2_S)*X_bf[mo,ea]  )
    a_tid_2_S_y = k_22E_div_r2s5_S*(  (2coeff1_2_S*r_star_S_2[2])-(ρ2s2_S+coeff3_2_S)*Y_bf[mo,ea]  )
    a_tid_2_S_z = k_22E_div_r2s5_S*(  -coeff3_2_S*Z_bf[mo,ea]  )

    # Numerical factors
    RE_div_r_p5 = (RE_au/r_p1d2[mo,ea])^5   # (Earth's radius / distance)^5
    aux_tidacc = tid_num_coeff*RE_div_r_p5  # 3*(m_E + m_M) * R_E^5 / 2 / m_E / r^5

    # Moon's and Sun's contributions * mass parameter
    a_tidal_coeff_M = μ[mo]*aux_tidacc
    a_tidal_coeff_S = μ[su]*aux_tidacc

    # Add contributions from long-period, diurnal and semi-diurnal tides (true-of-date coordinates)
    a_tidal_tod_x = (a_tidal_coeff_M*((a_tid_0_M_x+a_tid_1_M_x)+a_tid_2_M_x)) + (a_tidal_coeff_S*((a_tid_0_S_x+a_tid_1_S_x)+a_tid_2_S_x))
    a_tidal_tod_y = (a_tidal_coeff_M*((a_tid_0_M_y+a_tid_1_M_y)+a_tid_2_M_y)) + (a_tidal_coeff_S*((a_tid_0_S_y+a_tid_1_S_y)+a_tid_2_S_y))
    a_tidal_tod_z = (a_tidal_coeff_M*((a_tid_0_M_z+a_tid_1_M_z)+a_tid_2_M_z)) + (a_tidal_coeff_S*((a_tid_0_S_z+a_tid_1_S_z)+a_tid_2_S_z))

    # Transform from geocentric Earth-true-equator-of-date coordinates to geocentric mean equator of J2000.0 coordinates
    a_tidal_x = ((RotM[1,1,ea]*a_tidal_tod_x)+(RotM[2,1,ea]*a_tidal_tod_y)) + (RotM[3,1,ea]*a_tidal_tod_z)
    a_tidal_y = ((RotM[1,2,ea]*a_tidal_tod_x)+(RotM[2,2,ea]*a_tidal_tod_y)) + (RotM[3,2,ea]*a_tidal_tod_z)
    a_tidal_z = ((RotM[1,3,ea]*a_tidal_tod_x)+(RotM[2,3,ea]*a_tidal_tod_y)) + (RotM[3,3,ea]*a_tidal_tod_z)

    # Add tidal acceleration to Moon's acceleration due to extended-body effects
    accX_mo_tides = accX[mo] + a_tidal_x
    accY_mo_tides = accY[mo] + a_tidal_y
    accZ_mo_tides = accZ[mo] + a_tidal_z
    # Total extended-body effects
    accX[mo] = accX_mo_tides
    accY[mo] = accY_mo_tides
    accZ[mo] = accZ_mo_tides

    # Fill accelerations (post-Newtonian and extended body accelerations)
    # postNewton -> post-Newtonian accelerations (all bodies)
    # accX/Y/Z -> extended body accelerations (only first N_ext bodies)
    Threads.@threads for i in 1:N_ext
        dq[3(N+i)-2] = postNewtonX[i] + accX[i]
        dq[3(N+i)-1] = postNewtonY[i] + accY[i]
        dq[3(N+i)  ] = postNewtonZ[i] + accZ[i]
    end
    Threads.@threads for i in N_ext+1:N
        dq[3(N+i)-2] = postNewtonX[i]
        dq[3(N+i)-1] = postNewtonY[i]
        dq[3(N+i)  ] = postNewtonZ[i]
    end

    #=
    Lunar physical librations
    See equations (33)-(35) in pages 15-16 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract
    =#

    # Lunar moment of intertia I times angular velocity ω: Iω
    Iω_x = (I_m_t[1,1]*q[6N+4]) + ((I_m_t[1,2]*q[6N+5]) + (I_m_t[1,3]*q[6N+6])) # x-axis component
    Iω_y = (I_m_t[2,1]*q[6N+4]) + ((I_m_t[2,2]*q[6N+5]) + (I_m_t[2,3]*q[6N+6])) # y-axis component
    Iω_z = (I_m_t[3,1]*q[6N+4]) + ((I_m_t[3,2]*q[6N+5]) + (I_m_t[3,3]*q[6N+6])) # z-axis component

    # Cross product of angular velocity and Iω: ω × (I*ω)
    ωxIω_x = (q[6N+5]*Iω_z) - (q[6N+6]*Iω_y)   # x-axis component
    ωxIω_y = (q[6N+6]*Iω_x) - (q[6N+4]*Iω_z)   # y-axis component
    ωxIω_z = (q[6N+4]*Iω_y) - (q[6N+5]*Iω_x)   # z-axis component

    # Time derivative of moment of inertia times angular velocity: (dI/dt)*ω
    dIω_x = (dI_m_t[1,1]*q[6N+4]) + ((dI_m_t[1,2]*q[6N+5]) + (dI_m_t[1,3]*q[6N+6])) # x-axis component
    dIω_y = (dI_m_t[2,1]*q[6N+4]) + ((dI_m_t[2,2]*q[6N+5]) + (dI_m_t[2,3]*q[6N+6])) # y-axis component
    dIω_z = (dI_m_t[3,1]*q[6N+4]) + ((dI_m_t[3,2]*q[6N+5]) + (dI_m_t[3,3]*q[6N+6])) # z-axis component

    # Moon -> Earth radial unit vector (inertial coordinates)
    er_EM_I_1 = X[ea,mo]/r_p1d2[ea,mo]
    er_EM_I_2 = Y[ea,mo]/r_p1d2[ea,mo]
    er_EM_I_3 = Z[ea,mo]/r_p1d2[ea,mo]

    # Earth pole unit vector (inertial coordinates)
    p_E_I_1 = RotM[3,1,ea]
    p_E_I_2 = RotM[3,2,ea]
    p_E_I_3 = RotM[3,3,ea]

    # Transform Moon -> Earth radial unit vector (inertial coordinates) er_EM_I_i and
    # Earth pole unit vector p_E_I_i to lunar mantle frame coordinates
    er_EM_1 = (RotM[1,1,mo]*er_EM_I_1) + ((RotM[1,2,mo]*er_EM_I_2) + (RotM[1,3,mo]*er_EM_I_3))
    er_EM_2 = (RotM[2,1,mo]*er_EM_I_1) + ((RotM[2,2,mo]*er_EM_I_2) + (RotM[2,3,mo]*er_EM_I_3))
    er_EM_3 = (RotM[3,1,mo]*er_EM_I_1) + ((RotM[3,2,mo]*er_EM_I_2) + (RotM[3,3,mo]*er_EM_I_3))
    p_E_1 = (RotM[1,1,mo]*p_E_I_1) + ((RotM[1,2,mo]*p_E_I_2) + (RotM[1,3,mo]*p_E_I_3))
    p_E_2 = (RotM[2,1,mo]*p_E_I_1) + ((RotM[2,2,mo]*p_E_I_2) + (RotM[2,3,mo]*p_E_I_3))
    p_E_3 = (RotM[3,1,mo]*p_E_I_1) + ((RotM[3,2,mo]*p_E_I_2) + (RotM[3,3,mo]*p_E_I_3))

    # Evaluate equation (44) https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract
    # in lunar mantle frame coords

    # I*e_r
    I_er_EM_1 = (I_m_t[1,1]*er_EM_1) + ((I_m_t[1,2]*er_EM_2) + (I_m_t[1,3]*er_EM_3))
    I_er_EM_2 = (I_m_t[2,1]*er_EM_1) + ((I_m_t[2,2]*er_EM_2) + (I_m_t[2,3]*er_EM_3))
    I_er_EM_3 = (I_m_t[3,1]*er_EM_1) + ((I_m_t[3,2]*er_EM_2) + (I_m_t[3,3]*er_EM_3))

    # I*p_E
    I_p_E_1 = (I_m_t[1,1]*p_E_1) + ((I_m_t[1,2]*p_E_2) + (I_m_t[1,3]*p_E_3))
    I_p_E_2 = (I_m_t[2,1]*p_E_1) + ((I_m_t[2,2]*p_E_2) + (I_m_t[2,3]*p_E_3))
    I_p_E_3 = (I_m_t[3,1]*p_E_1) + ((I_m_t[3,2]*p_E_2) + (I_m_t[3,3]*p_E_3))

    # e_r × (I*e_r)
    er_EM_cross_I_er_EM_1 = (er_EM_2*I_er_EM_3) - (er_EM_3*I_er_EM_2)
    er_EM_cross_I_er_EM_2 = (er_EM_3*I_er_EM_1) - (er_EM_1*I_er_EM_3)
    er_EM_cross_I_er_EM_3 = (er_EM_1*I_er_EM_2) - (er_EM_2*I_er_EM_1)

    # e_r × (I*p_E)
    er_EM_cross_I_p_E_1 = (er_EM_2*I_p_E_3) - (er_EM_3*I_p_E_2)
    er_EM_cross_I_p_E_2 = (er_EM_3*I_p_E_1) - (er_EM_1*I_p_E_3)
    er_EM_cross_I_p_E_3 = (er_EM_1*I_p_E_2) - (er_EM_2*I_p_E_1)

    # p_E × (I*e_r)
    p_E_cross_I_er_EM_1 = (p_E_2*I_er_EM_3) - (p_E_3*I_er_EM_2)
    p_E_cross_I_er_EM_2 = (p_E_3*I_er_EM_1) - (p_E_1*I_er_EM_3)
    p_E_cross_I_er_EM_3 = (p_E_1*I_er_EM_2) - (p_E_2*I_er_EM_1)

    # p_E × (I*p_E)
    p_E_cross_I_p_E_1 = (p_E_2*I_p_E_3) - (p_E_3*I_p_E_2)
    p_E_cross_I_p_E_2 = (p_E_3*I_p_E_1) - (p_E_1*I_p_E_3)
    p_E_cross_I_p_E_3 = (p_E_1*I_p_E_2) - (p_E_2*I_p_E_1)

    # Coefficients of first and second terms inside {} in eq. (44)
    one_minus_7sin2ϕEM = one_t - (7((sin_ϕ[ea,mo])^2))
    two_sinϕEM = 2sin_ϕ[ea,mo]

    # Overall numerical factor in eq. (44) / r_{EM}^5
    N_MfigM_figE_factor_div_rEMp5 = (N_MfigM_figE_factor/(r_p1d2[mo,ea]^5))
    # Evaluation of eq. (44)
    N_MfigM_figE_1 = N_MfigM_figE_factor_div_rEMp5*( (one_minus_7sin2ϕEM*er_EM_cross_I_er_EM_1) + (two_sinϕEM*(er_EM_cross_I_p_E_1+p_E_cross_I_er_EM_1)) - (0.4p_E_cross_I_p_E_1))
    N_MfigM_figE_2 = N_MfigM_figE_factor_div_rEMp5*( (one_minus_7sin2ϕEM*er_EM_cross_I_er_EM_2) + (two_sinϕEM*(er_EM_cross_I_p_E_2+p_E_cross_I_er_EM_2)) - (0.4p_E_cross_I_p_E_2))
    N_MfigM_figE_3 = N_MfigM_figE_factor_div_rEMp5*( (one_minus_7sin2ϕEM*er_EM_cross_I_er_EM_3) + (two_sinϕEM*(er_EM_cross_I_p_E_3+p_E_cross_I_er_EM_3)) - (0.4p_E_cross_I_p_E_3))

    # RHS of equation (34) in page 16 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract

    # Torques acting upon lunar body-figure due to external point masses: transform coordinates from inertial frame to lunar mantle frame
    N_1_LMF = (RotM[1,1,mo]*N_MfigM[1]) + ((RotM[1,2,mo]*N_MfigM[2]) + (RotM[1,3,mo]*N_MfigM[3]))
    N_2_LMF = (RotM[2,1,mo]*N_MfigM[1]) + ((RotM[2,2,mo]*N_MfigM[2]) + (RotM[2,3,mo]*N_MfigM[3]))
    N_3_LMF = (RotM[3,1,mo]*N_MfigM[1]) + ((RotM[3,2,mo]*N_MfigM[2]) + (RotM[3,3,mo]*N_MfigM[3]))

    # Torque on the mantle due to the interaction between core and mantle (evaluated in mantle frame)
    # See equation (45) in page 18 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract
    N_cmb_1 = (k_ν*(q[6N+10]-q[6N+4])) - (C_c_m_A_c*(q[6N+12]*q[6N+11]))
    N_cmb_2 = (k_ν*(q[6N+11]-q[6N+5])) + (C_c_m_A_c*(q[6N+12]*q[6N+10]))
    N_cmb_3 = (k_ν*(q[6N+12]-q[6N+6]))

    # I*(dω/dt); i.e., I times RHS of equation (34) in page 16 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract
    I_dω_1 = ((N_MfigM_figE_1 + (μ[mo]*N_1_LMF)) + N_cmb_1) - (dIω_x + ωxIω_x)
    I_dω_2 = ((N_MfigM_figE_2 + (μ[mo]*N_2_LMF)) + N_cmb_2) - (dIω_y + ωxIω_y)
    I_dω_3 = ((N_MfigM_figE_3 + (μ[mo]*N_3_LMF)) + N_cmb_3) - (dIω_z + ωxIω_z)

    # RHS of equation (35) in page 16 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract

    # I_c * ω_c
    Ic_ωc_1 = I_c_t[1,1]*q[6N+10] # + ((I_c_t[1,2]*q[6N+11]) + (I_c_t[1,3]*q[6N+12]))
    Ic_ωc_2 = I_c_t[2,2]*q[6N+11] # + ((I_c_t[2,1]*q[6N+10]) + (I_c_t[2,3]*q[6N+12]))
    Ic_ωc_3 = I_c_t[3,3]*q[6N+12] # + ((I_c_t[3,1]*q[6N+10]) + (I_c_t[3,2]*q[6N+11]))

    # - ω_m × (I_c * ω_c)
    m_ωm_x_Icωc_1 = (q[6N+6]*Ic_ωc_2) - (q[6N+5]*Ic_ωc_3)
    m_ωm_x_Icωc_2 = (q[6N+4]*Ic_ωc_3) - (q[6N+6]*Ic_ωc_1)
    m_ωm_x_Icωc_3 = (q[6N+5]*Ic_ωc_1) - (q[6N+4]*Ic_ωc_2)

    # I_c*(dω_c/dt); i.e., I_c times RHS of of equation (35) in page 16 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract
    Ic_dωc_1 = m_ωm_x_Icωc_1 - N_cmb_1
    Ic_dωc_2 = m_ωm_x_Icωc_2 - N_cmb_2
    Ic_dωc_3 = m_ωm_x_Icωc_3 - N_cmb_3

    # Lunar mantle physical librations

    # Euler angles
    # See equation (14) in page 9 of  https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract
    dq[6N+1] = ((q[6N+4]*sin(q[6N+3])) + (q[6N+5]*cos(q[6N+3])) )/sin(q[6N+2])
    dq[6N+2] = (q[6N+4]*cos(q[6N+3])) - (q[6N+5]*sin(q[6N+3]))
    dq[6N+3] = q[6N+6] - (dq[6N+1]*cos(q[6N+2]))
    # Angular velocitiy
    # See equation (34) in page 16 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract
    dq[6N+4] = (inv_I_m_t[1,1]*I_dω_1) + ( (inv_I_m_t[1,2]*I_dω_2) + (inv_I_m_t[1,3]*I_dω_3) )
    dq[6N+5] = (inv_I_m_t[2,1]*I_dω_1) + ( (inv_I_m_t[2,2]*I_dω_2) + (inv_I_m_t[2,3]*I_dω_3) )
    dq[6N+6] = (inv_I_m_t[3,1]*I_dω_1) + ( (inv_I_m_t[3,2]*I_dω_2) + (inv_I_m_t[3,3]*I_dω_3) )

    # Lunar core physical librations

    # Euler angles
    # See equation (15) in page 9 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract
    # (core angular velocity components ω_c_CE_i represent lunar core-equator frame coordinates)
    dq[6N+9] = -(ω_c_CE_2/sin(q[6N+8])) ### evaluated first, since it's used below
    dq[6N+7] = ω_c_CE_3-(dq[6N+9]*cos(q[6N+8]))
    dq[6N+8] = ω_c_CE_1
    # Angular velocity of the core expressed in the mantle frame
    # See equation (35) in page 16 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract
    dq[6N+10] = inv_I_c_t[1,1]*Ic_dωc_1 # + ( (inv_I_c_t[1,2]*Ic_dωc_2) + (inv_I_c_t[1,3]*Ic_dωc_3) )
    dq[6N+11] = inv_I_c_t[2,2]*Ic_dωc_2 # + ( (inv_I_c_t[2,1]*Ic_dωc_1) + (inv_I_c_t[2,3]*Ic_dωc_3) )
    dq[6N+12] = inv_I_c_t[3,3]*Ic_dωc_3 # + ( (inv_I_c_t[3,1]*Ic_dωc_1) + (inv_I_c_t[3,2]*Ic_dωc_2) )

    # TT-TDB

    # Contribution to TT-TDB due to Solar J2
    # See equation (7) in page 7 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract
    w_LE = (newtonianCoeff[su,ea]*J2_t[su])*((one_t - (3(sin_ϕ[su,ea]^2)))/2)

    # Contributions of order 1/c^2 to TT-TDB
    # See equation (2) of https://ui.adsabs.harvard.edu/abs/2009A%26A...507.1675F/abstract
    α_TTmTDB = ((0.5v2[ea]) + newtonianNb_Potential[ea]) + w_LE

    # Contributions of order 1/c^4 to TT-TDB
    # See equation (3) of https://ui.adsabs.harvard.edu/abs/2009A%26A...507.1675F/abstract
    v4E = v2[ea]^2 # v_Earth^4
    ϕ_Earth_Newtonian_sq = newtonianNb_Potential[ea]^2
    β_TTmTDB = ( ϕ_Earth_Newtonian_sq / 2 ) - ( v4E / 8 )
    for i in 1:N
        if i == ea
            continue
        else
            β_TTmTDB_i_1 = 4vi_dot_vj[i,ea]
            β_TTmTDB_i_2 = newtonianNb_Potential[i] - ( (1.5v2[ea]) + (2v2[i]) )
            β_TTmTDB_i_3 = (  ((dq[3(N+i)-2]*X[i,ea]) + (dq[3(N+i)-1]*Y[i,ea])) + (dq[3(N+i)]*Z[i,ea])  ) /2
            β_TTmTDB_i_4 = rij_dot_vi_div_rij_sq[i,ea]/2
            β_TTmTDB_i = ((β_TTmTDB_i_1 + β_TTmTDB_i_2) + ( β_TTmTDB_i_3 + β_TTmTDB_i_4))
            temp_β_TTmTDB = β_TTmTDB + (newtonian1b_Potential[i,ea]*β_TTmTDB_i)
            β_TTmTDB = temp_β_TTmTDB
        end # else (i != j)
    end

    # d(TT-TDB)/d(TDB)
    # See equation (10) of https://ui.adsabs.harvard.edu/abs/2009A%26A...507.1675F/abstract
    # Since the initial condition is in seconds, and the independent variable is in days,
    # an additional factor daysec=86400 is needed for unit consistency
    dq[6N+13] = daysec*(( L_B - (c_m2*α_TTmTDB))*one_plus_L_B_minus_L_G + ( (c_m4*β_TTmTDB) - L_G ))

    nothing
end
