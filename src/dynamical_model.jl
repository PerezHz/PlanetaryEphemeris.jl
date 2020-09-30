macro taylorize_pleph(ex)
    nex = TaylorIntegration._make_parsed_jetcoeffs(ex)
    return quote
        $(esc(ex))
        $(esc(nex))
    end
end

# Solar System ( JPL DE430/431) dynamical model
# Bodies considered in the model are: the Sun, the eight planets, the Moon and
# the 343 main-belt asteroids included in the JPL DE 430 ephemeris
# effects considered are:
# - post-Newtonian point-mass accelerations between all bodies,
# - figure-effects (oblateness) of the Earth (J2 and J3)
# - J2 effect of the Sun
# - J2 and J3 effect of the Moon
# - Kinematic model for the precession and nutation of the Earth's orientation (IAU 1976/1980 Earth orientation model)
# - Kinematic model for the Moons's orientation (Seidelmann et al., 2006)
# - tidal secular acceleration of Moon due to rides raised on Earth by both the Moon and the Sun
function NBP_pN_A_J23E_J23M_J2S!(dq, q, params, t)
    # N: number of bodies
    # S: auxiliary variable =eltype(q0)
    # eulang_de430_: Taylor interpolant for DE430 lunar orientation Euler angles
    local N, S, eulang_de430_, jd0 = params
    local eulang_t = eulang_de430_( (t+(jd0-J2000))*daysec )
    #local eulang_t_del = eulang_de430_( ((t-τ_M)+(jd0-J2000))*daysec )

    #TODO: handle appropiately @taylorize'd version with postnewton_iter>1
    local postnewton_iter = 1 # number of iterations of post-Newtonian subroutine

    # parameters related to speed of light, c
    local c_p2 = 29979.063823897606 # c^2 = 29979.063823897606 au^2/d^2
    local c_m2 = 3.3356611996764786e-5 # c^-2 = 3.3356611996764786e-5 d^2/au^2

    local zero_q_1 = zero(q[1])
    local one_t = one(t)

    X = Array{Taylor1{S}}(undef, N, N)
    Y = Array{Taylor1{S}}(undef, N, N)
    Z = Array{Taylor1{S}}(undef, N, N)

    r_p2 = Array{Taylor1{S}}(undef, N, N)
    r_p3d2 = Array{Taylor1{S}}(undef, N, N)
    r_p7d2 = Array{Taylor1{S}}(undef, N, N)

    #Newtonian accelerations
    newtonX = Array{Taylor1{S}}(undef, N)
    newtonY = Array{Taylor1{S}}(undef, N)
    newtonZ = Array{Taylor1{S}}(undef, N)

    newtonianCoeff = Array{Taylor1{S}}(undef, N, N)

    #post-Newtonian stuff
    U = Array{Taylor1{S}}(undef, N, N)
    V = Array{Taylor1{S}}(undef, N, N)
    W = Array{Taylor1{S}}(undef, N, N)

    _4U_m_3X = Array{Taylor1{S}}(undef, N, N)
    _4V_m_3Y = Array{Taylor1{S}}(undef, N, N)
    _4W_m_3Z = Array{Taylor1{S}}(undef, N, N)

    UU = Array{Taylor1{S}}(undef, N, N)
    VV = Array{Taylor1{S}}(undef, N, N)
    WW = Array{Taylor1{S}}(undef, N, N)

    r_p1d2 = Array{Taylor1{S}}(undef, N, N)

    newtonianNb_Potential = Array{Taylor1{S}}(undef, N)
    newtonian1b_Potential = Array{Taylor1{S}}(undef, N, N)
    newtonianCoeff = Array{Taylor1{S}}(undef, N, N)
    newton_acc_X = Array{Taylor1{S}}(undef, N, N)
    newton_acc_Y = Array{Taylor1{S}}(undef, N, N)
    newton_acc_Z = Array{Taylor1{S}}(undef, N, N)

    v2 = Array{Taylor1{S}}(undef, N)
    vi_dot_vj = Array{Taylor1{S}}(undef, N, N)
    pn2 = Array{Taylor1{S}}(undef, N, N)
    pn3 = Array{Taylor1{S}}(undef, N, N)
    _4ϕj = Array{Taylor1{S}}(undef, N, N)
    ϕi_plus_4ϕj = Array{Taylor1{S}}(undef, N, N)
    sj2_plus_2si2_minus_4vivj = Array{Taylor1{S}}(undef, N, N)
    ϕs_and_vs = Array{Taylor1{S}}(undef, N, N)
    U_t_pn2 = Array{Taylor1{S}}(undef, N, N)
    V_t_pn2 = Array{Taylor1{S}}(undef, N, N)
    W_t_pn2 = Array{Taylor1{S}}(undef, N, N)
    pn1t1_7 = Array{Taylor1{S}}(undef, N, N)

    pntempX = Array{Taylor1{S}}(undef, N, postnewton_iter)
    pntempY = Array{Taylor1{S}}(undef, N, postnewton_iter)
    pntempZ = Array{Taylor1{S}}(undef, N, postnewton_iter)
    pn1 = Array{Taylor1{S}}(undef, N, N, postnewton_iter)
    X_t_pn1 = Array{Taylor1{S}}(undef, N, N, postnewton_iter)
    Y_t_pn1 = Array{Taylor1{S}}(undef, N, N, postnewton_iter)
    Z_t_pn1 = Array{Taylor1{S}}(undef, N, N, postnewton_iter)
    pNX_t_pn3 = Array{Taylor1{S}}(undef, N, N, postnewton_iter)
    pNY_t_pn3 = Array{Taylor1{S}}(undef, N, N, postnewton_iter)
    pNZ_t_pn3 = Array{Taylor1{S}}(undef, N, N, postnewton_iter)
    pNX_t_X = Array{Taylor1{S}}(undef, N, N, postnewton_iter)
    pNY_t_Y = Array{Taylor1{S}}(undef, N, N, postnewton_iter)
    pNZ_t_Z = Array{Taylor1{S}}(undef, N, N, postnewton_iter)
    postNewtonX = Array{Taylor1{S}}(undef, N, postnewton_iter+1)
    postNewtonY = Array{Taylor1{S}}(undef, N, postnewton_iter+1)
    postNewtonZ = Array{Taylor1{S}}(undef, N, postnewton_iter+1)

    # (Jn, Cmn, Smn) acceleration auxiliaries
    X_bf_1 = Array{Taylor1{S}}(undef, N, N)
    Y_bf_1 = Array{Taylor1{S}}(undef, N, N)
    Z_bf_1 = Array{Taylor1{S}}(undef, N, N)
    X_bf_2 = Array{Taylor1{S}}(undef, N, N)
    Y_bf_2 = Array{Taylor1{S}}(undef, N, N)
    Z_bf_2 = Array{Taylor1{S}}(undef, N, N)
    X_bf_3 = Array{Taylor1{S}}(undef, N, N)
    Y_bf_3 = Array{Taylor1{S}}(undef, N, N)
    Z_bf_3 = Array{Taylor1{S}}(undef, N, N)
    X_bf = Array{Taylor1{S}}(undef, N, N)
    Y_bf = Array{Taylor1{S}}(undef, N, N)
    Z_bf = Array{Taylor1{S}}(undef, N, N)
    F_JCS_x = Array{Taylor1{S}}(undef, N, N)
    F_JCS_y = Array{Taylor1{S}}(undef, N, N)
    F_JCS_z = Array{Taylor1{S}}(undef, N, N)
    temp_accX_j = Array{Taylor1{S}}(undef, N, N)
    temp_accY_j = Array{Taylor1{S}}(undef, N, N)
    temp_accZ_j = Array{Taylor1{S}}(undef, N, N)
    temp_accX_i = Array{Taylor1{S}}(undef, N, N)
    temp_accY_i = Array{Taylor1{S}}(undef, N, N)
    temp_accZ_i = Array{Taylor1{S}}(undef, N, N)
    sin_ϕ = Array{Taylor1{S}}(undef, N, N)
    cos_ϕ = Array{Taylor1{S}}(undef, N, N)
    sin_λ = Array{Taylor1{S}}(undef, N, N)
    cos_λ = Array{Taylor1{S}}(undef, N, N)
    r_xy = Array{Taylor1{S}}(undef, N, N)
    r_p4 = Array{Taylor1{S}}(undef, N, N)
    P_n = Array{Taylor1{S}}(undef, N, N, maximum(n1SEM)+1)
    dP_n = Array{Taylor1{S}}(undef, N, N, maximum(n1SEM)+1)
    temp_fjξ = Array{Taylor1{S}}(undef, N, N, maximum(n1SEM)+1)
    temp_fjζ = Array{Taylor1{S}}(undef, N, N, maximum(n1SEM)+1)
    temp_rn = Array{Taylor1{S}}(undef, N, N, maximum(n1SEM)+1)
    F_CS_ξ_36 = Array{Taylor1{S}}(undef, N, N)
    F_CS_η_36 = Array{Taylor1{S}}(undef, N, N)
    F_CS_ζ_36 = Array{Taylor1{S}}(undef, N, N)
    F_J_ξ_36 = Array{Taylor1{S}}(undef, N, N)
    F_J_ζ_36 = Array{Taylor1{S}}(undef, N, N)
    sin_mλ = Array{Taylor1{S}}(undef, N, N, n1SEM[mo])
    cos_mλ = Array{Taylor1{S}}(undef, N, N, n1SEM[mo])
    secϕ_P_nm = Array{Taylor1{S}}(undef, N, N, n1SEM[mo]+1, n1SEM[mo]+1)
    P_nm = Array{Taylor1{S}}(undef, N, N, n1SEM[mo]+1, n1SEM[mo]+1)
    cosϕ_dP_nm = Array{Taylor1{S}}(undef, N, N, n1SEM[mo]+1, n1SEM[mo]+1)
    F_J_ξ = Array{Taylor1{S}}(undef, N, N)
    F_J_η = Array{Taylor1{S}}(undef, N, N)
    F_J_ζ = Array{Taylor1{S}}(undef, N, N)
    F_CS_ξ = Array{Taylor1{S}}(undef, N, N)
    F_CS_η = Array{Taylor1{S}}(undef, N, N)
    F_CS_ζ = Array{Taylor1{S}}(undef, N, N)
    F_JCS_ξ = Array{Taylor1{S}}(undef, N, N)
    F_JCS_η = Array{Taylor1{S}}(undef, N, N)
    F_JCS_ζ = Array{Taylor1{S}}(undef, N, N)
    Rb2p = Array{Taylor1{S}}(undef, N, N, 3, 3) #R matrix body-fixed to "primed" ξηζ frame (Moyer, 1971, eq. 161)
    Gc2p = Array{Taylor1{S}}(undef, N, N, 3, 3) #G matrix "space-fixed" to "primed" ξηζ frame (Moyer, 1971, eq. 163)

    # extended-body accelerations
    accX = Array{Taylor1{S}}(undef, N)
    accY = Array{Taylor1{S}}(undef, N)
    accZ = Array{Taylor1{S}}(undef, N)

    # rotations to and from Earth, Sun and Moon pole-oriented frames
    local dsj2k = t+(jd0-2.451545e6) # J2000.0 = 2.451545e6
    local αs = deg2rad(α_p_sun*one_t)
    local δs = deg2rad(δ_p_sun*one_t)
    local αm = eulang_t[1] - (pi/2)
    local δm = (pi/2) - eulang_t[2]
    local Wm = eulang_t[3]
    local M_ = Array{Taylor1{S}}(undef, 3, 3, N)
    local M_[:,:,ea] = c2t_jpl_de430(dsj2k)
    local M_[:,:,su] = pole_rotation(αs, δs)
    local M_[:,:,mo] = pole_rotation(αm, δm, Wm)
    local ITM_t = ITM_und.*one_t #+ ITM2(eulang_t_del[1], eulang_t_del[2], eulang_t_del[3])
    local fact_num = -4.5257273867882326e-36 # == -k_2M*μ[ea]*(R_moon^5)
    local fact1_jsem = [(2n-1)/n for n in 1:maximum(n1SEM)]
    local fact2_jsem = [(n-1)/n for n in 1:maximum(n1SEM)]
    local fact3_jsem = [n for n in 1:maximum(n1SEM)]
    local fact4_jsem = [n+1 for n in 1:maximum(n1SEM)]
    local fact5_jsem = [(n+2) for n in 1:maximum(n1SEM)]
    local lnm1 = [(2n-1)/(n-m) for n in 1:6, m in 1:6]
    local lnm2 = [-(n+m-1)/(n-m) for n in 1:6, m in 1:6]
    local lnm3 = [-n for n in 1:6]
    local lnm4 = [n+m for n in 1:6, m in 1:6]
    local lnm5 = [2n-1 for n in 1:6]
    local lnm6 = [-(n+1) for n in 1:6]
    local lnm7 = [m for m in 1:6]
    # TODO: solve differences between parsed and non-parsed
    local J2E_t = (J2E + J2EDOT*(dsj2k/yr))*((RE/au)^2)
    local J2S_t = JSEM[su,2]*one_t

    for j in 1:N
        newtonX[j] = zero_q_1
        newtonY[j] = zero_q_1
        newtonZ[j] = zero_q_1

        newtonianNb_Potential[j] = zero_q_1

        accX[j] = zero_q_1
        accY[j] = zero_q_1
        accZ[j] = zero_q_1

        dq[3j-2] = q[3(N+j)-2]
        dq[3j-1] = q[3(N+j)-1]
        dq[3j  ] = q[3(N+j)  ]
    end

    #compute point-mass Newtonian accelerations, all bodies
    for j in 1:N
        for i in 1:N
            # i == j && continue
            if i == j
                continue
            else
                X[i,j] = q[3i-2]-q[3j-2]
                Y[i,j] = q[3i-1]-q[3j-1]
                Z[i,j] = q[3i]-q[3j]

                U[i,j] = dq[3i-2]-dq[3j-2]
                V[i,j] = dq[3i-1]-dq[3j-1]
                W[i,j] = dq[3i  ]-dq[3j  ]

                _4U_m_3X[i,j] = (4dq[3j-2])-(3dq[3i-2])
                _4V_m_3Y[i,j] = (4dq[3j-1])-(3dq[3i-1])
                _4W_m_3Z[i,j] = (4dq[3j  ])-(3dq[3i  ])

                pn2x = X[i,j]*_4U_m_3X[i,j]
                pn2y = Y[i,j]*_4V_m_3Y[i,j]
                pn2z = Z[i,j]*_4W_m_3Z[i,j]

                UU[i,j] = dq[3i-2]*dq[3j-2]
                VV[i,j] = dq[3i-1]*dq[3j-1]
                WW[i,j] = dq[3i  ]*dq[3j  ]

                vi_dot_vj[i,j] = ( UU[i,j]+VV[i,j] ) + WW[i,j]

                r_p2[i,j] = ( (X[i,j]^2)+(Y[i,j]^2) ) + (Z[i,j]^2)

                r_p1d2[i,j] = sqrt(r_p2[i,j])
                r_p3d2[i,j] = r_p2[i,j]^1.5
                r_p7d2[i,j] = r_p2[i,j]^3.5

                newtonianCoeff[i,j] =  μ[i]/r_p3d2[i,j]

                pn2[i,j] = newtonianCoeff[i,j]*(( pn2x+pn2y ) + pn2z)

                newton_acc_X[i,j] = X[i,j]*newtonianCoeff[i,j]
                newton_acc_Y[i,j] = Y[i,j]*newtonianCoeff[i,j]
                newton_acc_Z[i,j] = Z[i,j]*newtonianCoeff[i,j]

                newtonian1b_Potential[i,j] = μ[i]/r_p1d2[i, j]
                pn3[i,j] = 3.5newtonian1b_Potential[i,j]
                U_t_pn2[i,j] = pn2[i,j]*U[i,j]
                V_t_pn2[i,j] = pn2[i,j]*V[i,j]
                W_t_pn2[i,j] = pn2[i,j]*W[i,j]

                temp_001 = newtonX[j] + (X[i,j]*newtonianCoeff[i,j])
                newtonX[j] = temp_001
                temp_002 = newtonY[j] + (Y[i,j]*newtonianCoeff[i,j])
                newtonY[j] = temp_002
                temp_003 = newtonZ[j] + (Z[i,j]*newtonianCoeff[i,j])
                newtonZ[j] = temp_003
                temp_004 = newtonianNb_Potential[j] + newtonian1b_Potential[i, j]
                newtonianNb_Potential[j] = temp_004
            end # else (i != j)
        end #for, i
        v2[j] = ( (dq[3j-2]^2)+(dq[3j-1]^2) ) + (dq[3j]^2)
    end #for, j

    J2M_t = ( ITM_t[3,3] - ((ITM_t[1,1]+ITM_t[2,2])/2) )/(μ[mo]) # J_{2,M}*R_M^2
    C22M_t = ((ITM_t[2,2] - ITM_t[1,1])/(μ[mo]))/4 # C_{22,M}*R_M^2
    C21M_t = (-ITM_t[1,3])/(μ[mo]) # C_{21,M}*R_M^2
    S21M_t = (-ITM_t[3,2])/(μ[mo]) # S_{21,M}*R_M^2
    S22M_t = ((-ITM_t[2,1])/(μ[mo]))/2 # S_{22,M}*R_M^2

    for j in 1:N
        for i in 1:N
            # i == j && continue
            if i == j
                continue
            else
                #Jn, Cnm, Snm accelerations, if j-th body is flattened
                if UJ_interaction[i,j]
                    # # rotate from inertial frame to extended-body frame
                    X_bf_1[i,j] = X[i,j]*M_[1,1,j]
                    X_bf_2[i,j] = Y[i,j]*M_[1,2,j]
                    X_bf_3[i,j] = Z[i,j]*M_[1,3,j]
                    Y_bf_1[i,j] = X[i,j]*M_[2,1,j]
                    Y_bf_2[i,j] = Y[i,j]*M_[2,2,j]
                    Y_bf_3[i,j] = Z[i,j]*M_[2,3,j]
                    Z_bf_1[i,j] = X[i,j]*M_[3,1,j]
                    Z_bf_2[i,j] = Y[i,j]*M_[3,2,j]
                    Z_bf_3[i,j] = Z[i,j]*M_[3,3,j]
                    X_bf[i,j] = (X_bf_1[i,j] + X_bf_2[i,j]) + (X_bf_3[i,j]) # x-coordinate in body-fixed frame
                    Y_bf[i,j] = (Y_bf_1[i,j] + Y_bf_2[i,j]) + (Y_bf_3[i,j]) # y-coordinate in body-fixed frame
                    Z_bf[i,j] = (Z_bf_1[i,j] + Z_bf_2[i,j]) + (Z_bf_3[i,j]) # z-coordinate in body-fixed frame

                    # compute cartesian coordinates of acceleration due to body figure in body frame
                    sin_ϕ[i,j] = Z_bf[i,j]/r_p1d2[i,j] # Moyer (1971), eq. (165)
                    r_xy[i,j] = sqrt( (X_bf[i,j]^2)+(Y_bf[i,j]^2) )
                    cos_ϕ[i,j] = r_xy[i,j]/r_p1d2[i,j] # Moyer (1971), eq. (166)
                    sin_λ[i,j] = Y_bf[i,j]/r_xy[i,j] # Moyer (1971), eq. (167)
                    cos_λ[i,j] = X_bf[i,j]/r_xy[i,j] # Moyer (1971), eq. (168)

                    # compute accelerations due to zonal harmonics J_{n}
                    P_n[i,j,1] = one_t
                    P_n[i,j,2] = sin_ϕ[i,j]
                    dP_n[i,j,1] = zero_q_1
                    dP_n[i,j,2] = one_t
                    for n in 2:n1SEM[j] #min(3,n1SEM[j])
                        P_n[i,j,n+1] = ((P_n[i,j,n]*sin_ϕ[i,j])*fact1_jsem[n]) - (P_n[i,j,n-1]*fact2_jsem[n])
                        dP_n[i,j,n+1] = (dP_n[i,j,n]*sin_ϕ[i,j]) + (P_n[i,j,n]*fact3_jsem[n])
                        temp_rn[i,j,n] = r_p1d2[i,j]^fact5_jsem[n]
                    end
                    r_p4[i,j] = r_p2[i,j]^2
                    if j == mo
                        F_J_ξ[i,j] = ((P_n[i,j,3]*fact4_jsem[2])*J2M_t)/r_p4[i,j]
                        F_J_ζ[i,j] = (((-dP_n[i,j,3])*cos_ϕ[i,j])*J2M_t)/r_p4[i,j]
                    elseif j == ea
                        F_J_ξ[i,j] = ((P_n[i,j,3]*fact4_jsem[2])*J2E_t)/r_p4[i,j]
                        F_J_ζ[i,j] = (((-dP_n[i,j,3])*cos_ϕ[i,j])*J2E_t)/r_p4[i,j]
                    elseif j == su
                        F_J_ξ[i,j] = ((P_n[i,j,3]*fact4_jsem[2])*J2S_t)/r_p4[i,j]
                        F_J_ζ[i,j] = (((-dP_n[i,j,3])*cos_ϕ[i,j])*J2S_t)/r_p4[i,j]
                    end
                    F_J_ξ_36[i,j] = zero_q_1
                    F_J_ζ_36[i,j] = zero_q_1
                    for n in 3:n1SEM[j] #min(3,n1SEM[j])
                        temp_fjξ[i,j,n] = (((P_n[i,j,n+1]*fact4_jsem[n])*JSEM[j,n])/temp_rn[i,j,n]) + F_J_ξ_36[i,j]
                        temp_fjζ[i,j,n] = ((((-dP_n[i,j,n+1])*cos_ϕ[i,j])*JSEM[j,n])/temp_rn[i,j,n]) + F_J_ζ_36[i,j]
                        F_J_ξ_36[i,j] = temp_fjξ[i,j,n]
                        F_J_ζ_36[i,j] = temp_fjζ[i,j,n]
                    end

                    # Moon: compute accelerations due to tesseral harmonics C_{nm}, S_{nm}
                    if j == mo
                        # Associate Legendre polynomials recursion
                        for m in 1:n1SEM[mo]
                            if m == 1
                                sin_mλ[i,j,1] = sin_λ[i,j] # Moyer (1971), eq. (167)
                                cos_mλ[i,j,1] = cos_λ[i,j] # Moyer (1971), eq. (168)
                                secϕ_P_nm[i,j,1,1] = one_t
                            else
                                sin_mλ[i,j,m] = (sin_mλ[i,j,1]*cos_mλ[i,j,m-1]) + (cos_mλ[i,j,1]*sin_mλ[i,j,m-1])
                                cos_mλ[i,j,m] = (cos_mλ[i,j,1]*cos_mλ[i,j,m-1]) - (sin_mλ[i,j,1]*sin_mλ[i,j,m-1])
                                secϕ_P_nm[i,j,m,m] = (secϕ_P_nm[i,j,m-1,m-1]*cos_ϕ[i,j])*lnm5[m]
                                P_nm[i,j,m,m] = secϕ_P_nm[i,j,m,m]*cos_ϕ[i,j]
                                cosϕ_dP_nm[i,j,m,m] = (secϕ_P_nm[i,j,m,m]*sin_ϕ[i,j])*lnm3[m] #+0 (second term in Eq. 183 from Moyer, 1971, vanishes when n=m)
                            end
                            for n in m+1:n1SEM[mo]
                                if n == m+1
                                    secϕ_P_nm[i,j,n,m] = (secϕ_P_nm[i,j,n-1,m]*sin_ϕ[i,j])*lnm1[n,m]
                                else
                                    secϕ_P_nm[i,j,n,m] = ((secϕ_P_nm[i,j,n-1,m]*sin_ϕ[i,j])*lnm1[n,m]) + (secϕ_P_nm[i,j,n-2,m]*lnm2[n,m])
                                end
                                P_nm[i,j,n,m] = secϕ_P_nm[i,j,n,m]*cos_ϕ[i,j]
                                cosϕ_dP_nm[i,j,n,m] = ((secϕ_P_nm[i,j,n,m]*sin_ϕ[i,j])*lnm3[n]) + (secϕ_P_nm[i,j,n-1,m]*lnm4[n,m])
                            end
                        end

                        F_CS_ξ[i,j] = (((P_nm[i,j,2,1]*lnm6[2])*((C21M_t*cos_mλ[i,j,1])+(S21M_t*sin_mλ[i,j,1]))) + ((P_nm[i,j,2,2]*lnm6[2])*((C22M_t*cos_mλ[i,j,2])+(S22M_t*sin_mλ[i,j,2]))))/r_p4[i,j]
                        F_CS_η[i,j] = (((secϕ_P_nm[i,j,2,1]*lnm7[1])*((S21M_t*cos_mλ[i,j,1])-(C21M_t*sin_mλ[i,j,1]))) + ((secϕ_P_nm[i,j,2,2]*lnm7[2])*((S22M_t*cos_mλ[i,j,2])-(C22M_t*sin_mλ[i,j,2]))))/r_p4[i,j]
                        F_CS_ζ[i,j] = (((cosϕ_dP_nm[i,j,2,1])*((C21M_t*cos_mλ[i,j,1])+(S21M_t*sin_mλ[i,j,1]))) + ((cosϕ_dP_nm[i,j,2,2])*((C22M_t*cos_mλ[i,j,2])+(S22M_t*sin_mλ[i,j,2]))))/r_p4[i,j]

                        F_CS_ξ_36[i,j] = zero_q_1
                        F_CS_η_36[i,j] = zero_q_1
                        F_CS_ζ_36[i,j] = zero_q_1
                        for n in 3:n1SEM[mo]
                            for m in 1:n
                                temp_CS_ξ = (((P_nm[i,j,n,m]*lnm6[n])*((cos_mλ[i,j,m]*CM[n,m])+(sin_mλ[i,j,m]*SM[n,m])))/temp_rn[i,j,n]) + F_CS_ξ_36[i,j]
                                temp_CS_η = (((secϕ_P_nm[i,j,n,m]*lnm7[m])*((cos_mλ[i,j,m]*SM[n,m])-(sin_mλ[i,j,m]*CM[n,m])))/temp_rn[i,j,n]) + F_CS_η_36[i,j]
                                temp_CS_ζ = (((cosϕ_dP_nm[i,j,n,m])*((cos_mλ[i,j,m]*CM[n,m])+(sin_mλ[i,j,m]*SM[n,m])))/temp_rn[i,j,n]) + F_CS_ζ_36[i,j]
                                F_CS_ξ_36[i,j] = temp_CS_ξ
                                F_CS_η_36[i,j] = temp_CS_η
                                F_CS_ζ_36[i,j] = temp_CS_ζ
                            end
                        end
                        F_JCS_ξ[i,j] = (F_J_ξ[i,j] + F_J_ξ_36[i,j]) + (F_CS_ξ[i,j]+F_CS_ξ_36[i,j])
                        F_JCS_η[i,j] = (F_CS_η[i,j]+F_CS_η_36[i,j])
                        F_JCS_ζ[i,j] = (F_J_ζ[i,j] + F_J_ζ_36[i,j]) + (F_CS_ζ[i,j]+F_CS_ζ_36[i,j])
                    else
                        F_JCS_ξ[i,j] = (F_J_ξ[i,j] + F_J_ξ_36[i,j])
                        F_JCS_η[i,j] = zero_q_1
                        F_JCS_ζ[i,j] = (F_J_ζ[i,j] + F_J_ζ_36[i,j])
                    end

                    # R matrix: body-fixed -> "primed" ξηζ system
                    Rb2p[i,j,1,1] = cos_ϕ[i,j]*cos_λ[i,j]
                    Rb2p[i,j,2,1] = -sin_λ[i,j]
                    Rb2p[i,j,3,1] = -sin_ϕ[i,j]*cos_λ[i,j]
                    Rb2p[i,j,1,2] = cos_ϕ[i,j]*sin_λ[i,j]
                    Rb2p[i,j,2,2] = cos_λ[i,j]
                    Rb2p[i,j,3,2] = -sin_ϕ[i,j]*sin_λ[i,j]
                    Rb2p[i,j,1,3] = sin_ϕ[i,j]
                    Rb2p[i,j,2,3] = zero_q_1
                    Rb2p[i,j,3,3] = cos_ϕ[i,j]
                    # G matrix: space-fixed -> body-fixed -> "primed" ξηζ system
                    # G_{i,j} = \sum_k R_{i,k} M_{k,j}
                    Gc2p[i,j,1,1] = ((Rb2p[i,j,1,1]*M_[1,1,j]) + (Rb2p[i,j,1,2]*M_[2,1,j])) + (Rb2p[i,j,1,3]*M_[3,1,j])
                    Gc2p[i,j,2,1] = ((Rb2p[i,j,2,1]*M_[1,1,j]) + (Rb2p[i,j,2,2]*M_[2,1,j])) + (Rb2p[i,j,2,3]*M_[3,1,j])
                    Gc2p[i,j,3,1] = ((Rb2p[i,j,3,1]*M_[1,1,j]) + (Rb2p[i,j,3,2]*M_[2,1,j])) + (Rb2p[i,j,3,3]*M_[3,1,j])
                    Gc2p[i,j,1,2] = ((Rb2p[i,j,1,1]*M_[1,2,j]) + (Rb2p[i,j,1,2]*M_[2,2,j])) + (Rb2p[i,j,1,3]*M_[3,2,j])
                    Gc2p[i,j,2,2] = ((Rb2p[i,j,2,1]*M_[1,2,j]) + (Rb2p[i,j,2,2]*M_[2,2,j])) + (Rb2p[i,j,2,3]*M_[3,2,j])
                    Gc2p[i,j,3,2] = ((Rb2p[i,j,3,1]*M_[1,2,j]) + (Rb2p[i,j,3,2]*M_[2,2,j])) + (Rb2p[i,j,3,3]*M_[3,2,j])
                    Gc2p[i,j,1,3] = ((Rb2p[i,j,1,1]*M_[1,3,j]) + (Rb2p[i,j,1,2]*M_[2,3,j])) + (Rb2p[i,j,1,3]*M_[3,3,j])
                    Gc2p[i,j,2,3] = ((Rb2p[i,j,2,1]*M_[1,3,j]) + (Rb2p[i,j,2,2]*M_[2,3,j])) + (Rb2p[i,j,2,3]*M_[3,3,j])
                    Gc2p[i,j,3,3] = ((Rb2p[i,j,3,1]*M_[1,3,j]) + (Rb2p[i,j,3,2]*M_[2,3,j])) + (Rb2p[i,j,3,3]*M_[3,3,j])
                    # compute cartesian coordinates of acceleration due to body figure in inertial frame
                    F_JCS_x[i,j] = ((F_JCS_ξ[i,j]*Gc2p[i,j,1,1]) + (F_JCS_η[i,j]*Gc2p[i,j,2,1])) + (F_JCS_ζ[i,j]*Gc2p[i,j,3,1])
                    F_JCS_y[i,j] = ((F_JCS_ξ[i,j]*Gc2p[i,j,1,2]) + (F_JCS_η[i,j]*Gc2p[i,j,2,2])) + (F_JCS_ζ[i,j]*Gc2p[i,j,3,2])
                    F_JCS_z[i,j] = ((F_JCS_ξ[i,j]*Gc2p[i,j,1,3]) + (F_JCS_η[i,j]*Gc2p[i,j,2,3])) + (F_JCS_ζ[i,j]*Gc2p[i,j,3,3])
                end #if UJ_interaction[i,j]
            end # else (i != j)
        end #for i in 1:N
    end #for j in 1:N

    for j in 1:N
        for i in 1:N
            # i == j && continue
            if i == j
                continue
            else
                #Jn, Cnm, Snm accelerations, if j-th body is flattened
                if UJ_interaction[i,j]
                    # # add result to total acceleration upon j-th body figure due to i-th point mass
                    temp_accX_j[i,j] = accX[j] - (μ[i]*F_JCS_x[i,j])
                    accX[j] = temp_accX_j[i,j]
                    temp_accY_j[i,j] = accY[j] - (μ[i]*F_JCS_y[i,j])
                    accY[j] = temp_accY_j[i,j]
                    temp_accZ_j[i,j] = accZ[j] - (μ[i]*F_JCS_z[i,j])
                    accZ[j] = temp_accZ_j[i,j]

                    # # reaction force on i-th body
                    temp_accX_i[i,j] = accX[i] + (μ[j]*F_JCS_x[i,j])
                    accX[i] = temp_accX_i[i,j]
                    temp_accY_i[i,j] = accY[i] + (μ[j]*F_JCS_y[i,j])
                    accY[i] = temp_accY_i[i,j]
                    temp_accZ_i[i,j] = accZ[i] + (μ[j]*F_JCS_z[i,j])
                    accZ[i] = temp_accZ_i[i,j]
                end
            end # else (i != j)
        end
    end

    #post-Newtonian corrections to gravitational acceleration
    #Moyer, 1971, page 7 eq. 35
    # post-Newtonian iterative procedure setup and initialization
    for j in 1:N
        for i in 1:N
            # i == j && continue
            if i == j
                continue
            else
                _4ϕj[i,j] = 4newtonianNb_Potential[j]
                ϕi_plus_4ϕj[i,j] = newtonianNb_Potential[i] + _4ϕj[i,j]
                sj2_plus_2si2_minus_4vivj[i,j] = (v2[j] + (2v2[i])) - (4vi_dot_vj[i,j])
                ϕs_and_vs[i,j] = sj2_plus_2si2_minus_4vivj[i,j] - ϕi_plus_4ϕj[i,j]
                Xij_t_Ui = X[i,j]*dq[3i-2]
                Yij_t_Vi = Y[i,j]*dq[3i-1]
                Zij_t_Wi = Z[i,j]*dq[3i]
                Rij_dot_Vi = ( Xij_t_Ui+Yij_t_Vi ) + Zij_t_Wi
                # the expression below inside the (...)^2 should have a minus sign in front of the numerator,
                # but upon squaring it is eliminated, so at the end of the day, it is irrelevant ;)
                pn1t7 = (Rij_dot_Vi^2)/r_p2[i,j]
                pn1t2_7 = ϕs_and_vs[i,j] - (1.5pn1t7)
                pn1t1_7[i,j] = c_p2+pn1t2_7
                for k in 1:postnewton_iter
                    pn1[i,j,k] = zero_q_1
                    X_t_pn1[i,j,k] = zero_q_1
                    Y_t_pn1[i,j,k] = zero_q_1
                    Z_t_pn1[i,j,k] = zero_q_1
                    pNX_t_pn3[i,j,k] = zero_q_1
                    pNY_t_pn3[i,j,k] = zero_q_1
                    pNZ_t_pn3[i,j,k] = zero_q_1
                    pNX_t_X[i,j,k] = zero_q_1
                    pNY_t_Y[i,j,k] = zero_q_1
                    pNZ_t_Z[i,j,k] = zero_q_1
                end
            end # else (i != j)
        end
        postNewtonX[j,1] = newtonX[j]
        postNewtonY[j,1] = newtonY[j]
        postNewtonZ[j,1] = newtonZ[j]
        for k in 1:postnewton_iter
            pntempX[j,k] = zero_q_1
            pntempY[j,k] = zero_q_1
            pntempZ[j,k] = zero_q_1
        end
    end

    # post-Newtonian iterations
    for k in 1:postnewton_iter
        for j in 1:N
            for i in 1:N
                # i == j && continue
                if i == j
                    continue
                else
                    pNX_t_X[i,j,k] = postNewtonX[i,k]*X[i,j]
                    pNY_t_Y[i,j,k] = postNewtonY[i,k]*Y[i,j]
                    pNZ_t_Z[i,j,k] = postNewtonZ[i,k]*Z[i,j]
                    pn1[i,j,k] = (  pn1t1_7[i,j]  +  0.5*( (pNX_t_X[i,j,k]+pNY_t_Y[i,j,k]) + pNZ_t_Z[i,j,k] )  )

                    X_t_pn1[i,j,k] = newton_acc_X[i,j]*pn1[i,j,k]
                    Y_t_pn1[i,j,k] = newton_acc_Y[i,j]*pn1[i,j,k]
                    Z_t_pn1[i,j,k] = newton_acc_Z[i,j]*pn1[i,j,k]

                    pNX_t_pn3[i,j,k] = postNewtonX[i,k]*pn3[i,j]
                    pNY_t_pn3[i,j,k] = postNewtonY[i,k]*pn3[i,j]
                    pNZ_t_pn3[i,j,k] = postNewtonZ[i,k]*pn3[i,j]

                    termpnx = ( X_t_pn1[i,j,k] + (U_t_pn2[i,j]+pNX_t_pn3[i,j,k]) )
                    sumpnx = pntempX[j,k] + termpnx
                    pntempX[j,k] = sumpnx
                    termpny = ( Y_t_pn1[i,j,k] + (V_t_pn2[i,j]+pNY_t_pn3[i,j,k]) )
                    sumpny = pntempY[j,k] + termpny
                    pntempY[j,k] = sumpny
                    termpnz = ( Z_t_pn1[i,j,k] + (W_t_pn2[i,j]+pNZ_t_pn3[i,j,k]) )
                    sumpnz = pntempZ[j,k] + termpnz
                    pntempZ[j,k] = sumpnz
                end # else (i != j)
            end
            postNewtonX[j,k+1] = pntempX[j,k]*c_m2
            postNewtonY[j,k+1] = pntempY[j,k]*c_m2
            postNewtonZ[j,k+1] = pntempZ[j,k]*c_m2
        end
    end #for k in 1:postnewton_iter # (post-Newtonian iterations)

    #fill accelerations (post-Newtonian and extended body accelerations)
    for i in 1:N
        dq[3(N+i)-2] = postNewtonX[i,postnewton_iter+1] + accX[i]
        dq[3(N+i)-1] = postNewtonY[i,postnewton_iter+1] + accY[i]
        dq[3(N+i)  ] = postNewtonZ[i,postnewton_iter+1] + accZ[i]
    end

    nothing
end

function NBP_pN_A_J23E_J23M_J2S_threads!(dq, q, params, t)
    # N: number of bodies
    # S: auxiliary variable =eltype(q0)
    # eulang_de430_: Taylor interpolant for DE430 lunar orientation Euler angles
    local N, S, eulang_de430_, jd0 = params
    local N_ext = 11 # number of bodies in extended-body accelerations
    local eulang_t = eulang_de430_( (t+(jd0-J2000))*daysec )
    #local eulang_t_del = eulang_de430_( ((t-τ_M)+(jd0-J2000))*daysec )

    #TODO: handle appropiately @taylorize'd version with postnewton_iter>1
    local postnewton_iter = 1 # number of iterations of post-Newtonian subroutine

    # parameters related to speed of light, c
    local c_p2 = 29979.063823897606 # c^2 = 29979.063823897606 au^2/d^2
    local c_m2 = 3.3356611996764786e-5 # c^-2 = 3.3356611996764786e-5 d^2/au^2

    local zero_q_1 = zero(q[1])
    local one_t = one(t)

    X = Array{Taylor1{S}}(undef, N, N)
    Y = Array{Taylor1{S}}(undef, N, N)
    Z = Array{Taylor1{S}}(undef, N, N)

    r_p2 = Array{Taylor1{S}}(undef, N, N)
    r_p3d2 = Array{Taylor1{S}}(undef, N, N)
    r_p7d2 = Array{Taylor1{S}}(undef, N, N)

    #Newtonian accelerations
    newtonX = Array{Taylor1{S}}(undef, N)
    newtonY = Array{Taylor1{S}}(undef, N)
    newtonZ = Array{Taylor1{S}}(undef, N)

    newtonianCoeff = Array{Taylor1{S}}(undef, N, N)

    #post-Newtonian stuff
    U = Array{Taylor1{S}}(undef, N, N)
    V = Array{Taylor1{S}}(undef, N, N)
    W = Array{Taylor1{S}}(undef, N, N)

    _4U_m_3X = Array{Taylor1{S}}(undef, N, N)
    _4V_m_3Y = Array{Taylor1{S}}(undef, N, N)
    _4W_m_3Z = Array{Taylor1{S}}(undef, N, N)

    UU = Array{Taylor1{S}}(undef, N, N)
    VV = Array{Taylor1{S}}(undef, N, N)
    WW = Array{Taylor1{S}}(undef, N, N)

    r_p1d2 = Array{Taylor1{S}}(undef, N, N)

    newtonianNb_Potential = Array{Taylor1{S}}(undef, N)
    newtonian1b_Potential = Array{Taylor1{S}}(undef, N, N)
    newtonianCoeff = Array{Taylor1{S}}(undef, N, N)
    newton_acc_X = Array{Taylor1{S}}(undef, N, N)
    newton_acc_Y = Array{Taylor1{S}}(undef, N, N)
    newton_acc_Z = Array{Taylor1{S}}(undef, N, N)

    v2 = Array{Taylor1{S}}(undef, N)
    vi_dot_vj = Array{Taylor1{S}}(undef, N, N)
    pn2 = Array{Taylor1{S}}(undef, N, N)
    pn3 = Array{Taylor1{S}}(undef, N, N)
    _4ϕj = Array{Taylor1{S}}(undef, N, N)
    ϕi_plus_4ϕj = Array{Taylor1{S}}(undef, N, N)
    sj2_plus_2si2_minus_4vivj = Array{Taylor1{S}}(undef, N, N)
    ϕs_and_vs = Array{Taylor1{S}}(undef, N, N)
    U_t_pn2 = Array{Taylor1{S}}(undef, N, N)
    V_t_pn2 = Array{Taylor1{S}}(undef, N, N)
    W_t_pn2 = Array{Taylor1{S}}(undef, N, N)
    pn1t1_7 = Array{Taylor1{S}}(undef, N, N)

    pntempX = Array{Taylor1{S}}(undef, N, postnewton_iter)
    pntempY = Array{Taylor1{S}}(undef, N, postnewton_iter)
    pntempZ = Array{Taylor1{S}}(undef, N, postnewton_iter)
    pn1 = Array{Taylor1{S}}(undef, N, N, postnewton_iter)
    X_t_pn1 = Array{Taylor1{S}}(undef, N, N, postnewton_iter)
    Y_t_pn1 = Array{Taylor1{S}}(undef, N, N, postnewton_iter)
    Z_t_pn1 = Array{Taylor1{S}}(undef, N, N, postnewton_iter)
    pNX_t_pn3 = Array{Taylor1{S}}(undef, N, N, postnewton_iter)
    pNY_t_pn3 = Array{Taylor1{S}}(undef, N, N, postnewton_iter)
    pNZ_t_pn3 = Array{Taylor1{S}}(undef, N, N, postnewton_iter)
    pNX_t_X = Array{Taylor1{S}}(undef, N, N, postnewton_iter)
    pNY_t_Y = Array{Taylor1{S}}(undef, N, N, postnewton_iter)
    pNZ_t_Z = Array{Taylor1{S}}(undef, N, N, postnewton_iter)
    postNewtonX = Array{Taylor1{S}}(undef, N, postnewton_iter+1)
    postNewtonY = Array{Taylor1{S}}(undef, N, postnewton_iter+1)
    postNewtonZ = Array{Taylor1{S}}(undef, N, postnewton_iter+1)

    # (Jn, Cmn, Smn) acceleration auxiliaries
    X_bf_1 = Array{Taylor1{S}}(undef, N_ext, N_ext)
    Y_bf_1 = Array{Taylor1{S}}(undef, N_ext, N_ext)
    Z_bf_1 = Array{Taylor1{S}}(undef, N_ext, N_ext)
    X_bf_2 = Array{Taylor1{S}}(undef, N_ext, N_ext)
    Y_bf_2 = Array{Taylor1{S}}(undef, N_ext, N_ext)
    Z_bf_2 = Array{Taylor1{S}}(undef, N_ext, N_ext)
    X_bf_3 = Array{Taylor1{S}}(undef, N_ext, N_ext)
    Y_bf_3 = Array{Taylor1{S}}(undef, N_ext, N_ext)
    Z_bf_3 = Array{Taylor1{S}}(undef, N_ext, N_ext)
    X_bf = Array{Taylor1{S}}(undef, N_ext, N_ext)
    Y_bf = Array{Taylor1{S}}(undef, N_ext, N_ext)
    Z_bf = Array{Taylor1{S}}(undef, N_ext, N_ext)
    F_JCS_x = Array{Taylor1{S}}(undef, N_ext, N_ext)
    F_JCS_y = Array{Taylor1{S}}(undef, N_ext, N_ext)
    F_JCS_z = Array{Taylor1{S}}(undef, N_ext, N_ext)
    temp_accX_j = Array{Taylor1{S}}(undef, N_ext, N_ext)
    temp_accY_j = Array{Taylor1{S}}(undef, N_ext, N_ext)
    temp_accZ_j = Array{Taylor1{S}}(undef, N_ext, N_ext)
    temp_accX_i = Array{Taylor1{S}}(undef, N_ext, N_ext)
    temp_accY_i = Array{Taylor1{S}}(undef, N_ext, N_ext)
    temp_accZ_i = Array{Taylor1{S}}(undef, N_ext, N_ext)
    sin_ϕ = Array{Taylor1{S}}(undef, N_ext, N_ext)
    cos_ϕ = Array{Taylor1{S}}(undef, N_ext, N_ext)
    sin_λ = Array{Taylor1{S}}(undef, N_ext, N_ext)
    cos_λ = Array{Taylor1{S}}(undef, N_ext, N_ext)
    r_xy = Array{Taylor1{S}}(undef, N_ext, N_ext)
    r_p4 = Array{Taylor1{S}}(undef, N_ext, N_ext)
    P_n = Array{Taylor1{S}}(undef, N_ext, N_ext, maximum(n1SEM)+1)
    dP_n = Array{Taylor1{S}}(undef, N_ext, N_ext, maximum(n1SEM)+1)
    temp_fjξ = Array{Taylor1{S}}(undef, N_ext, N_ext, maximum(n1SEM)+1)
    temp_fjζ = Array{Taylor1{S}}(undef, N_ext, N_ext, maximum(n1SEM)+1)
    temp_rn = Array{Taylor1{S}}(undef, N_ext, N_ext, maximum(n1SEM)+1)
    F_CS_ξ_36 = Array{Taylor1{S}}(undef, N_ext, N_ext)
    F_CS_η_36 = Array{Taylor1{S}}(undef, N_ext, N_ext)
    F_CS_ζ_36 = Array{Taylor1{S}}(undef, N_ext, N_ext)
    F_J_ξ_36 = Array{Taylor1{S}}(undef, N_ext, N_ext)
    F_J_ζ_36 = Array{Taylor1{S}}(undef, N_ext, N_ext)
    sin_mλ = Array{Taylor1{S}}(undef, N_ext, N_ext, n1SEM[mo])
    cos_mλ = Array{Taylor1{S}}(undef, N_ext, N_ext, n1SEM[mo])
    secϕ_P_nm = Array{Taylor1{S}}(undef, N_ext, N_ext, n1SEM[mo]+1, n1SEM[mo]+1)
    P_nm = Array{Taylor1{S}}(undef, N_ext, N_ext, n1SEM[mo]+1, n1SEM[mo]+1)
    cosϕ_dP_nm = Array{Taylor1{S}}(undef, N_ext, N_ext, n1SEM[mo]+1, n1SEM[mo]+1)
    F_J_ξ = Array{Taylor1{S}}(undef, N_ext, N_ext)
    F_J_η = Array{Taylor1{S}}(undef, N_ext, N_ext)
    F_J_ζ = Array{Taylor1{S}}(undef, N_ext, N_ext)
    F_CS_ξ = Array{Taylor1{S}}(undef, N_ext, N_ext)
    F_CS_η = Array{Taylor1{S}}(undef, N_ext, N_ext)
    F_CS_ζ = Array{Taylor1{S}}(undef, N_ext, N_ext)
    F_JCS_ξ = Array{Taylor1{S}}(undef, N_ext, N_ext)
    F_JCS_η = Array{Taylor1{S}}(undef, N_ext, N_ext)
    F_JCS_ζ = Array{Taylor1{S}}(undef, N_ext, N_ext)
    Rb2p = Array{Taylor1{S}}(undef, N_ext, N_ext, 3, 3) #R matrix body-fixed to "primed" ξηζ frame (Moyer, 1971, eq. 161)
    Gc2p = Array{Taylor1{S}}(undef, N_ext, N_ext, 3, 3) #G matrix "space-fixed" to "primed" ξηζ frame (Moyer, 1971, eq. 163)

    # extended-body accelerations
    accX = Array{Taylor1{S}}(undef, N_ext)
    accY = Array{Taylor1{S}}(undef, N_ext)
    accZ = Array{Taylor1{S}}(undef, N_ext)

    # rotations to and from Earth, Sun and Moon pole-oriented frames
    local dsj2k = t+(jd0-2.451545e6) # J2000.0 = 2.451545e6
    local αs = deg2rad(α_p_sun*one_t)
    local δs = deg2rad(δ_p_sun*one_t)
    local αm = eulang_t[1] - (pi/2)
    local δm = (pi/2) - eulang_t[2]
    local Wm = eulang_t[3]
    local M_ = Array{Taylor1{S}}(undef, 3, 3, 5)
    local M_[:,:,ea] = c2t_jpl_de430(dsj2k)
    local M_[:,:,su] = pole_rotation(αs, δs)
    local M_[:,:,mo] = pole_rotation(αm, δm, Wm)
    local ITM_t = ITM_und.*one_t #+ ITM2(eulang_t_del[1], eulang_t_del[2], eulang_t_del[3])
    local fact_num = -4.5257273867882326e-36 # == -k_2M*μ[ea]*(R_moon^5)
    local fact1_jsem = [(2n-1)/n for n in 1:maximum(n1SEM)]
    local fact2_jsem = [(n-1)/n for n in 1:maximum(n1SEM)]
    local fact3_jsem = [n for n in 1:maximum(n1SEM)]
    local fact4_jsem = [n+1 for n in 1:maximum(n1SEM)]
    local fact5_jsem = [(n+2) for n in 1:maximum(n1SEM)]
    local lnm1 = [(2n-1)/(n-m) for n in 1:6, m in 1:6]
    local lnm2 = [-(n+m-1)/(n-m) for n in 1:6, m in 1:6]
    local lnm3 = [-n for n in 1:6]
    local lnm4 = [n+m for n in 1:6, m in 1:6]
    local lnm5 = [2n-1 for n in 1:6]
    local lnm6 = [-(n+1) for n in 1:6]
    local lnm7 = [m for m in 1:6]
    # TODO: solve differences between parsed and non-parsed
    local RE_au = (RE/au)
    local J2E_t = (J2E + J2EDOT*(dsj2k/yr))*(RE_au^2)
    local J2S_t = JSEM[su,2]*one_t
    local J2_t = Array{Taylor1{S}}(undef, 5)
    J2_t[su] = J2S_t
    J2_t[ea] = J2E_t

    Threads.@threads for j in 1:N
        newtonX[j] = zero_q_1
        newtonY[j] = zero_q_1
        newtonZ[j] = zero_q_1
        newtonianNb_Potential[j] = zero_q_1
        dq[3j-2] = q[3(N+j)-2]
        dq[3j-1] = q[3(N+j)-1]
        dq[3j  ] = q[3(N+j)  ]
    end

    Threads.@threads for j in 1:N_ext
        accX[j] = zero_q_1
        accY[j] = zero_q_1
        accZ[j] = zero_q_1
    end

    #compute point-mass Newtonian accelerations, all bodies
    Threads.@threads for j in 1:N
        for i in 1:N
            # i == j && continue
            if i == j
                continue
            else
                X[i,j] = q[3i-2]-q[3j-2]
                Y[i,j] = q[3i-1]-q[3j-1]
                Z[i,j] = q[3i]-q[3j]

                U[i,j] = dq[3i-2]-dq[3j-2]
                V[i,j] = dq[3i-1]-dq[3j-1]
                W[i,j] = dq[3i  ]-dq[3j  ]

                _4U_m_3X[i,j] = (4dq[3j-2])-(3dq[3i-2])
                _4V_m_3Y[i,j] = (4dq[3j-1])-(3dq[3i-1])
                _4W_m_3Z[i,j] = (4dq[3j  ])-(3dq[3i  ])

                pn2x = X[i,j]*_4U_m_3X[i,j]
                pn2y = Y[i,j]*_4V_m_3Y[i,j]
                pn2z = Z[i,j]*_4W_m_3Z[i,j]

                UU[i,j] = dq[3i-2]*dq[3j-2]
                VV[i,j] = dq[3i-1]*dq[3j-1]
                WW[i,j] = dq[3i  ]*dq[3j  ]

                vi_dot_vj[i,j] = ( UU[i,j]+VV[i,j] ) + WW[i,j]

                r_p2[i,j] = ( (X[i,j]^2)+(Y[i,j]^2) ) + (Z[i,j]^2)

                r_p1d2[i,j] = sqrt(r_p2[i,j])
                r_p3d2[i,j] = r_p2[i,j]^1.5
                r_p7d2[i,j] = r_p2[i,j]^3.5

                newtonianCoeff[i,j] =  μ[i]/r_p3d2[i,j]

                pn2[i,j] = newtonianCoeff[i,j]*(( pn2x+pn2y ) + pn2z)

                newton_acc_X[i,j] = X[i,j]*newtonianCoeff[i,j]
                newton_acc_Y[i,j] = Y[i,j]*newtonianCoeff[i,j]
                newton_acc_Z[i,j] = Z[i,j]*newtonianCoeff[i,j]

                newtonian1b_Potential[i,j] = μ[i]/r_p1d2[i, j]
                pn3[i,j] = 3.5newtonian1b_Potential[i,j]
                U_t_pn2[i,j] = pn2[i,j]*U[i,j]
                V_t_pn2[i,j] = pn2[i,j]*V[i,j]
                W_t_pn2[i,j] = pn2[i,j]*W[i,j]

                temp_001 = newtonX[j] + (X[i,j]*newtonianCoeff[i,j])
                newtonX[j] = temp_001
                temp_002 = newtonY[j] + (Y[i,j]*newtonianCoeff[i,j])
                newtonY[j] = temp_002
                temp_003 = newtonZ[j] + (Z[i,j]*newtonianCoeff[i,j])
                newtonZ[j] = temp_003
                temp_004 = newtonianNb_Potential[j] + newtonian1b_Potential[i, j]
                newtonianNb_Potential[j] = temp_004
            end # else (i != j)
        end #for, i
        v2[j] = ( (dq[3j-2]^2)+(dq[3j-1]^2) ) + (dq[3j]^2)
    end #for, j

    J2M_t = ( ITM_t[3,3] - ((ITM_t[1,1]+ITM_t[2,2])/2) )/(μ[mo]) # J_{2,M}*R_M^2
    C22M_t = ((ITM_t[2,2] - ITM_t[1,1])/(μ[mo]))/4 # C_{22,M}*R_M^2
    C21M_t = (-ITM_t[1,3])/(μ[mo]) # C_{21,M}*R_M^2
    S21M_t = (-ITM_t[3,2])/(μ[mo]) # S_{21,M}*R_M^2
    S22M_t = ((-ITM_t[2,1])/(μ[mo]))/2 # S_{22,M}*R_M^2
    J2_t[mo] = J2M_t

    Threads.@threads for j in 1:N_ext
        for i in 1:N_ext
            # i == j && continue
            if i == j
                continue
            else
                #Jn, Cnm, Snm accelerations, if j-th body is flattened
                if UJ_interaction[i,j]
                    # rotate from inertial frame to extended-body frame
                    X_bf_1[i,j] = X[i,j]*M_[1,1,j]
                    X_bf_2[i,j] = Y[i,j]*M_[1,2,j]
                    X_bf_3[i,j] = Z[i,j]*M_[1,3,j]
                    Y_bf_1[i,j] = X[i,j]*M_[2,1,j]
                    Y_bf_2[i,j] = Y[i,j]*M_[2,2,j]
                    Y_bf_3[i,j] = Z[i,j]*M_[2,3,j]
                    Z_bf_1[i,j] = X[i,j]*M_[3,1,j]
                    Z_bf_2[i,j] = Y[i,j]*M_[3,2,j]
                    Z_bf_3[i,j] = Z[i,j]*M_[3,3,j]
                    X_bf[i,j] = (X_bf_1[i,j] + X_bf_2[i,j]) + (X_bf_3[i,j]) # x-coordinate in body-fixed frame
                    Y_bf[i,j] = (Y_bf_1[i,j] + Y_bf_2[i,j]) + (Y_bf_3[i,j]) # y-coordinate in body-fixed frame
                    Z_bf[i,j] = (Z_bf_1[i,j] + Z_bf_2[i,j]) + (Z_bf_3[i,j]) # z-coordinate in body-fixed frame

                    # compute cartesian coordinates of acceleration due to body figure in body frame
                    sin_ϕ[i,j] = Z_bf[i,j]/r_p1d2[i,j] # Moyer (1971), eq. (165)
                    r_xy[i,j] = sqrt( (X_bf[i,j]^2)+(Y_bf[i,j]^2) )
                    cos_ϕ[i,j] = r_xy[i,j]/r_p1d2[i,j] # Moyer (1971), eq. (166)
                    sin_λ[i,j] = Y_bf[i,j]/r_xy[i,j] # Moyer (1971), eq. (167)
                    cos_λ[i,j] = X_bf[i,j]/r_xy[i,j] # Moyer (1971), eq. (168)

                    # compute accelerations due to zonal harmonics J_{n}
                    P_n[i,j,1] = one_t
                    P_n[i,j,2] = sin_ϕ[i,j]
                    dP_n[i,j,1] = zero_q_1
                    dP_n[i,j,2] = one_t
                    for n in 2:n1SEM[j] #min(3,n1SEM[j])
                        P_n[i,j,n+1] = ((P_n[i,j,n]*sin_ϕ[i,j])*fact1_jsem[n]) - (P_n[i,j,n-1]*fact2_jsem[n])
                        dP_n[i,j,n+1] = (dP_n[i,j,n]*sin_ϕ[i,j]) + (P_n[i,j,n]*fact3_jsem[n])
                        temp_rn[i,j,n] = r_p1d2[i,j]^fact5_jsem[n]
                    end
                    r_p4[i,j] = r_p2[i,j]^2
                    # if j == mo
                    #     F_J_ξ[i,j] = ((P_n[i,j,3]*fact4_jsem[2])*J2M_t)/r_p4[i,j]
                    #     F_J_ζ[i,j] = (((-dP_n[i,j,3])*cos_ϕ[i,j])*J2M_t)/r_p4[i,j]
                    # elseif j == ea
                    #     F_J_ξ[i,j] = ((P_n[i,j,3]*fact4_jsem[2])*J2E_t)/r_p4[i,j]
                    #     F_J_ζ[i,j] = (((-dP_n[i,j,3])*cos_ϕ[i,j])*J2E_t)/r_p4[i,j]
                    # elseif j == su
                    #     F_J_ξ[i,j] = ((P_n[i,j,3]*fact4_jsem[2])*J2S_t)/r_p4[i,j]
                    #     F_J_ζ[i,j] = (((-dP_n[i,j,3])*cos_ϕ[i,j])*J2S_t)/r_p4[i,j]
                    # end
                    F_J_ξ[i,j] = ((P_n[i,j,3]*fact4_jsem[2])*J2_t[j])/r_p4[i,j]
                    F_J_ζ[i,j] = (((-dP_n[i,j,3])*cos_ϕ[i,j])*J2_t[j])/r_p4[i,j]
                    F_J_ξ_36[i,j] = zero_q_1
                    F_J_ζ_36[i,j] = zero_q_1
                    for n in 3:n1SEM[j] #min(3,n1SEM[j])
                        temp_fjξ[i,j,n] = (((P_n[i,j,n+1]*fact4_jsem[n])*JSEM[j,n])/temp_rn[i,j,n]) + F_J_ξ_36[i,j]
                        temp_fjζ[i,j,n] = ((((-dP_n[i,j,n+1])*cos_ϕ[i,j])*JSEM[j,n])/temp_rn[i,j,n]) + F_J_ζ_36[i,j]
                        F_J_ξ_36[i,j] = temp_fjξ[i,j,n]
                        F_J_ζ_36[i,j] = temp_fjζ[i,j,n]
                    end

                    # Moon: compute accelerations due to tesseral harmonics C_{nm}, S_{nm}
                    if j == mo
                        # Associate Legendre polynomials recursion
                        for m in 1:n1SEM[mo]
                            if m == 1
                                sin_mλ[i,j,1] = sin_λ[i,j] # Moyer (1971), eq. (167)
                                cos_mλ[i,j,1] = cos_λ[i,j] # Moyer (1971), eq. (168)
                                secϕ_P_nm[i,j,1,1] = one_t
                            else
                                sin_mλ[i,j,m] = (sin_mλ[i,j,1]*cos_mλ[i,j,m-1]) + (cos_mλ[i,j,1]*sin_mλ[i,j,m-1])
                                cos_mλ[i,j,m] = (cos_mλ[i,j,1]*cos_mλ[i,j,m-1]) - (sin_mλ[i,j,1]*sin_mλ[i,j,m-1])
                                secϕ_P_nm[i,j,m,m] = (secϕ_P_nm[i,j,m-1,m-1]*cos_ϕ[i,j])*lnm5[m]
                                P_nm[i,j,m,m] = secϕ_P_nm[i,j,m,m]*cos_ϕ[i,j]
                                cosϕ_dP_nm[i,j,m,m] = (secϕ_P_nm[i,j,m,m]*sin_ϕ[i,j])*lnm3[m] #+0 (second term in Eq. 183 from Moyer, 1971, vanishes when n=m)
                            end
                            for n in m+1:n1SEM[mo]
                                if n == m+1
                                    secϕ_P_nm[i,j,n,m] = (secϕ_P_nm[i,j,n-1,m]*sin_ϕ[i,j])*lnm1[n,m]
                                else
                                    secϕ_P_nm[i,j,n,m] = ((secϕ_P_nm[i,j,n-1,m]*sin_ϕ[i,j])*lnm1[n,m]) + (secϕ_P_nm[i,j,n-2,m]*lnm2[n,m])
                                end
                                P_nm[i,j,n,m] = secϕ_P_nm[i,j,n,m]*cos_ϕ[i,j]
                                cosϕ_dP_nm[i,j,n,m] = ((secϕ_P_nm[i,j,n,m]*sin_ϕ[i,j])*lnm3[n]) + (secϕ_P_nm[i,j,n-1,m]*lnm4[n,m])
                            end
                        end

                        F_CS_ξ[i,j] = (((P_nm[i,j,2,1]*lnm6[2])*((C21M_t*cos_mλ[i,j,1])+(S21M_t*sin_mλ[i,j,1]))) + ((P_nm[i,j,2,2]*lnm6[2])*((C22M_t*cos_mλ[i,j,2])+(S22M_t*sin_mλ[i,j,2]))))/r_p4[i,j]
                        F_CS_η[i,j] = (((secϕ_P_nm[i,j,2,1]*lnm7[1])*((S21M_t*cos_mλ[i,j,1])-(C21M_t*sin_mλ[i,j,1]))) + ((secϕ_P_nm[i,j,2,2]*lnm7[2])*((S22M_t*cos_mλ[i,j,2])-(C22M_t*sin_mλ[i,j,2]))))/r_p4[i,j]
                        F_CS_ζ[i,j] = (((cosϕ_dP_nm[i,j,2,1])*((C21M_t*cos_mλ[i,j,1])+(S21M_t*sin_mλ[i,j,1]))) + ((cosϕ_dP_nm[i,j,2,2])*((C22M_t*cos_mλ[i,j,2])+(S22M_t*sin_mλ[i,j,2]))))/r_p4[i,j]

                        F_CS_ξ_36[i,j] = zero_q_1
                        F_CS_η_36[i,j] = zero_q_1
                        F_CS_ζ_36[i,j] = zero_q_1
                        for n in 3:n1SEM[mo]
                            for m in 1:n
                                temp_CS_ξ = (((P_nm[i,j,n,m]*lnm6[n])*((cos_mλ[i,j,m]*CM[n,m])+(sin_mλ[i,j,m]*SM[n,m])))/temp_rn[i,j,n]) + F_CS_ξ_36[i,j]
                                temp_CS_η = (((secϕ_P_nm[i,j,n,m]*lnm7[m])*((cos_mλ[i,j,m]*SM[n,m])-(sin_mλ[i,j,m]*CM[n,m])))/temp_rn[i,j,n]) + F_CS_η_36[i,j]
                                temp_CS_ζ = (((cosϕ_dP_nm[i,j,n,m])*((cos_mλ[i,j,m]*CM[n,m])+(sin_mλ[i,j,m]*SM[n,m])))/temp_rn[i,j,n]) + F_CS_ζ_36[i,j]
                                F_CS_ξ_36[i,j] = temp_CS_ξ
                                F_CS_η_36[i,j] = temp_CS_η
                                F_CS_ζ_36[i,j] = temp_CS_ζ
                            end
                        end
                        F_JCS_ξ[i,j] = (F_J_ξ[i,j] + F_J_ξ_36[i,j]) + (F_CS_ξ[i,j]+F_CS_ξ_36[i,j])
                        F_JCS_η[i,j] = (F_CS_η[i,j]+F_CS_η_36[i,j])
                        F_JCS_ζ[i,j] = (F_J_ζ[i,j] + F_J_ζ_36[i,j]) + (F_CS_ζ[i,j]+F_CS_ζ_36[i,j])
                    else
                        F_JCS_ξ[i,j] = (F_J_ξ[i,j] + F_J_ξ_36[i,j])
                        F_JCS_η[i,j] = zero_q_1
                        F_JCS_ζ[i,j] = (F_J_ζ[i,j] + F_J_ζ_36[i,j])
                    end

                    # R matrix: body-fixed -> "primed" ξηζ system
                    Rb2p[i,j,1,1] = cos_ϕ[i,j]*cos_λ[i,j]
                    Rb2p[i,j,2,1] = -sin_λ[i,j]
                    Rb2p[i,j,3,1] = -sin_ϕ[i,j]*cos_λ[i,j]
                    Rb2p[i,j,1,2] = cos_ϕ[i,j]*sin_λ[i,j]
                    Rb2p[i,j,2,2] = cos_λ[i,j]
                    Rb2p[i,j,3,2] = -sin_ϕ[i,j]*sin_λ[i,j]
                    Rb2p[i,j,1,3] = sin_ϕ[i,j]
                    Rb2p[i,j,2,3] = zero_q_1
                    Rb2p[i,j,3,3] = cos_ϕ[i,j]
                    # G matrix: space-fixed -> body-fixed -> "primed" ξηζ system
                    # G_{i,j} = \sum_k R_{i,k} M_{k,j}
                    Gc2p[i,j,1,1] = ((Rb2p[i,j,1,1]*M_[1,1,j]) + (Rb2p[i,j,1,2]*M_[2,1,j])) + (Rb2p[i,j,1,3]*M_[3,1,j])
                    Gc2p[i,j,2,1] = ((Rb2p[i,j,2,1]*M_[1,1,j]) + (Rb2p[i,j,2,2]*M_[2,1,j])) + (Rb2p[i,j,2,3]*M_[3,1,j])
                    Gc2p[i,j,3,1] = ((Rb2p[i,j,3,1]*M_[1,1,j]) + (Rb2p[i,j,3,2]*M_[2,1,j])) + (Rb2p[i,j,3,3]*M_[3,1,j])
                    Gc2p[i,j,1,2] = ((Rb2p[i,j,1,1]*M_[1,2,j]) + (Rb2p[i,j,1,2]*M_[2,2,j])) + (Rb2p[i,j,1,3]*M_[3,2,j])
                    Gc2p[i,j,2,2] = ((Rb2p[i,j,2,1]*M_[1,2,j]) + (Rb2p[i,j,2,2]*M_[2,2,j])) + (Rb2p[i,j,2,3]*M_[3,2,j])
                    Gc2p[i,j,3,2] = ((Rb2p[i,j,3,1]*M_[1,2,j]) + (Rb2p[i,j,3,2]*M_[2,2,j])) + (Rb2p[i,j,3,3]*M_[3,2,j])
                    Gc2p[i,j,1,3] = ((Rb2p[i,j,1,1]*M_[1,3,j]) + (Rb2p[i,j,1,2]*M_[2,3,j])) + (Rb2p[i,j,1,3]*M_[3,3,j])
                    Gc2p[i,j,2,3] = ((Rb2p[i,j,2,1]*M_[1,3,j]) + (Rb2p[i,j,2,2]*M_[2,3,j])) + (Rb2p[i,j,2,3]*M_[3,3,j])
                    Gc2p[i,j,3,3] = ((Rb2p[i,j,3,1]*M_[1,3,j]) + (Rb2p[i,j,3,2]*M_[2,3,j])) + (Rb2p[i,j,3,3]*M_[3,3,j])
                    # compute cartesian coordinates of acceleration due to body figure in inertial frame
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
                #Jn, Cnm, Snm accelerations, if j-th body is flattened
                if UJ_interaction[i,j]
                    # # add result to total acceleration upon j-th body figure due to i-th point mass
                    temp_accX_j[i,j] = accX[j] - (μ[i]*F_JCS_x[i,j])
                    accX[j] = temp_accX_j[i,j]
                    temp_accY_j[i,j] = accY[j] - (μ[i]*F_JCS_y[i,j])
                    accY[j] = temp_accY_j[i,j]
                    temp_accZ_j[i,j] = accZ[j] - (μ[i]*F_JCS_z[i,j])
                    accZ[j] = temp_accZ_j[i,j]

                    # # reaction force on i-th body
                    temp_accX_i[i,j] = accX[i] + (μ[j]*F_JCS_x[i,j])
                    accX[i] = temp_accX_i[i,j]
                    temp_accY_i[i,j] = accY[i] + (μ[j]*F_JCS_y[i,j])
                    accY[i] = temp_accY_i[i,j]
                    temp_accZ_i[i,j] = accZ[i] + (μ[j]*F_JCS_z[i,j])
                    accZ[i] = temp_accZ_i[i,j]
                end
            end # else (i != j)
        end
    end

    #post-Newtonian corrections to gravitational acceleration
    #Moyer, 1971, page 7 eq. 35
    # post-Newtonian iterative procedure setup and initialization
    Threads.@threads for j in 1:N
        for i in 1:N
            # i == j && continue
            if i == j
                continue
            else
                _4ϕj[i,j] = 4newtonianNb_Potential[j]
                ϕi_plus_4ϕj[i,j] = newtonianNb_Potential[i] + _4ϕj[i,j]
                sj2_plus_2si2_minus_4vivj[i,j] = (v2[j] + (2v2[i])) - (4vi_dot_vj[i,j])
                ϕs_and_vs[i,j] = sj2_plus_2si2_minus_4vivj[i,j] - ϕi_plus_4ϕj[i,j]
                Xij_t_Ui = X[i,j]*dq[3i-2]
                Yij_t_Vi = Y[i,j]*dq[3i-1]
                Zij_t_Wi = Z[i,j]*dq[3i]
                Rij_dot_Vi = ( Xij_t_Ui+Yij_t_Vi ) + Zij_t_Wi
                # the expression below inside the (...)^2 should have a minus sign in front of the numerator,
                # but upon squaring it is eliminated, so at the end of the day, it is irrelevant ;)
                pn1t7 = (Rij_dot_Vi^2)/r_p2[i,j]
                pn1t2_7 = ϕs_and_vs[i,j] - (1.5pn1t7)
                pn1t1_7[i,j] = c_p2+pn1t2_7
                for k in 1:postnewton_iter
                    pn1[i,j,k] = zero_q_1
                    X_t_pn1[i,j,k] = zero_q_1
                    Y_t_pn1[i,j,k] = zero_q_1
                    Z_t_pn1[i,j,k] = zero_q_1
                    pNX_t_pn3[i,j,k] = zero_q_1
                    pNY_t_pn3[i,j,k] = zero_q_1
                    pNZ_t_pn3[i,j,k] = zero_q_1
                    pNX_t_X[i,j,k] = zero_q_1
                    pNY_t_Y[i,j,k] = zero_q_1
                    pNZ_t_Z[i,j,k] = zero_q_1
                end
            end # else (i != j)
        end
        postNewtonX[j,1] = newtonX[j]
        postNewtonY[j,1] = newtonY[j]
        postNewtonZ[j,1] = newtonZ[j]
        for k in 1:postnewton_iter
            pntempX[j,k] = zero_q_1
            pntempY[j,k] = zero_q_1
            pntempZ[j,k] = zero_q_1
        end
    end

    # post-Newtonian iterations
    for k in 1:postnewton_iter
        Threads.@threads for j in 1:N
            for i in 1:N
                # i == j && continue
                if i == j
                    continue
                else
                    pNX_t_X[i,j,k] = postNewtonX[i,k]*X[i,j]
                    pNY_t_Y[i,j,k] = postNewtonY[i,k]*Y[i,j]
                    pNZ_t_Z[i,j,k] = postNewtonZ[i,k]*Z[i,j]
                    pn1[i,j,k] = (  pn1t1_7[i,j]  +  0.5*( (pNX_t_X[i,j,k]+pNY_t_Y[i,j,k]) + pNZ_t_Z[i,j,k] )  )

                    X_t_pn1[i,j,k] = newton_acc_X[i,j]*pn1[i,j,k]
                    Y_t_pn1[i,j,k] = newton_acc_Y[i,j]*pn1[i,j,k]
                    Z_t_pn1[i,j,k] = newton_acc_Z[i,j]*pn1[i,j,k]

                    pNX_t_pn3[i,j,k] = postNewtonX[i,k]*pn3[i,j]
                    pNY_t_pn3[i,j,k] = postNewtonY[i,k]*pn3[i,j]
                    pNZ_t_pn3[i,j,k] = postNewtonZ[i,k]*pn3[i,j]

                    termpnx = ( X_t_pn1[i,j,k] + (U_t_pn2[i,j]+pNX_t_pn3[i,j,k]) )
                    sumpnx = pntempX[j,k] + termpnx
                    pntempX[j,k] = sumpnx
                    termpny = ( Y_t_pn1[i,j,k] + (V_t_pn2[i,j]+pNY_t_pn3[i,j,k]) )
                    sumpny = pntempY[j,k] + termpny
                    pntempY[j,k] = sumpny
                    termpnz = ( Z_t_pn1[i,j,k] + (W_t_pn2[i,j]+pNZ_t_pn3[i,j,k]) )
                    sumpnz = pntempZ[j,k] + termpnz
                    pntempZ[j,k] = sumpnz
                end # else (i != j)
            end
            postNewtonX[j,k+1] = pntempX[j,k]*c_m2
            postNewtonY[j,k+1] = pntempY[j,k]*c_m2
            postNewtonZ[j,k+1] = pntempZ[j,k]*c_m2
        end
    end #for k in 1:postnewton_iter # (post-Newtonian iterations)

    #fill accelerations (post-Newtonian and extended body accelerations)
    Threads.@threads for i in 1:N_ext
        dq[3(N+i)-2] = postNewtonX[i,postnewton_iter+1] + accX[i]
        dq[3(N+i)-1] = postNewtonY[i,postnewton_iter+1] + accY[i]
        dq[3(N+i)  ] = postNewtonZ[i,postnewton_iter+1] + accZ[i]
    end
    Threads.@threads for i in N_ext+1:N
        dq[3(N+i)-2] = postNewtonX[i,postnewton_iter+1]
        dq[3(N+i)-1] = postNewtonY[i,postnewton_iter+1]
        dq[3(N+i)  ] = postNewtonZ[i,postnewton_iter+1]
    end

    nothing
end

function DE430!(dq, q, params, t)
    # N: number of bodies
    # S: auxiliary variable =eltype(q0)
    # eulang_de430_: Taylor interpolant for DE430 lunar orientation Euler angles
    local N, S, eulang_de430_, jd0 = params
    local N_ext = 11 # number of bodies in extended-body accelerations
    local N_back = 11 # number of bodies in backward integration
    local params_back = (N_back, S, eulang_de430_, jd0)
    local qq_ = Taylor1.(constant_term.(q[union(1:3N_back, 3N+1:3N+3N_back)]), t.order)
    local dqq_ = similar(qq_)
    # local xaux_ = similar(qq_)
    # local jtcffs = TaylorIntegration.__jetcoeffs!(Val(false), NBP_pN_A_J23E_J23M_J2S_threads!, t, qq_, dqq_, xaux_, params_back)
    # local jtcffs = TaylorIntegration.__jetcoeffs!(Val(true), NBP_pN_A_J23E_J23M_J2S_threads!, t, qq_, dqq_, xaux_, params_back)
    # local jtcffs = TaylorIntegration.jetcoeffs!(NBP_pN_A_J23E_J23M_J2S_threads!, t, qq_, dqq_, xaux_, params_back)
    local jtcffs = TaylorIntegration.jetcoeffs!(Val(NBP_pN_A_J23E_J23M_J2S_threads!), t, qq_, dqq_, params_back)
    local __t = Taylor1(t.order)
    local q_del_τ_M = qq_(__t-τ_M)
    local q_del_τ_0 = qq_(__t-τ_0p)
    local q_del_τ_1 = qq_(__t-τ_1p)
    local q_del_τ_2 = qq_(__t-τ_2p)

    local dsj2k = t+(jd0-2.451545e6) # days since J2000.0 (TDB)
    local eulang_t = eulang_de430_( dsj2k*daysec )
    local eulang_t_del = eulang_de430_( (dsj2k-τ_M)*daysec )

    #TODO: handle appropiately @taylorize'd version with postnewton_iter>1
    local postnewton_iter = 1 # number of iterations of post-Newtonian subroutine

    # parameters related to speed of light, c
    local c_p2 = 29979.063823897606 # c^2 = 29979.063823897606 au^2/d^2
    local c_m2 = 3.3356611996764786e-5 # c^-2 = 3.3356611996764786e-5 d^2/au^2

    local zero_q_1 = zero(q[1])
    local one_t = one(t)

    X = Array{Taylor1{S}}(undef, N, N)
    Y = Array{Taylor1{S}}(undef, N, N)
    Z = Array{Taylor1{S}}(undef, N, N)

    r_p2 = Array{Taylor1{S}}(undef, N, N)
    r_p3d2 = Array{Taylor1{S}}(undef, N, N)
    r_p7d2 = Array{Taylor1{S}}(undef, N, N)

    #Newtonian accelerations
    newtonX = Array{Taylor1{S}}(undef, N)
    newtonY = Array{Taylor1{S}}(undef, N)
    newtonZ = Array{Taylor1{S}}(undef, N)

    newtonianCoeff = Array{Taylor1{S}}(undef, N, N)

    #post-Newtonian stuff
    U = Array{Taylor1{S}}(undef, N, N)
    V = Array{Taylor1{S}}(undef, N, N)
    W = Array{Taylor1{S}}(undef, N, N)

    _4U_m_3X = Array{Taylor1{S}}(undef, N, N)
    _4V_m_3Y = Array{Taylor1{S}}(undef, N, N)
    _4W_m_3Z = Array{Taylor1{S}}(undef, N, N)

    UU = Array{Taylor1{S}}(undef, N, N)
    VV = Array{Taylor1{S}}(undef, N, N)
    WW = Array{Taylor1{S}}(undef, N, N)

    r_p1d2 = Array{Taylor1{S}}(undef, N, N)

    newtonianNb_Potential = Array{Taylor1{S}}(undef, N)
    newtonian1b_Potential = Array{Taylor1{S}}(undef, N, N)
    newtonianCoeff = Array{Taylor1{S}}(undef, N, N)
    newton_acc_X = Array{Taylor1{S}}(undef, N, N)
    newton_acc_Y = Array{Taylor1{S}}(undef, N, N)
    newton_acc_Z = Array{Taylor1{S}}(undef, N, N)

    v2 = Array{Taylor1{S}}(undef, N)
    vi_dot_vj = Array{Taylor1{S}}(undef, N, N)
    pn2 = Array{Taylor1{S}}(undef, N, N)
    pn3 = Array{Taylor1{S}}(undef, N, N)
    _4ϕj = Array{Taylor1{S}}(undef, N, N)
    ϕi_plus_4ϕj = Array{Taylor1{S}}(undef, N, N)
    sj2_plus_2si2_minus_4vivj = Array{Taylor1{S}}(undef, N, N)
    ϕs_and_vs = Array{Taylor1{S}}(undef, N, N)
    U_t_pn2 = Array{Taylor1{S}}(undef, N, N)
    V_t_pn2 = Array{Taylor1{S}}(undef, N, N)
    W_t_pn2 = Array{Taylor1{S}}(undef, N, N)
    pn1t1_7 = Array{Taylor1{S}}(undef, N, N)

    pntempX = Array{Taylor1{S}}(undef, N, postnewton_iter)
    pntempY = Array{Taylor1{S}}(undef, N, postnewton_iter)
    pntempZ = Array{Taylor1{S}}(undef, N, postnewton_iter)
    pn1 = Array{Taylor1{S}}(undef, N, N, postnewton_iter)
    X_t_pn1 = Array{Taylor1{S}}(undef, N, N, postnewton_iter)
    Y_t_pn1 = Array{Taylor1{S}}(undef, N, N, postnewton_iter)
    Z_t_pn1 = Array{Taylor1{S}}(undef, N, N, postnewton_iter)
    pNX_t_pn3 = Array{Taylor1{S}}(undef, N, N, postnewton_iter)
    pNY_t_pn3 = Array{Taylor1{S}}(undef, N, N, postnewton_iter)
    pNZ_t_pn3 = Array{Taylor1{S}}(undef, N, N, postnewton_iter)
    pNX_t_X = Array{Taylor1{S}}(undef, N, N, postnewton_iter)
    pNY_t_Y = Array{Taylor1{S}}(undef, N, N, postnewton_iter)
    pNZ_t_Z = Array{Taylor1{S}}(undef, N, N, postnewton_iter)
    postNewtonX = Array{Taylor1{S}}(undef, N, postnewton_iter+1)
    postNewtonY = Array{Taylor1{S}}(undef, N, postnewton_iter+1)
    postNewtonZ = Array{Taylor1{S}}(undef, N, postnewton_iter+1)

    # (Jn, Cmn, Smn) acceleration auxiliaries
    X_bf_1 = Array{Taylor1{S}}(undef, N_ext, N_ext)
    Y_bf_1 = Array{Taylor1{S}}(undef, N_ext, N_ext)
    Z_bf_1 = Array{Taylor1{S}}(undef, N_ext, N_ext)
    X_bf_2 = Array{Taylor1{S}}(undef, N_ext, N_ext)
    Y_bf_2 = Array{Taylor1{S}}(undef, N_ext, N_ext)
    Z_bf_2 = Array{Taylor1{S}}(undef, N_ext, N_ext)
    X_bf_3 = Array{Taylor1{S}}(undef, N_ext, N_ext)
    Y_bf_3 = Array{Taylor1{S}}(undef, N_ext, N_ext)
    Z_bf_3 = Array{Taylor1{S}}(undef, N_ext, N_ext)
    X_bf = Array{Taylor1{S}}(undef, N_ext, N_ext)
    Y_bf = Array{Taylor1{S}}(undef, N_ext, N_ext)
    Z_bf = Array{Taylor1{S}}(undef, N_ext, N_ext)
    F_JCS_x = Array{Taylor1{S}}(undef, N_ext, N_ext)
    F_JCS_y = Array{Taylor1{S}}(undef, N_ext, N_ext)
    F_JCS_z = Array{Taylor1{S}}(undef, N_ext, N_ext)
    temp_accX_j = Array{Taylor1{S}}(undef, N_ext, N_ext)
    temp_accY_j = Array{Taylor1{S}}(undef, N_ext, N_ext)
    temp_accZ_j = Array{Taylor1{S}}(undef, N_ext, N_ext)
    temp_accX_i = Array{Taylor1{S}}(undef, N_ext, N_ext)
    temp_accY_i = Array{Taylor1{S}}(undef, N_ext, N_ext)
    temp_accZ_i = Array{Taylor1{S}}(undef, N_ext, N_ext)
    sin_ϕ = Array{Taylor1{S}}(undef, N_ext, N_ext)
    cos_ϕ = Array{Taylor1{S}}(undef, N_ext, N_ext)
    sin_λ = Array{Taylor1{S}}(undef, N_ext, N_ext)
    cos_λ = Array{Taylor1{S}}(undef, N_ext, N_ext)
    r_xy = Array{Taylor1{S}}(undef, N_ext, N_ext)
    r_p4 = Array{Taylor1{S}}(undef, N_ext, N_ext)
    P_n = Array{Taylor1{S}}(undef, N_ext, N_ext, maximum(n1SEM)+1)
    dP_n = Array{Taylor1{S}}(undef, N_ext, N_ext, maximum(n1SEM)+1)
    temp_fjξ = Array{Taylor1{S}}(undef, N_ext, N_ext, maximum(n1SEM)+1)
    temp_fjζ = Array{Taylor1{S}}(undef, N_ext, N_ext, maximum(n1SEM)+1)
    temp_rn = Array{Taylor1{S}}(undef, N_ext, N_ext, maximum(n1SEM)+1)
    F_CS_ξ_36 = Array{Taylor1{S}}(undef, N_ext, N_ext)
    F_CS_η_36 = Array{Taylor1{S}}(undef, N_ext, N_ext)
    F_CS_ζ_36 = Array{Taylor1{S}}(undef, N_ext, N_ext)
    F_J_ξ_36 = Array{Taylor1{S}}(undef, N_ext, N_ext)
    F_J_ζ_36 = Array{Taylor1{S}}(undef, N_ext, N_ext)
    sin_mλ = Array{Taylor1{S}}(undef, N_ext, N_ext, n1SEM[mo])
    cos_mλ = Array{Taylor1{S}}(undef, N_ext, N_ext, n1SEM[mo])
    secϕ_P_nm = Array{Taylor1{S}}(undef, N_ext, N_ext, n1SEM[mo]+1, n1SEM[mo]+1)
    P_nm = Array{Taylor1{S}}(undef, N_ext, N_ext, n1SEM[mo]+1, n1SEM[mo]+1)
    cosϕ_dP_nm = Array{Taylor1{S}}(undef, N_ext, N_ext, n1SEM[mo]+1, n1SEM[mo]+1)
    F_J_ξ = Array{Taylor1{S}}(undef, N_ext, N_ext)
    F_J_η = Array{Taylor1{S}}(undef, N_ext, N_ext)
    F_J_ζ = Array{Taylor1{S}}(undef, N_ext, N_ext)
    F_CS_ξ = Array{Taylor1{S}}(undef, N_ext, N_ext)
    F_CS_η = Array{Taylor1{S}}(undef, N_ext, N_ext)
    F_CS_ζ = Array{Taylor1{S}}(undef, N_ext, N_ext)
    F_JCS_ξ = Array{Taylor1{S}}(undef, N_ext, N_ext)
    F_JCS_η = Array{Taylor1{S}}(undef, N_ext, N_ext)
    F_JCS_ζ = Array{Taylor1{S}}(undef, N_ext, N_ext)
    Rb2p = Array{Taylor1{S}}(undef, N_ext, N_ext, 3, 3) #R matrix body-fixed to "primed" ξηζ frame (Moyer, 1971, eq. 161)
    Gc2p = Array{Taylor1{S}}(undef, N_ext, N_ext, 3, 3) #G matrix "space-fixed" to "primed" ξηζ frame (Moyer, 1971, eq. 163)

    # extended-body accelerations
    accX = Array{Taylor1{S}}(undef, N_ext)
    accY = Array{Taylor1{S}}(undef, N_ext)
    accZ = Array{Taylor1{S}}(undef, N_ext)

    # tidal accelerations
    r_star_M_0 = Array{Taylor1{S}}(undef, 3)
    r_star_S_0 = Array{Taylor1{S}}(undef, 3)
    r_star_M_1 = Array{Taylor1{S}}(undef, 3)
    r_star_S_1 = Array{Taylor1{S}}(undef, 3)
    r_star_M_2 = Array{Taylor1{S}}(undef, 3)
    r_star_S_2 = Array{Taylor1{S}}(undef, 3)

    # rotations to and from Earth, Sun and Moon pole-oriented frames
    local αs = deg2rad(α_p_sun*one_t)
    local δs = deg2rad(δ_p_sun*one_t)
    local αm = eulang_t[1] - (pi/2)
    local δm = (pi/2) - eulang_t[2]
    local Wm = eulang_t[3]
    local M_ = Array{Taylor1{S}}(undef, 3, 3, 5)
    local M_[:,:,ea] = c2t_jpl_de430(dsj2k)
    local M_[:,:,su] = pole_rotation(αs, δs)
    local M_[:,:,mo] = pole_rotation(αm, δm, Wm)
    local M_del_mo = pole_rotation(eulang_t_del[1] - (pi/2), (pi/2) - eulang_t_del[2], eulang_t_del[3])
    ITM_t = Array{Taylor1{S}}(undef, 3, 3)
    ITM2_t = Array{Taylor1{S}}(undef, 3, 3)
    local ITM2_t = ITM_und.*one_t + ITM2(eulang_t_del[1], eulang_t_del[2], eulang_t_del[3])
    local fact_num = -4.5257273867882326e-36 # == -k_2M*μ[ea]*(R_moon^5)
    local fact1_jsem = [(2n-1)/n for n in 1:maximum(n1SEM)]
    local fact2_jsem = [(n-1)/n for n in 1:maximum(n1SEM)]
    local fact3_jsem = [n for n in 1:maximum(n1SEM)]
    local fact4_jsem = [n+1 for n in 1:maximum(n1SEM)]
    local fact5_jsem = [(n+2) for n in 1:maximum(n1SEM)]
    local lnm1 = [(2n-1)/(n-m) for n in 1:6, m in 1:6]
    local lnm2 = [-(n+m-1)/(n-m) for n in 1:6, m in 1:6]
    local lnm3 = [-n for n in 1:6]
    local lnm4 = [n+m for n in 1:6, m in 1:6]
    local lnm5 = [2n-1 for n in 1:6]
    local lnm6 = [-(n+1) for n in 1:6]
    local lnm7 = [m for m in 1:6]
    # TODO: solve differences between parsed and non-parsed
    local RE_au = (RE/au)
    local J2E_t = (J2E + J2EDOT*(dsj2k/yr))*(RE_au^2)
    local J2S_t = JSEM[su,2]*one_t
    local J2_t = Array{Taylor1{S}}(undef, 5)
    J2_t[su] = J2S_t
    J2_t[ea] = J2E_t
    # Moon tidal acc: geocentric space-fixed -> rotational time-delay -> geocentric Earth true-equator-of-date frame
    local R30 = c2t_jpl_de430(dsj2k-τ_0p) #M_[:,:,ea] # #Rz(-ω_E*τ_0) == Id(3x3), since τ_0=0
    local R31 = Rz(-ω_E*τ_1)*c2t_jpl_de430(dsj2k-τ_1p) # *R30
    local R32 = Rz(-ω_E*τ_2)*c2t_jpl_de430(dsj2k-τ_2p) # *R30
    local tid_num_coeff = 1.5*(1.0 + μ[mo]/μ[ea])

    Threads.@threads for j in 1:N
        newtonX[j] = zero_q_1
        newtonY[j] = zero_q_1
        newtonZ[j] = zero_q_1
        newtonianNb_Potential[j] = zero_q_1
        dq[3j-2] = q[3(N+j)-2]
        dq[3j-1] = q[3(N+j)-1]
        dq[3j  ] = q[3(N+j)  ]
    end

    Threads.@threads for j in 1:N_ext
        accX[j] = zero_q_1
        accY[j] = zero_q_1
        accZ[j] = zero_q_1
    end

    #compute point-mass Newtonian accelerations, all bodies
    Threads.@threads for j in 1:N
        for i in 1:N
            # i == j && continue
            if i == j
                continue
            else
                X[i,j] = q[3i-2]-q[3j-2]
                Y[i,j] = q[3i-1]-q[3j-1]
                Z[i,j] = q[3i]-q[3j]

                U[i,j] = dq[3i-2]-dq[3j-2]
                V[i,j] = dq[3i-1]-dq[3j-1]
                W[i,j] = dq[3i  ]-dq[3j  ]

                _4U_m_3X[i,j] = (4dq[3j-2])-(3dq[3i-2])
                _4V_m_3Y[i,j] = (4dq[3j-1])-(3dq[3i-1])
                _4W_m_3Z[i,j] = (4dq[3j  ])-(3dq[3i  ])

                pn2x = X[i,j]*_4U_m_3X[i,j]
                pn2y = Y[i,j]*_4V_m_3Y[i,j]
                pn2z = Z[i,j]*_4W_m_3Z[i,j]

                UU[i,j] = dq[3i-2]*dq[3j-2]
                VV[i,j] = dq[3i-1]*dq[3j-1]
                WW[i,j] = dq[3i  ]*dq[3j  ]

                vi_dot_vj[i,j] = ( UU[i,j]+VV[i,j] ) + WW[i,j]

                r_p2[i,j] = ( (X[i,j]^2)+(Y[i,j]^2) ) + (Z[i,j]^2)

                r_p1d2[i,j] = sqrt(r_p2[i,j])
                r_p3d2[i,j] = r_p2[i,j]^1.5
                r_p7d2[i,j] = r_p2[i,j]^3.5

                newtonianCoeff[i,j] =  μ[i]/r_p3d2[i,j]

                pn2[i,j] = newtonianCoeff[i,j]*(( pn2x+pn2y ) + pn2z)

                newton_acc_X[i,j] = X[i,j]*newtonianCoeff[i,j]
                newton_acc_Y[i,j] = Y[i,j]*newtonianCoeff[i,j]
                newton_acc_Z[i,j] = Z[i,j]*newtonianCoeff[i,j]

                newtonian1b_Potential[i,j] = μ[i]/r_p1d2[i, j]
                pn3[i,j] = 3.5newtonian1b_Potential[i,j]
                U_t_pn2[i,j] = pn2[i,j]*U[i,j]
                V_t_pn2[i,j] = pn2[i,j]*V[i,j]
                W_t_pn2[i,j] = pn2[i,j]*W[i,j]

                temp_001 = newtonX[j] + (X[i,j]*newtonianCoeff[i,j])
                newtonX[j] = temp_001
                temp_002 = newtonY[j] + (Y[i,j]*newtonianCoeff[i,j])
                newtonY[j] = temp_002
                temp_003 = newtonZ[j] + (Z[i,j]*newtonianCoeff[i,j])
                newtonZ[j] = temp_003
                temp_004 = newtonianNb_Potential[j] + newtonian1b_Potential[i, j]
                newtonianNb_Potential[j] = temp_004
            end # else (i != j)
        end #for, i
        v2[j] = ( (dq[3j-2]^2)+(dq[3j-1]^2) ) + (dq[3j]^2)
    end #for, j

    # transform delayed geocentric position of Moon (space-fixed->lunar mantle frame)
    X_me_del_τ_M = q_del_τ_M[3mo-2] - q_del_τ_M[3ea-2]
    Y_me_del_τ_M = q_del_τ_M[3mo-1] - q_del_τ_M[3ea-1]
    Z_me_del_τ_M = q_del_τ_M[3mo  ] - q_del_τ_M[3ea  ]
    xmed = ((M_del_mo[1,1]*X_me_del_τ_M)+(M_del_mo[1,2]*Y_me_del_τ_M)) + (M_del_mo[1,3]*Z_me_del_τ_M)
    ymed = ((M_del_mo[2,1]*X_me_del_τ_M)+(M_del_mo[2,2]*Y_me_del_τ_M)) + (M_del_mo[2,3]*Z_me_del_τ_M)
    zmed = ((M_del_mo[3,1]*X_me_del_τ_M)+(M_del_mo[3,2]*Y_me_del_τ_M)) + (M_del_mo[3,3]*Z_me_del_τ_M)
    # compute matrix elements of lunar moment of inertia (Folkner et al. 2014, eq. 41)
    rmed2 = ((xmed^2)+(ymed^2))+(zmed^2)
    factmed = fact_num/(rmed2^2.5)
    ITM_t[1,1] = ITM2_t[1,1] + ( factmed*((xmed^2)-(rmed2/3)) )
    ITM_t[2,2] = ITM2_t[2,2] + ( factmed*((ymed^2)-(rmed2/3)) )
    ITM_t[3,3] = ITM2_t[3,3] + ( factmed*((zmed^2)-(rmed2/3)) )
    ITM_t[1,2] = ITM2_t[1,2] + (factmed*(xmed*ymed))
    ITM_t[2,1] = ITM_t[1,2]
    ITM_t[1,3] = ITM2_t[1,3] + (factmed*(xmed*zmed))
    ITM_t[3,1] = ITM_t[1,3]
    ITM_t[2,3] = ITM2_t[2,3] + (factmed*(ymed*zmed))
    ITM_t[3,2] = ITM_t[2,3]
    J2M_t = ( ITM_t[3,3] - ((ITM_t[1,1]+ITM_t[2,2])/2) )/(μ[mo]) # J_{2,M}*R_M^2
    C22M_t = ((ITM_t[2,2] - ITM_t[1,1])/(μ[mo]))/4 # C_{22,M}*R_M^2
    C21M_t = (-ITM_t[1,3])/(μ[mo]) # C_{21,M}*R_M^2
    S21M_t = (-ITM_t[3,2])/(μ[mo]) # S_{21,M}*R_M^2
    S22M_t = ((-ITM_t[2,1])/(μ[mo]))/2 # S_{22,M}*R_M^2
    J2_t[mo] = J2M_t

    Threads.@threads for j in 1:N_ext
        for i in 1:N_ext
            # i == j && continue
            if i == j
                continue
            else
                #Jn, Cnm, Snm accelerations, if j-th body is flattened
                if UJ_interaction[i,j]
                    # rotate from inertial frame to extended-body frame
                    X_bf_1[i,j] = X[i,j]*M_[1,1,j]
                    X_bf_2[i,j] = Y[i,j]*M_[1,2,j]
                    X_bf_3[i,j] = Z[i,j]*M_[1,3,j]
                    Y_bf_1[i,j] = X[i,j]*M_[2,1,j]
                    Y_bf_2[i,j] = Y[i,j]*M_[2,2,j]
                    Y_bf_3[i,j] = Z[i,j]*M_[2,3,j]
                    Z_bf_1[i,j] = X[i,j]*M_[3,1,j]
                    Z_bf_2[i,j] = Y[i,j]*M_[3,2,j]
                    Z_bf_3[i,j] = Z[i,j]*M_[3,3,j]
                    X_bf[i,j] = (X_bf_1[i,j] + X_bf_2[i,j]) + (X_bf_3[i,j]) # x-coordinate in body-fixed frame
                    Y_bf[i,j] = (Y_bf_1[i,j] + Y_bf_2[i,j]) + (Y_bf_3[i,j]) # y-coordinate in body-fixed frame
                    Z_bf[i,j] = (Z_bf_1[i,j] + Z_bf_2[i,j]) + (Z_bf_3[i,j]) # z-coordinate in body-fixed frame

                    # compute cartesian coordinates of acceleration due to body figure in body frame
                    sin_ϕ[i,j] = Z_bf[i,j]/r_p1d2[i,j] # Moyer (1971), eq. (165)
                    r_xy[i,j] = sqrt( (X_bf[i,j]^2)+(Y_bf[i,j]^2) )
                    cos_ϕ[i,j] = r_xy[i,j]/r_p1d2[i,j] # Moyer (1971), eq. (166)
                    sin_λ[i,j] = Y_bf[i,j]/r_xy[i,j] # Moyer (1971), eq. (167)
                    cos_λ[i,j] = X_bf[i,j]/r_xy[i,j] # Moyer (1971), eq. (168)

                    # compute accelerations due to zonal harmonics J_{n}
                    P_n[i,j,1] = one_t
                    P_n[i,j,2] = sin_ϕ[i,j]
                    dP_n[i,j,1] = zero_q_1
                    dP_n[i,j,2] = one_t
                    for n in 2:n1SEM[j] #min(3,n1SEM[j])
                        P_n[i,j,n+1] = ((P_n[i,j,n]*sin_ϕ[i,j])*fact1_jsem[n]) - (P_n[i,j,n-1]*fact2_jsem[n])
                        dP_n[i,j,n+1] = (dP_n[i,j,n]*sin_ϕ[i,j]) + (P_n[i,j,n]*fact3_jsem[n])
                        temp_rn[i,j,n] = r_p1d2[i,j]^fact5_jsem[n]
                    end
                    r_p4[i,j] = r_p2[i,j]^2
                    # if j == mo
                    #     F_J_ξ[i,j] = ((P_n[i,j,3]*fact4_jsem[2])*J2M_t)/r_p4[i,j]
                    #     F_J_ζ[i,j] = (((-dP_n[i,j,3])*cos_ϕ[i,j])*J2M_t)/r_p4[i,j]
                    # elseif j == ea
                    #     F_J_ξ[i,j] = ((P_n[i,j,3]*fact4_jsem[2])*J2E_t)/r_p4[i,j]
                    #     F_J_ζ[i,j] = (((-dP_n[i,j,3])*cos_ϕ[i,j])*J2E_t)/r_p4[i,j]
                    # elseif j == su
                    #     F_J_ξ[i,j] = ((P_n[i,j,3]*fact4_jsem[2])*J2S_t)/r_p4[i,j]
                    #     F_J_ζ[i,j] = (((-dP_n[i,j,3])*cos_ϕ[i,j])*J2S_t)/r_p4[i,j]
                    # end
                    F_J_ξ[i,j] = ((P_n[i,j,3]*fact4_jsem[2])*J2_t[j])/r_p4[i,j]
                    F_J_ζ[i,j] = (((-dP_n[i,j,3])*cos_ϕ[i,j])*J2_t[j])/r_p4[i,j]
                    F_J_ξ_36[i,j] = zero_q_1
                    F_J_ζ_36[i,j] = zero_q_1
                    for n in 3:n1SEM[j] #min(3,n1SEM[j])
                        temp_fjξ[i,j,n] = (((P_n[i,j,n+1]*fact4_jsem[n])*JSEM[j,n])/temp_rn[i,j,n]) + F_J_ξ_36[i,j]
                        temp_fjζ[i,j,n] = ((((-dP_n[i,j,n+1])*cos_ϕ[i,j])*JSEM[j,n])/temp_rn[i,j,n]) + F_J_ζ_36[i,j]
                        F_J_ξ_36[i,j] = temp_fjξ[i,j,n]
                        F_J_ζ_36[i,j] = temp_fjζ[i,j,n]
                    end

                    # Moon: compute accelerations due to tesseral harmonics C_{nm}, S_{nm}
                    if j == mo
                        # Associate Legendre polynomials recursion
                        for m in 1:n1SEM[mo]
                            if m == 1
                                sin_mλ[i,j,1] = sin_λ[i,j] # Moyer (1971), eq. (167)
                                cos_mλ[i,j,1] = cos_λ[i,j] # Moyer (1971), eq. (168)
                                secϕ_P_nm[i,j,1,1] = one_t
                            else
                                sin_mλ[i,j,m] = (sin_mλ[i,j,1]*cos_mλ[i,j,m-1]) + (cos_mλ[i,j,1]*sin_mλ[i,j,m-1])
                                cos_mλ[i,j,m] = (cos_mλ[i,j,1]*cos_mλ[i,j,m-1]) - (sin_mλ[i,j,1]*sin_mλ[i,j,m-1])
                                secϕ_P_nm[i,j,m,m] = (secϕ_P_nm[i,j,m-1,m-1]*cos_ϕ[i,j])*lnm5[m]
                                P_nm[i,j,m,m] = secϕ_P_nm[i,j,m,m]*cos_ϕ[i,j]
                                cosϕ_dP_nm[i,j,m,m] = (secϕ_P_nm[i,j,m,m]*sin_ϕ[i,j])*lnm3[m] #+0 (second term in Eq. 183 from Moyer, 1971, vanishes when n=m)
                            end
                            for n in m+1:n1SEM[mo]
                                if n == m+1
                                    secϕ_P_nm[i,j,n,m] = (secϕ_P_nm[i,j,n-1,m]*sin_ϕ[i,j])*lnm1[n,m]
                                else
                                    secϕ_P_nm[i,j,n,m] = ((secϕ_P_nm[i,j,n-1,m]*sin_ϕ[i,j])*lnm1[n,m]) + (secϕ_P_nm[i,j,n-2,m]*lnm2[n,m])
                                end
                                P_nm[i,j,n,m] = secϕ_P_nm[i,j,n,m]*cos_ϕ[i,j]
                                cosϕ_dP_nm[i,j,n,m] = ((secϕ_P_nm[i,j,n,m]*sin_ϕ[i,j])*lnm3[n]) + (secϕ_P_nm[i,j,n-1,m]*lnm4[n,m])
                            end
                        end

                        F_CS_ξ[i,j] = (((P_nm[i,j,2,1]*lnm6[2])*((C21M_t*cos_mλ[i,j,1])+(S21M_t*sin_mλ[i,j,1]))) + ((P_nm[i,j,2,2]*lnm6[2])*((C22M_t*cos_mλ[i,j,2])+(S22M_t*sin_mλ[i,j,2]))))/r_p4[i,j]
                        F_CS_η[i,j] = (((secϕ_P_nm[i,j,2,1]*lnm7[1])*((S21M_t*cos_mλ[i,j,1])-(C21M_t*sin_mλ[i,j,1]))) + ((secϕ_P_nm[i,j,2,2]*lnm7[2])*((S22M_t*cos_mλ[i,j,2])-(C22M_t*sin_mλ[i,j,2]))))/r_p4[i,j]
                        F_CS_ζ[i,j] = (((cosϕ_dP_nm[i,j,2,1])*((C21M_t*cos_mλ[i,j,1])+(S21M_t*sin_mλ[i,j,1]))) + ((cosϕ_dP_nm[i,j,2,2])*((C22M_t*cos_mλ[i,j,2])+(S22M_t*sin_mλ[i,j,2]))))/r_p4[i,j]

                        F_CS_ξ_36[i,j] = zero_q_1
                        F_CS_η_36[i,j] = zero_q_1
                        F_CS_ζ_36[i,j] = zero_q_1
                        for n in 3:n1SEM[mo]
                            for m in 1:n
                                temp_CS_ξ = (((P_nm[i,j,n,m]*lnm6[n])*((cos_mλ[i,j,m]*CM[n,m])+(sin_mλ[i,j,m]*SM[n,m])))/temp_rn[i,j,n]) + F_CS_ξ_36[i,j]
                                temp_CS_η = (((secϕ_P_nm[i,j,n,m]*lnm7[m])*((cos_mλ[i,j,m]*SM[n,m])-(sin_mλ[i,j,m]*CM[n,m])))/temp_rn[i,j,n]) + F_CS_η_36[i,j]
                                temp_CS_ζ = (((cosϕ_dP_nm[i,j,n,m])*((cos_mλ[i,j,m]*CM[n,m])+(sin_mλ[i,j,m]*SM[n,m])))/temp_rn[i,j,n]) + F_CS_ζ_36[i,j]
                                F_CS_ξ_36[i,j] = temp_CS_ξ
                                F_CS_η_36[i,j] = temp_CS_η
                                F_CS_ζ_36[i,j] = temp_CS_ζ
                            end
                        end
                        F_JCS_ξ[i,j] = (F_J_ξ[i,j] + F_J_ξ_36[i,j]) + (F_CS_ξ[i,j]+F_CS_ξ_36[i,j])
                        F_JCS_η[i,j] = (F_CS_η[i,j]+F_CS_η_36[i,j])
                        F_JCS_ζ[i,j] = (F_J_ζ[i,j] + F_J_ζ_36[i,j]) + (F_CS_ζ[i,j]+F_CS_ζ_36[i,j])
                    else
                        F_JCS_ξ[i,j] = (F_J_ξ[i,j] + F_J_ξ_36[i,j])
                        F_JCS_η[i,j] = zero_q_1
                        F_JCS_ζ[i,j] = (F_J_ζ[i,j] + F_J_ζ_36[i,j])
                    end

                    # R matrix: body-fixed -> "primed" ξηζ system
                    Rb2p[i,j,1,1] = cos_ϕ[i,j]*cos_λ[i,j]
                    Rb2p[i,j,2,1] = -sin_λ[i,j]
                    Rb2p[i,j,3,1] = -sin_ϕ[i,j]*cos_λ[i,j]
                    Rb2p[i,j,1,2] = cos_ϕ[i,j]*sin_λ[i,j]
                    Rb2p[i,j,2,2] = cos_λ[i,j]
                    Rb2p[i,j,3,2] = -sin_ϕ[i,j]*sin_λ[i,j]
                    Rb2p[i,j,1,3] = sin_ϕ[i,j]
                    Rb2p[i,j,2,3] = zero_q_1
                    Rb2p[i,j,3,3] = cos_ϕ[i,j]
                    # G matrix: space-fixed -> body-fixed -> "primed" ξηζ system
                    # G_{i,j} = \sum_k R_{i,k} M_{k,j}
                    Gc2p[i,j,1,1] = ((Rb2p[i,j,1,1]*M_[1,1,j]) + (Rb2p[i,j,1,2]*M_[2,1,j])) + (Rb2p[i,j,1,3]*M_[3,1,j])
                    Gc2p[i,j,2,1] = ((Rb2p[i,j,2,1]*M_[1,1,j]) + (Rb2p[i,j,2,2]*M_[2,1,j])) + (Rb2p[i,j,2,3]*M_[3,1,j])
                    Gc2p[i,j,3,1] = ((Rb2p[i,j,3,1]*M_[1,1,j]) + (Rb2p[i,j,3,2]*M_[2,1,j])) + (Rb2p[i,j,3,3]*M_[3,1,j])
                    Gc2p[i,j,1,2] = ((Rb2p[i,j,1,1]*M_[1,2,j]) + (Rb2p[i,j,1,2]*M_[2,2,j])) + (Rb2p[i,j,1,3]*M_[3,2,j])
                    Gc2p[i,j,2,2] = ((Rb2p[i,j,2,1]*M_[1,2,j]) + (Rb2p[i,j,2,2]*M_[2,2,j])) + (Rb2p[i,j,2,3]*M_[3,2,j])
                    Gc2p[i,j,3,2] = ((Rb2p[i,j,3,1]*M_[1,2,j]) + (Rb2p[i,j,3,2]*M_[2,2,j])) + (Rb2p[i,j,3,3]*M_[3,2,j])
                    Gc2p[i,j,1,3] = ((Rb2p[i,j,1,1]*M_[1,3,j]) + (Rb2p[i,j,1,2]*M_[2,3,j])) + (Rb2p[i,j,1,3]*M_[3,3,j])
                    Gc2p[i,j,2,3] = ((Rb2p[i,j,2,1]*M_[1,3,j]) + (Rb2p[i,j,2,2]*M_[2,3,j])) + (Rb2p[i,j,2,3]*M_[3,3,j])
                    Gc2p[i,j,3,3] = ((Rb2p[i,j,3,1]*M_[1,3,j]) + (Rb2p[i,j,3,2]*M_[2,3,j])) + (Rb2p[i,j,3,3]*M_[3,3,j])
                    # compute cartesian coordinates of acceleration due to body figure in inertial frame
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
                #Jn, Cnm, Snm accelerations, if j-th body is flattened
                if UJ_interaction[i,j]
                    # # add result to total acceleration upon j-th body figure due to i-th point mass
                    temp_accX_j[i,j] = accX[j] - (μ[i]*F_JCS_x[i,j])
                    accX[j] = temp_accX_j[i,j]
                    temp_accY_j[i,j] = accY[j] - (μ[i]*F_JCS_y[i,j])
                    accY[j] = temp_accY_j[i,j]
                    temp_accZ_j[i,j] = accZ[j] - (μ[i]*F_JCS_z[i,j])
                    accZ[j] = temp_accZ_j[i,j]

                    # # reaction force on i-th body
                    temp_accX_i[i,j] = accX[i] + (μ[j]*F_JCS_x[i,j])
                    accX[i] = temp_accX_i[i,j]
                    temp_accY_i[i,j] = accY[i] + (μ[j]*F_JCS_y[i,j])
                    accY[i] = temp_accY_i[i,j]
                    temp_accZ_i[i,j] = accZ[i] + (μ[j]*F_JCS_z[i,j])
                    accZ[i] = temp_accZ_i[i,j]
                end
            end # else (i != j)
        end
    end

    #post-Newtonian corrections to gravitational acceleration
    #Moyer, 1971, page 7 eq. 35
    # post-Newtonian iterative procedure setup and initialization
    Threads.@threads for j in 1:N
        for i in 1:N
            # i == j && continue
            if i == j
                continue
            else
                _4ϕj[i,j] = 4newtonianNb_Potential[j]
                ϕi_plus_4ϕj[i,j] = newtonianNb_Potential[i] + _4ϕj[i,j]
                sj2_plus_2si2_minus_4vivj[i,j] = (v2[j] + (2v2[i])) - (4vi_dot_vj[i,j])
                ϕs_and_vs[i,j] = sj2_plus_2si2_minus_4vivj[i,j] - ϕi_plus_4ϕj[i,j]
                Xij_t_Ui = X[i,j]*dq[3i-2]
                Yij_t_Vi = Y[i,j]*dq[3i-1]
                Zij_t_Wi = Z[i,j]*dq[3i]
                Rij_dot_Vi = ( Xij_t_Ui+Yij_t_Vi ) + Zij_t_Wi
                # the expression below inside the (...)^2 should have a minus sign in front of the numerator,
                # but upon squaring it is eliminated, so at the end of the day, it is irrelevant ;)
                pn1t7 = (Rij_dot_Vi^2)/r_p2[i,j]
                pn1t2_7 = ϕs_and_vs[i,j] - (1.5pn1t7)
                pn1t1_7[i,j] = c_p2+pn1t2_7
                for k in 1:postnewton_iter
                    pn1[i,j,k] = zero_q_1
                    X_t_pn1[i,j,k] = zero_q_1
                    Y_t_pn1[i,j,k] = zero_q_1
                    Z_t_pn1[i,j,k] = zero_q_1
                    pNX_t_pn3[i,j,k] = zero_q_1
                    pNY_t_pn3[i,j,k] = zero_q_1
                    pNZ_t_pn3[i,j,k] = zero_q_1
                    pNX_t_X[i,j,k] = zero_q_1
                    pNY_t_Y[i,j,k] = zero_q_1
                    pNZ_t_Z[i,j,k] = zero_q_1
                end
            end # else (i != j)
        end
        postNewtonX[j,1] = newtonX[j]
        postNewtonY[j,1] = newtonY[j]
        postNewtonZ[j,1] = newtonZ[j]
        for k in 1:postnewton_iter
            pntempX[j,k] = zero_q_1
            pntempY[j,k] = zero_q_1
            pntempZ[j,k] = zero_q_1
        end
    end

    # post-Newtonian iterations
    for k in 1:postnewton_iter
        Threads.@threads for j in 1:N
            for i in 1:N
                # i == j && continue
                if i == j
                    continue
                else
                    pNX_t_X[i,j,k] = postNewtonX[i,k]*X[i,j]
                    pNY_t_Y[i,j,k] = postNewtonY[i,k]*Y[i,j]
                    pNZ_t_Z[i,j,k] = postNewtonZ[i,k]*Z[i,j]
                    pn1[i,j,k] = (  pn1t1_7[i,j]  +  0.5*( (pNX_t_X[i,j,k]+pNY_t_Y[i,j,k]) + pNZ_t_Z[i,j,k] )  )

                    X_t_pn1[i,j,k] = newton_acc_X[i,j]*pn1[i,j,k]
                    Y_t_pn1[i,j,k] = newton_acc_Y[i,j]*pn1[i,j,k]
                    Z_t_pn1[i,j,k] = newton_acc_Z[i,j]*pn1[i,j,k]

                    pNX_t_pn3[i,j,k] = postNewtonX[i,k]*pn3[i,j]
                    pNY_t_pn3[i,j,k] = postNewtonY[i,k]*pn3[i,j]
                    pNZ_t_pn3[i,j,k] = postNewtonZ[i,k]*pn3[i,j]

                    termpnx = ( X_t_pn1[i,j,k] + (U_t_pn2[i,j]+pNX_t_pn3[i,j,k]) )
                    sumpnx = pntempX[j,k] + termpnx
                    pntempX[j,k] = sumpnx
                    termpny = ( Y_t_pn1[i,j,k] + (V_t_pn2[i,j]+pNY_t_pn3[i,j,k]) )
                    sumpny = pntempY[j,k] + termpny
                    pntempY[j,k] = sumpny
                    termpnz = ( Z_t_pn1[i,j,k] + (W_t_pn2[i,j]+pNZ_t_pn3[i,j,k]) )
                    sumpnz = pntempZ[j,k] + termpnz
                    pntempZ[j,k] = sumpnz
                end # else (i != j)
            end
            postNewtonX[j,k+1] = pntempX[j,k]*c_m2
            postNewtonY[j,k+1] = pntempY[j,k]*c_m2
            postNewtonZ[j,k+1] = pntempZ[j,k]*c_m2
        end
    end #for k in 1:postnewton_iter # (post-Newtonian iterations)

    # compute acceleration of the Moon due to tides raised on Earth by the Sun and Moon
    # time-delayed barycentric Earth position
    X_E_τ_0 = q_del_τ_0[3ea-2]
    Y_E_τ_0 = q_del_τ_0[3ea-1]
    Z_E_τ_0 = q_del_τ_0[3ea  ]
    X_E_τ_1 = q_del_τ_1[3ea-2]
    Y_E_τ_1 = q_del_τ_1[3ea-1]
    Z_E_τ_1 = q_del_τ_1[3ea  ]
    X_E_τ_2 = q_del_τ_2[3ea-2]
    Y_E_τ_2 = q_del_τ_2[3ea-1]
    Z_E_τ_2 = q_del_τ_2[3ea  ]
    # time-delayed geocentric Moon position
    X_ME_τ_0 = q_del_τ_0[3mo-2] - X_E_τ_0
    Y_ME_τ_0 = q_del_τ_0[3mo-1] - Y_E_τ_0
    Z_ME_τ_0 = q_del_τ_0[3mo  ] - Z_E_τ_0
    X_ME_τ_1 = q_del_τ_1[3mo-2] - X_E_τ_1
    Y_ME_τ_1 = q_del_τ_1[3mo-1] - Y_E_τ_1
    Z_ME_τ_1 = q_del_τ_1[3mo  ] - Z_E_τ_1
    X_ME_τ_2 = q_del_τ_2[3mo-2] - X_E_τ_2
    Y_ME_τ_2 = q_del_τ_2[3mo-1] - Y_E_τ_2
    Z_ME_τ_2 = q_del_τ_2[3mo  ] - Z_E_τ_2
    # time-delayed geocentric Sun position
    X_SE_τ_0 = q_del_τ_0[3su-2] - X_E_τ_0
    Y_SE_τ_0 = q_del_τ_0[3su-1] - Y_E_τ_0
    Z_SE_τ_0 = q_del_τ_0[3su  ] - Z_E_τ_0
    X_SE_τ_1 = q_del_τ_1[3su-2] - X_E_τ_1
    Y_SE_τ_1 = q_del_τ_1[3su-1] - Y_E_τ_1
    Z_SE_τ_1 = q_del_τ_1[3su  ] - Z_E_τ_1
    X_SE_τ_2 = q_del_τ_2[3su-2] - X_E_τ_2
    Y_SE_τ_2 = q_del_τ_2[3su-1] - Y_E_τ_2
    Z_SE_τ_2 = q_del_τ_2[3su  ] - Z_E_τ_2

    # Folkner et al. (2014), Eq. (31)
    # note: "starred" vectors should be initialized to zero
    # r-star 0, Moon
    r_star_M_0[1] = ((R30[1,1]*X_ME_τ_0) + (R30[1,2]*Y_ME_τ_0)) + (R30[1,3]*Z_ME_τ_0)
    r_star_M_0[2] = ((R30[2,1]*X_ME_τ_0) + (R30[2,2]*Y_ME_τ_0)) + (R30[2,3]*Z_ME_τ_0)
    r_star_M_0[3] = ((R30[3,1]*X_ME_τ_0) + (R30[3,2]*Y_ME_τ_0)) + (R30[3,3]*Z_ME_τ_0)
    # r-star 1, Moon
    r_star_M_1[1] = ((R31[1,1]*X_ME_τ_1) + (R31[1,2]*Y_ME_τ_1)) + (R31[1,3]*Z_ME_τ_1)
    r_star_M_1[2] = ((R31[2,1]*X_ME_τ_1) + (R31[2,2]*Y_ME_τ_1)) + (R31[2,3]*Z_ME_τ_1)
    r_star_M_1[3] = ((R31[3,1]*X_ME_τ_1) + (R31[3,2]*Y_ME_τ_1)) + (R31[3,3]*Z_ME_τ_1)
    # r-star 2, Moon
    r_star_M_2[1] = ((R32[1,1]*X_ME_τ_2) + (R32[1,2]*Y_ME_τ_2)) + (R32[1,3]*Z_ME_τ_2)
    r_star_M_2[2] = ((R32[2,1]*X_ME_τ_2) + (R32[2,2]*Y_ME_τ_2)) + (R32[2,3]*Z_ME_τ_2)
    r_star_M_2[3] = ((R32[3,1]*X_ME_τ_2) + (R32[3,2]*Y_ME_τ_2)) + (R32[3,3]*Z_ME_τ_2)
    # r-star 0, Sun
    r_star_S_0[1] = ((R30[1,1]*X_SE_τ_0) + (R30[1,2]*Y_SE_τ_0)) + (R30[1,3]*Z_SE_τ_0)
    r_star_S_0[2] = ((R30[2,1]*X_SE_τ_0) + (R30[2,2]*Y_SE_τ_0)) + (R30[2,3]*Z_SE_τ_0)
    r_star_S_0[3] = ((R30[3,1]*X_SE_τ_0) + (R30[3,2]*Y_SE_τ_0)) + (R30[3,3]*Z_SE_τ_0)
    # r-star 1, Sun
    r_star_S_1[1] = ((R31[1,1]*X_SE_τ_1) + (R31[1,2]*Y_SE_τ_1)) + (R31[1,3]*Z_SE_τ_1)
    r_star_S_1[2] = ((R31[2,1]*X_SE_τ_1) + (R31[2,2]*Y_SE_τ_1)) + (R31[2,3]*Z_SE_τ_1)
    r_star_S_1[3] = ((R31[3,1]*X_SE_τ_1) + (R31[3,2]*Y_SE_τ_1)) + (R31[3,3]*Z_SE_τ_1)
    # r-star 2, Sun
    r_star_S_2[1] = ((R32[1,1]*X_SE_τ_2) + (R32[1,2]*Y_SE_τ_2)) + (R32[1,3]*Z_SE_τ_2)
    r_star_S_2[2] = ((R32[2,1]*X_SE_τ_2) + (R32[2,2]*Y_SE_τ_2)) + (R32[2,3]*Z_SE_τ_2)
    r_star_S_2[3] = ((R32[3,1]*X_SE_τ_2) + (R32[3,2]*Y_SE_τ_2)) + (R32[3,3]*Z_SE_τ_2)

    # X_bf[mo,ea] are geocentric, Earth-fixed "unprimed" position of perturbed body (Moon) in cylindrical coordinates

    ρ0s2_M = (r_star_M_0[1]^2) + (r_star_M_0[2]^2)
    ρ0s_M = sqrt(ρ0s2_M)
    z0s2_M = r_star_M_0[3]^2
    r0s2_M = ρ0s2_M + z0s2_M
    r0s_M = sqrt(r0s2_M)
    r0s5_M = r0s_M^5

    ρ0s2_S = (r_star_S_0[1]^2) + (r_star_S_0[2]^2)
    ρ0s_S = sqrt(ρ0s2_S)
    z0s2_S = r_star_S_0[3]^2
    r0s2_S = ρ0s2_S + z0s2_S
    r0s_S = sqrt(r0s2_S)
    r0s5_S = r0s_S^5

    coeff0_M = r0s2_M - 5( ( ((Z_bf[mo,ea]*r_star_M_0[3])^2) + 0.5((r_xy[mo,ea]*ρ0s_M)^2) )/r_p2[mo,ea] )
    coeff0_S = r0s2_S - 5( ( ((Z_bf[mo,ea]*r_star_S_0[3])^2) + 0.5((r_xy[mo,ea]*ρ0s_S)^2) )/r_p2[mo,ea] )

    k_20E_div_r0s5_M = k_20E/r0s5_M
    k_20E_div_r0s5_S = k_20E/r0s5_S

    aux0_M_x = k_20E_div_r0s5_M*((ρ0s2_M + coeff0_M)*X_bf[mo,ea])
    aux0_M_y = k_20E_div_r0s5_M*((ρ0s2_M + coeff0_M)*Y_bf[mo,ea])
    aux0_M_z = k_20E_div_r0s5_M*(((2z0s2_M) + coeff0_M)*Z_bf[mo,ea])
    aux0_S_x = k_20E_div_r0s5_S*((ρ0s2_S + coeff0_S)*X_bf[mo,ea])
    aux0_S_y = k_20E_div_r0s5_S*((ρ0s2_S + coeff0_S)*Y_bf[mo,ea])
    aux0_S_z = k_20E_div_r0s5_S*(((2z0s2_S) + coeff0_S)*Z_bf[mo,ea])

    ρ1s2_M = (r_star_M_1[1]^2) + (r_star_M_1[2]^2)
    ρ1s_M = sqrt(ρ1s2_M)
    z1s2_M = r_star_M_1[3]^2
    r1s2_M = ρ1s2_M + z1s2_M
    r1s_M = sqrt(r1s2_M)
    r1s5_M = r1s_M^5

    ρ1s2_S = (r_star_S_1[1]^2) + (r_star_S_1[2]^2)
    ρ1s_S = sqrt(ρ1s2_S)
    z1s2_S = r_star_S_1[3]^2
    r1s2_S = ρ1s2_S + z1s2_S
    r1s_S = sqrt(r1s2_S)
    r1s5_S = r1s_S^5

    coeff1_1_M = (X_bf[mo,ea]*r_star_M_1[1]) + (Y_bf[mo,ea]*r_star_M_1[2])
    coeff1_1_S = (X_bf[mo,ea]*r_star_S_1[1]) + (Y_bf[mo,ea]*r_star_S_1[2])
    coeff2_1_M = Z_bf[mo,ea]*r_star_M_1[3]
    coeff2_1_S = Z_bf[mo,ea]*r_star_S_1[3]

    coeff3_1_M = ((10coeff1_1_M)*coeff2_1_M)/r_p2[mo,ea]
    coeff3_1_S = ((10coeff1_1_S)*coeff2_1_S)/r_p2[mo,ea]

    k_21E_div_r1s5_M = k_21E/r1s5_M
    k_21E_div_r1s5_S = k_21E/r1s5_S

    aux1_M_x = k_21E_div_r1s5_M*((2coeff2_1_M*r_star_M_1[1]) - (coeff3_1_M*X_bf[mo,ea]))
    aux1_M_y = k_21E_div_r1s5_M*((2coeff2_1_M*r_star_M_1[2]) - (coeff3_1_M*Y_bf[mo,ea]))
    aux1_M_z = k_21E_div_r1s5_M*((2coeff1_1_M*r_star_M_1[3]) - (coeff3_1_M*Z_bf[mo,ea]))
    aux1_S_x = k_21E_div_r1s5_S*((2coeff2_1_S*r_star_S_1[1]) - (coeff3_1_S*X_bf[mo,ea]))
    aux1_S_y = k_21E_div_r1s5_S*((2coeff2_1_S*r_star_S_1[2]) - (coeff3_1_S*Y_bf[mo,ea]))
    aux1_S_z = k_21E_div_r1s5_S*((2coeff1_1_S*r_star_S_1[3]) - (coeff3_1_S*Z_bf[mo,ea]))

    ρ2s2_M = (r_star_M_2[1]^2) + (r_star_M_2[2]^2)
    ρ2s_M = sqrt(ρ2s2_M)
    z2s2_M = r_star_M_2[3]^2
    r2s2_M = ρ2s2_M + z2s2_M
    r2s_M = sqrt(r2s2_M)
    r2s5_M = r2s_M^5

    ρ2s2_S = (r_star_S_2[1]^2) + (r_star_S_2[2]^2)
    ρ2s_S = sqrt(ρ2s2_S)
    z2s2_S = r_star_S_2[3]^2
    r2s2_S = ρ2s2_S + z2s2_S
    r2s_S = sqrt(r2s2_S)
    r2s5_S = r2s_S^5

    coeff1_2_M = (X_bf[mo,ea]*r_star_M_2[1]) + (Y_bf[mo,ea]*r_star_M_2[2])
    coeff1_2_S = (X_bf[mo,ea]*r_star_S_2[1]) + (Y_bf[mo,ea]*r_star_S_2[2])

    coeff3_2_M = 5( (coeff1_2_M^2) - 0.5((r_xy[mo,ea]^2)*ρ2s2_M) )/r_p2[mo,ea]
    coeff3_2_S = 5( (coeff1_2_S^2) - 0.5((r_xy[mo,ea]^2)*ρ2s2_S) )/r_p2[mo,ea]

    k_22E_div_r2s5_M = k_22E/r2s5_M
    k_22E_div_r2s5_S = k_22E/r2s5_S

    aux2_M_x = k_22E_div_r2s5_M*(  (2coeff1_2_M*r_star_M_2[1])-(ρ2s2_M+coeff3_2_M)*X_bf[mo,ea]  )
    aux2_M_y = k_22E_div_r2s5_M*(  (2coeff1_2_M*r_star_M_2[2])-(ρ2s2_M+coeff3_2_M)*Y_bf[mo,ea]  )
    aux2_M_z = k_22E_div_r2s5_M*(  -coeff3_2_M*Z_bf[mo,ea]  )
    aux2_S_x = k_22E_div_r2s5_S*(  (2coeff1_2_S*r_star_S_2[1])-(ρ2s2_S+coeff3_2_S)*X_bf[mo,ea]  )
    aux2_S_y = k_22E_div_r2s5_S*(  (2coeff1_2_S*r_star_S_2[2])-(ρ2s2_S+coeff3_2_S)*Y_bf[mo,ea]  )
    aux2_S_z = k_22E_div_r2s5_S*(  -coeff3_2_S*Z_bf[mo,ea]  )

    RE_div_r_p5 = (RE_au/r_p1d2[mo,ea])^5
    aux_tidacc = tid_num_coeff*RE_div_r_p5
    tide_acc_coeff_M = μ[mo]*aux_tidacc
    tide_acc_coeff_S = μ[su]*aux_tidacc

    # add contributions from long-period, diurnal and semi-diurnal tides
    tidal_bf_x = (tide_acc_coeff_M*((aux0_M_x+aux1_M_x)+aux2_M_x)) + (tide_acc_coeff_S*((aux0_S_x+aux1_S_x)+aux2_S_x))
    tidal_bf_y = (tide_acc_coeff_M*((aux0_M_y+aux1_M_y)+aux2_M_y)) + (tide_acc_coeff_S*((aux0_S_y+aux1_S_y)+aux2_S_y))
    tidal_bf_z = (tide_acc_coeff_M*((aux0_M_z+aux1_M_z)+aux2_M_z)) + (tide_acc_coeff_S*((aux0_S_z+aux1_S_z)+aux2_S_z))

    # transform from geocentric Earth-true-equator-of-date coordinates to geocentric mean equator of J2000.0 coordinates
    tidal_x = ((M_[1,1,ea]*tidal_bf_x)+(M_[2,1,ea]*tidal_bf_y)) + (M_[3,1,ea]*tidal_bf_z)
    tidal_y = ((M_[1,2,ea]*tidal_bf_x)+(M_[2,2,ea]*tidal_bf_y)) + (M_[3,2,ea]*tidal_bf_z)
    tidal_z = ((M_[1,3,ea]*tidal_bf_x)+(M_[2,3,ea]*tidal_bf_y)) + (M_[3,3,ea]*tidal_bf_z)

    # add tidal acceleration to Moon's acceleration due to extended-body effects
    accX_mo_tides = accX[mo] + tidal_x
    accY_mo_tides = accY[mo] + tidal_y
    accZ_mo_tides = accZ[mo] + tidal_z
    accX[mo] = accX_mo_tides
    accY[mo] = accY_mo_tides
    accZ[mo] = accZ_mo_tides

    #fill accelerations (post-Newtonian and extended body accelerations)
    Threads.@threads for i in 1:N_ext
        dq[3(N+i)-2] = postNewtonX[i,postnewton_iter+1] + accX[i]
        dq[3(N+i)-1] = postNewtonY[i,postnewton_iter+1] + accY[i]
        dq[3(N+i)  ] = postNewtonZ[i,postnewton_iter+1] + accZ[i]
    end
    Threads.@threads for i in N_ext+1:N
        dq[3(N+i)-2] = postNewtonX[i,postnewton_iter+1]
        dq[3(N+i)-1] = postNewtonY[i,postnewton_iter+1]
        dq[3(N+i)  ] = postNewtonZ[i,postnewton_iter+1]
    end

    nothing
end

function TaylorIntegration.jetcoeffs!(::Val{NBP_pN_A_J23E_J23M_J2S!}, t::Taylor1{_T}, q::AbstractVector{Taylor1{_S}}, dq::AbstractVector{Taylor1{_S}}, params) where {_T <: Real, _S <: Number}
    order = t.order
    local (N, S, eulang_de430_, jd0) = params
    local eulang_t = eulang_de430_((t + (jd0 - J2000)) * daysec)
    local postnewton_iter = 1
    local c_p2 = 29979.063823897606
    local c_m2 = 3.3356611996764786e-5
    local zero_q_1 = zero(q[1])
    local one_t = one(t)
    X = Array{Taylor1{S}}(undef, N, N)
    Y = Array{Taylor1{S}}(undef, N, N)
    Z = Array{Taylor1{S}}(undef, N, N)
    r_p2 = Array{Taylor1{S}}(undef, N, N)
    r_p3d2 = Array{Taylor1{S}}(undef, N, N)
    r_p7d2 = Array{Taylor1{S}}(undef, N, N)
    newtonX = Array{Taylor1{S}}(undef, N)
    newtonY = Array{Taylor1{S}}(undef, N)
    newtonZ = Array{Taylor1{S}}(undef, N)
    newtonianCoeff = Array{Taylor1{S}}(undef, N, N)
    U = Array{Taylor1{S}}(undef, N, N)
    V = Array{Taylor1{S}}(undef, N, N)
    W = Array{Taylor1{S}}(undef, N, N)
    _4U_m_3X = Array{Taylor1{S}}(undef, N, N)
    _4V_m_3Y = Array{Taylor1{S}}(undef, N, N)
    _4W_m_3Z = Array{Taylor1{S}}(undef, N, N)
    UU = Array{Taylor1{S}}(undef, N, N)
    VV = Array{Taylor1{S}}(undef, N, N)
    WW = Array{Taylor1{S}}(undef, N, N)
    r_p1d2 = Array{Taylor1{S}}(undef, N, N)
    newtonianNb_Potential = Array{Taylor1{S}}(undef, N)
    newtonian1b_Potential = Array{Taylor1{S}}(undef, N, N)
    newton_acc_X = Array{Taylor1{S}}(undef, N, N)
    newton_acc_Y = Array{Taylor1{S}}(undef, N, N)
    newton_acc_Z = Array{Taylor1{S}}(undef, N, N)
    v2 = Array{Taylor1{S}}(undef, N)
    vi_dot_vj = Array{Taylor1{S}}(undef, N, N)
    pn2 = Array{Taylor1{S}}(undef, N, N)
    pn3 = Array{Taylor1{S}}(undef, N, N)
    _4ϕj = Array{Taylor1{S}}(undef, N, N)
    ϕi_plus_4ϕj = Array{Taylor1{S}}(undef, N, N)
    sj2_plus_2si2_minus_4vivj = Array{Taylor1{S}}(undef, N, N)
    ϕs_and_vs = Array{Taylor1{S}}(undef, N, N)
    U_t_pn2 = Array{Taylor1{S}}(undef, N, N)
    V_t_pn2 = Array{Taylor1{S}}(undef, N, N)
    W_t_pn2 = Array{Taylor1{S}}(undef, N, N)
    pn1t1_7 = Array{Taylor1{S}}(undef, N, N)
    pntempX = Array{Taylor1{S}}(undef, N, postnewton_iter)
    pntempY = Array{Taylor1{S}}(undef, N, postnewton_iter)
    pntempZ = Array{Taylor1{S}}(undef, N, postnewton_iter)
    pn1 = Array{Taylor1{S}}(undef, N, N, postnewton_iter)
    X_t_pn1 = Array{Taylor1{S}}(undef, N, N, postnewton_iter)
    Y_t_pn1 = Array{Taylor1{S}}(undef, N, N, postnewton_iter)
    Z_t_pn1 = Array{Taylor1{S}}(undef, N, N, postnewton_iter)
    pNX_t_pn3 = Array{Taylor1{S}}(undef, N, N, postnewton_iter)
    pNY_t_pn3 = Array{Taylor1{S}}(undef, N, N, postnewton_iter)
    pNZ_t_pn3 = Array{Taylor1{S}}(undef, N, N, postnewton_iter)
    pNX_t_X = Array{Taylor1{S}}(undef, N, N, postnewton_iter)
    pNY_t_Y = Array{Taylor1{S}}(undef, N, N, postnewton_iter)
    pNZ_t_Z = Array{Taylor1{S}}(undef, N, N, postnewton_iter)
    postNewtonX = Array{Taylor1{S}}(undef, N, postnewton_iter + 1)
    postNewtonY = Array{Taylor1{S}}(undef, N, postnewton_iter + 1)
    postNewtonZ = Array{Taylor1{S}}(undef, N, postnewton_iter + 1)
    X_bf_1 = Array{Taylor1{S}}(undef, N, N)
    Y_bf_1 = Array{Taylor1{S}}(undef, N, N)
    Z_bf_1 = Array{Taylor1{S}}(undef, N, N)
    X_bf_2 = Array{Taylor1{S}}(undef, N, N)
    Y_bf_2 = Array{Taylor1{S}}(undef, N, N)
    Z_bf_2 = Array{Taylor1{S}}(undef, N, N)
    X_bf_3 = Array{Taylor1{S}}(undef, N, N)
    Y_bf_3 = Array{Taylor1{S}}(undef, N, N)
    Z_bf_3 = Array{Taylor1{S}}(undef, N, N)
    X_bf = Array{Taylor1{S}}(undef, N, N)
    Y_bf = Array{Taylor1{S}}(undef, N, N)
    Z_bf = Array{Taylor1{S}}(undef, N, N)
    F_JCS_x = Array{Taylor1{S}}(undef, N, N)
    F_JCS_y = Array{Taylor1{S}}(undef, N, N)
    F_JCS_z = Array{Taylor1{S}}(undef, N, N)
    temp_accX_j = Array{Taylor1{S}}(undef, N, N)
    temp_accY_j = Array{Taylor1{S}}(undef, N, N)
    temp_accZ_j = Array{Taylor1{S}}(undef, N, N)
    temp_accX_i = Array{Taylor1{S}}(undef, N, N)
    temp_accY_i = Array{Taylor1{S}}(undef, N, N)
    temp_accZ_i = Array{Taylor1{S}}(undef, N, N)
    sin_ϕ = Array{Taylor1{S}}(undef, N, N)
    cos_ϕ = Array{Taylor1{S}}(undef, N, N)
    sin_λ = Array{Taylor1{S}}(undef, N, N)
    cos_λ = Array{Taylor1{S}}(undef, N, N)
    r_xy = Array{Taylor1{S}}(undef, N, N)
    r_p4 = Array{Taylor1{S}}(undef, N, N)
    P_n = Array{Taylor1{S}}(undef, N, N, maximum(n1SEM) + 1)
    dP_n = Array{Taylor1{S}}(undef, N, N, maximum(n1SEM) + 1)
    temp_fjξ = Array{Taylor1{S}}(undef, N, N, maximum(n1SEM) + 1)
    temp_fjζ = Array{Taylor1{S}}(undef, N, N, maximum(n1SEM) + 1)
    temp_rn = Array{Taylor1{S}}(undef, N, N, maximum(n1SEM) + 1)
    F_CS_ξ_36 = Array{Taylor1{S}}(undef, N, N)
    F_CS_η_36 = Array{Taylor1{S}}(undef, N, N)
    F_CS_ζ_36 = Array{Taylor1{S}}(undef, N, N)
    F_J_ξ_36 = Array{Taylor1{S}}(undef, N, N)
    F_J_ζ_36 = Array{Taylor1{S}}(undef, N, N)
    sin_mλ = Array{Taylor1{S}}(undef, N, N, n1SEM[mo])
    cos_mλ = Array{Taylor1{S}}(undef, N, N, n1SEM[mo])
    secϕ_P_nm = Array{Taylor1{S}}(undef, N, N, n1SEM[mo] + 1, n1SEM[mo] + 1)
    P_nm = Array{Taylor1{S}}(undef, N, N, n1SEM[mo] + 1, n1SEM[mo] + 1)
    cosϕ_dP_nm = Array{Taylor1{S}}(undef, N, N, n1SEM[mo] + 1, n1SEM[mo] + 1)
    F_J_ξ = Array{Taylor1{S}}(undef, N, N)
    F_J_η = Array{Taylor1{S}}(undef, N, N)
    F_J_ζ = Array{Taylor1{S}}(undef, N, N)
    F_CS_ξ = Array{Taylor1{S}}(undef, N, N)
    F_CS_η = Array{Taylor1{S}}(undef, N, N)
    F_CS_ζ = Array{Taylor1{S}}(undef, N, N)
    F_JCS_ξ = Array{Taylor1{S}}(undef, N, N)
    F_JCS_η = Array{Taylor1{S}}(undef, N, N)
    F_JCS_ζ = Array{Taylor1{S}}(undef, N, N)
    Rb2p = Array{Taylor1{S}}(undef, N, N, 3, 3)
    Gc2p = Array{Taylor1{S}}(undef, N, N, 3, 3)
    accX = Array{Taylor1{S}}(undef, N)
    accY = Array{Taylor1{S}}(undef, N)
    accZ = Array{Taylor1{S}}(undef, N)
    local dsj2k = t + (jd0 - 2.451545e6)
    local αs = deg2rad(α_p_sun * one_t)
    local δs = deg2rad(δ_p_sun * one_t)
    local αm = eulang_t[1] - pi / 2
    local δm = pi / 2 - eulang_t[2]
    local Wm = eulang_t[3]
    local M_ = Array{Taylor1{S}}(undef, 3, 3, N)
    local M_[:, :, ea] = c2t_jpl_de430(dsj2k)
    local M_[:, :, su] = pole_rotation(αs, δs)
    local M_[:, :, mo] = pole_rotation(αm, δm, Wm)
    local ITM_t = ITM_und .* one_t
    local fact_num = -4.5257273867882326e-36
    local fact1_jsem = [(2n - 1) / n for n = 1:maximum(n1SEM)]
    local fact2_jsem = [(n - 1) / n for n = 1:maximum(n1SEM)]
    local fact3_jsem = [n for n = 1:maximum(n1SEM)]
    local fact4_jsem = [n + 1 for n = 1:maximum(n1SEM)]
    local fact5_jsem = [n + 2 for n = 1:maximum(n1SEM)]
    local lnm1 = [(2n - 1) / (n - m) for n = 1:6, m = 1:6]
    local lnm2 = [-(((n + m) - 1)) / (n - m) for n = 1:6, m = 1:6]
    local lnm3 = [-n for n = 1:6]
    local lnm4 = [n + m for n = 1:6, m = 1:6]
    local lnm5 = [2n - 1 for n = 1:6]
    local lnm6 = [-((n + 1)) for n = 1:6]
    local lnm7 = [m for m = 1:6]
    local J2E_t = (J2E + J2EDOT * (dsj2k / yr)) * (RE / au) ^ 2
    local J2S_t = JSEM[su, 2] * one_t
    for j = 1:N
        newtonX[j] = Taylor1(identity(constant_term(zero_q_1)), order)
        newtonY[j] = Taylor1(identity(constant_term(zero_q_1)), order)
        newtonZ[j] = Taylor1(identity(constant_term(zero_q_1)), order)
        newtonianNb_Potential[j] = Taylor1(identity(constant_term(zero_q_1)), order)
        accX[j] = Taylor1(identity(constant_term(zero_q_1)), order)
        accY[j] = Taylor1(identity(constant_term(zero_q_1)), order)
        accZ[j] = Taylor1(identity(constant_term(zero_q_1)), order)
        dq[3j - 2] = Taylor1(identity(constant_term(q[3 * (N + j) - 2])), order)
        dq[3j - 1] = Taylor1(identity(constant_term(q[3 * (N + j) - 1])), order)
        dq[3j] = Taylor1(identity(constant_term(q[3 * (N + j)])), order)
    end
    tmp894 = Array{Taylor1{_S}}(undef, size(dq))
    tmp894 .= Taylor1(zero(_S), order)
    tmp896 = Array{Taylor1{_S}}(undef, size(dq))
    tmp896 .= Taylor1(zero(_S), order)
    tmp899 = Array{Taylor1{_S}}(undef, size(dq))
    tmp899 .= Taylor1(zero(_S), order)
    tmp901 = Array{Taylor1{_S}}(undef, size(dq))
    tmp901 .= Taylor1(zero(_S), order)
    tmp904 = Array{Taylor1{_S}}(undef, size(dq))
    tmp904 .= Taylor1(zero(_S), order)
    tmp906 = Array{Taylor1{_S}}(undef, size(dq))
    tmp906 .= Taylor1(zero(_S), order)
    pn2x = Array{Taylor1{_S}}(undef, size(X))
    pn2x .= Taylor1(zero(_S), order)
    pn2y = Array{Taylor1{_S}}(undef, size(Y))
    pn2y .= Taylor1(zero(_S), order)
    pn2z = Array{Taylor1{_S}}(undef, size(Z))
    pn2z .= Taylor1(zero(_S), order)
    tmp914 = Array{Taylor1{_S}}(undef, size(UU))
    tmp914 .= Taylor1(zero(_S), order)
    tmp917 = Array{Taylor1{_S}}(undef, size(X))
    tmp917 .= Taylor1(zero(_S), order)
    tmp919 = Array{Taylor1{_S}}(undef, size(Y))
    tmp919 .= Taylor1(zero(_S), order)
    tmp920 = Array{Taylor1{_S}}(undef, size(tmp917))
    tmp920 .= Taylor1(zero(_S), order)
    tmp922 = Array{Taylor1{_S}}(undef, size(Z))
    tmp922 .= Taylor1(zero(_S), order)
    tmp930 = Array{Taylor1{_S}}(undef, size(pn2x))
    tmp930 .= Taylor1(zero(_S), order)
    tmp931 = Array{Taylor1{_S}}(undef, size(tmp930))
    tmp931 .= Taylor1(zero(_S), order)
    tmp942 = Array{Taylor1{_S}}(undef, size(X))
    tmp942 .= Taylor1(zero(_S), order)
    temp_001 = Array{Taylor1{_S}}(undef, size(tmp942))
    temp_001 .= Taylor1(zero(_S), order)
    tmp944 = Array{Taylor1{_S}}(undef, size(Y))
    tmp944 .= Taylor1(zero(_S), order)
    temp_002 = Array{Taylor1{_S}}(undef, size(tmp944))
    temp_002 .= Taylor1(zero(_S), order)
    tmp946 = Array{Taylor1{_S}}(undef, size(Z))
    tmp946 .= Taylor1(zero(_S), order)
    temp_003 = Array{Taylor1{_S}}(undef, size(tmp946))
    temp_003 .= Taylor1(zero(_S), order)
    temp_004 = Array{Taylor1{_S}}(undef, size(newtonian1b_Potential))
    temp_004 .= Taylor1(zero(_S), order)
    tmp950 = Array{Taylor1{_S}}(undef, size(dq))
    tmp950 .= Taylor1(zero(_S), order)
    tmp952 = Array{Taylor1{_S}}(undef, size(dq))
    tmp952 .= Taylor1(zero(_S), order)
    tmp953 = Array{Taylor1{_S}}(undef, size(tmp950))
    tmp953 .= Taylor1(zero(_S), order)
    tmp955 = Array{Taylor1{_S}}(undef, size(dq))
    tmp955 .= Taylor1(zero(_S), order)
    for j = 1:N
        for i = 1:N
            if i == j
                continue
            else
                X[i, j] = Taylor1(constant_term(q[3i - 2]) - constant_term(q[3j - 2]), order)
                Y[i, j] = Taylor1(constant_term(q[3i - 1]) - constant_term(q[3j - 1]), order)
                Z[i, j] = Taylor1(constant_term(q[3i]) - constant_term(q[3j]), order)
                U[i, j] = Taylor1(constant_term(dq[3i - 2]) - constant_term(dq[3j - 2]), order)
                V[i, j] = Taylor1(constant_term(dq[3i - 1]) - constant_term(dq[3j - 1]), order)
                W[i, j] = Taylor1(constant_term(dq[3i]) - constant_term(dq[3j]), order)
                tmp894[3j - 2] = Taylor1(constant_term(4) * constant_term(dq[3j - 2]), order)
                tmp896[3i - 2] = Taylor1(constant_term(3) * constant_term(dq[3i - 2]), order)
                _4U_m_3X[i, j] = Taylor1(constant_term(tmp894[3j - 2]) - constant_term(tmp896[3i - 2]), order)
                tmp899[3j - 1] = Taylor1(constant_term(4) * constant_term(dq[3j - 1]), order)
                tmp901[3i - 1] = Taylor1(constant_term(3) * constant_term(dq[3i - 1]), order)
                _4V_m_3Y[i, j] = Taylor1(constant_term(tmp899[3j - 1]) - constant_term(tmp901[3i - 1]), order)
                tmp904[3j] = Taylor1(constant_term(4) * constant_term(dq[3j]), order)
                tmp906[3i] = Taylor1(constant_term(3) * constant_term(dq[3i]), order)
                _4W_m_3Z[i, j] = Taylor1(constant_term(tmp904[3j]) - constant_term(tmp906[3i]), order)
                pn2x[i, j] = Taylor1(constant_term(X[i, j]) * constant_term(_4U_m_3X[i, j]), order)
                pn2y[i, j] = Taylor1(constant_term(Y[i, j]) * constant_term(_4V_m_3Y[i, j]), order)
                pn2z[i, j] = Taylor1(constant_term(Z[i, j]) * constant_term(_4W_m_3Z[i, j]), order)
                UU[i, j] = Taylor1(constant_term(dq[3i - 2]) * constant_term(dq[3j - 2]), order)
                VV[i, j] = Taylor1(constant_term(dq[3i - 1]) * constant_term(dq[3j - 1]), order)
                WW[i, j] = Taylor1(constant_term(dq[3i]) * constant_term(dq[3j]), order)
                tmp914[i, j] = Taylor1(constant_term(UU[i, j]) + constant_term(VV[i, j]), order)
                vi_dot_vj[i, j] = Taylor1(constant_term(tmp914[i, j]) + constant_term(WW[i, j]), order)
                tmp917[i, j] = Taylor1(constant_term(X[i, j]) ^ constant_term(2), order)
                tmp919[i, j] = Taylor1(constant_term(Y[i, j]) ^ constant_term(2), order)
                tmp920[i, j] = Taylor1(constant_term(tmp917[i, j]) + constant_term(tmp919[i, j]), order)
                tmp922[i, j] = Taylor1(constant_term(Z[i, j]) ^ constant_term(2), order)
                r_p2[i, j] = Taylor1(constant_term(tmp920[i, j]) + constant_term(tmp922[i, j]), order)
                r_p1d2[i, j] = Taylor1(sqrt(constant_term(r_p2[i, j])), order)
                r_p3d2[i, j] = Taylor1(constant_term(r_p2[i, j]) ^ constant_term(1.5), order)
                r_p7d2[i, j] = Taylor1(constant_term(r_p2[i, j]) ^ constant_term(3.5), order)
                newtonianCoeff[i, j] = Taylor1(constant_term(μ[i]) / constant_term(r_p3d2[i, j]), order)
                tmp930[i, j] = Taylor1(constant_term(pn2x[i, j]) + constant_term(pn2y[i, j]), order)
                tmp931[i, j] = Taylor1(constant_term(tmp930[i, j]) + constant_term(pn2z[i, j]), order)
                pn2[i, j] = Taylor1(constant_term(newtonianCoeff[i, j]) * constant_term(tmp931[i, j]), order)
                newton_acc_X[i, j] = Taylor1(constant_term(X[i, j]) * constant_term(newtonianCoeff[i, j]), order)
                newton_acc_Y[i, j] = Taylor1(constant_term(Y[i, j]) * constant_term(newtonianCoeff[i, j]), order)
                newton_acc_Z[i, j] = Taylor1(constant_term(Z[i, j]) * constant_term(newtonianCoeff[i, j]), order)
                newtonian1b_Potential[i, j] = Taylor1(constant_term(μ[i]) / constant_term(r_p1d2[i, j]), order)
                pn3[i, j] = Taylor1(constant_term(3.5) * constant_term(newtonian1b_Potential[i, j]), order)
                U_t_pn2[i, j] = Taylor1(constant_term(pn2[i, j]) * constant_term(U[i, j]), order)
                V_t_pn2[i, j] = Taylor1(constant_term(pn2[i, j]) * constant_term(V[i, j]), order)
                W_t_pn2[i, j] = Taylor1(constant_term(pn2[i, j]) * constant_term(W[i, j]), order)
                tmp942[i, j] = Taylor1(constant_term(X[i, j]) * constant_term(newtonianCoeff[i, j]), order)
                temp_001[i, j] = Taylor1(constant_term(newtonX[j]) + constant_term(tmp942[i, j]), order)
                newtonX[j] = Taylor1(identity(constant_term(temp_001[i, j])), order)
                tmp944[i, j] = Taylor1(constant_term(Y[i, j]) * constant_term(newtonianCoeff[i, j]), order)
                temp_002[i, j] = Taylor1(constant_term(newtonY[j]) + constant_term(tmp944[i, j]), order)
                newtonY[j] = Taylor1(identity(constant_term(temp_002[i, j])), order)
                tmp946[i, j] = Taylor1(constant_term(Z[i, j]) * constant_term(newtonianCoeff[i, j]), order)
                temp_003[i, j] = Taylor1(constant_term(newtonZ[j]) + constant_term(tmp946[i, j]), order)
                newtonZ[j] = Taylor1(identity(constant_term(temp_003[i, j])), order)
                temp_004[i, j] = Taylor1(constant_term(newtonianNb_Potential[j]) + constant_term(newtonian1b_Potential[i, j]), order)
                newtonianNb_Potential[j] = Taylor1(identity(constant_term(temp_004[i, j])), order)
            end
        end
        tmp950[3j - 2] = Taylor1(constant_term(dq[3j - 2]) ^ constant_term(2), order)
        tmp952[3j - 1] = Taylor1(constant_term(dq[3j - 1]) ^ constant_term(2), order)
        tmp953[3j - 2] = Taylor1(constant_term(tmp950[3j - 2]) + constant_term(tmp952[3j - 1]), order)
        tmp955[3j] = Taylor1(constant_term(dq[3j]) ^ constant_term(2), order)
        v2[j] = Taylor1(constant_term(tmp953[3j - 2]) + constant_term(tmp955[3j]), order)
    end
    tmp957 = Taylor1(constant_term(ITM_t[1, 1]) + constant_term(ITM_t[2, 2]), order)
    tmp959 = Taylor1(constant_term(tmp957) / constant_term(2), order)
    tmp960 = Taylor1(constant_term(ITM_t[3, 3]) - constant_term(tmp959), order)
    J2M_t = Taylor1(constant_term(tmp960) / constant_term(μ[mo]), order)
    tmp962 = Taylor1(constant_term(ITM_t[2, 2]) - constant_term(ITM_t[1, 1]), order)
    tmp963 = Taylor1(constant_term(tmp962) / constant_term(μ[mo]), order)
    C22M_t = Taylor1(constant_term(tmp963) / constant_term(4), order)
    tmp966 = Taylor1(-(constant_term(ITM_t[1, 3])), order)
    C21M_t = Taylor1(constant_term(tmp966) / constant_term(μ[mo]), order)
    tmp968 = Taylor1(-(constant_term(ITM_t[3, 2])), order)
    S21M_t = Taylor1(constant_term(tmp968) / constant_term(μ[mo]), order)
    tmp970 = Taylor1(-(constant_term(ITM_t[2, 1])), order)
    tmp971 = Taylor1(constant_term(tmp970) / constant_term(μ[mo]), order)
    S22M_t = Taylor1(constant_term(tmp971) / constant_term(2), order)
    tmp983 = Array{Taylor1{_S}}(undef, size(X_bf_1))
    tmp983 .= Taylor1(zero(_S), order)
    tmp985 = Array{Taylor1{_S}}(undef, size(Y_bf_1))
    tmp985 .= Taylor1(zero(_S), order)
    tmp987 = Array{Taylor1{_S}}(undef, size(Z_bf_1))
    tmp987 .= Taylor1(zero(_S), order)
    tmp991 = Array{Taylor1{_S}}(undef, size(X_bf))
    tmp991 .= Taylor1(zero(_S), order)
    tmp993 = Array{Taylor1{_S}}(undef, size(Y_bf))
    tmp993 .= Taylor1(zero(_S), order)
    tmp994 = Array{Taylor1{_S}}(undef, size(tmp991))
    tmp994 .= Taylor1(zero(_S), order)
    tmp999 = Array{Taylor1{_S}}(undef, size(P_n))
    tmp999 .= Taylor1(zero(_S), order)
    tmp1000 = Array{Taylor1{_S}}(undef, size(tmp999))
    tmp1000 .= Taylor1(zero(_S), order)
    tmp1001 = Array{Taylor1{_S}}(undef, size(P_n))
    tmp1001 .= Taylor1(zero(_S), order)
    tmp1003 = Array{Taylor1{_S}}(undef, size(dP_n))
    tmp1003 .= Taylor1(zero(_S), order)
    tmp1004 = Array{Taylor1{_S}}(undef, size(P_n))
    tmp1004 .= Taylor1(zero(_S), order)
    tmp1009 = Array{Taylor1{_S}}(undef, size(P_n))
    tmp1009 .= Taylor1(zero(_S), order)
    tmp1010 = Array{Taylor1{_S}}(undef, size(tmp1009))
    tmp1010 .= Taylor1(zero(_S), order)
    tmp1012 = Array{Taylor1{_S}}(undef, size(dP_n))
    tmp1012 .= Taylor1(zero(_S), order)
    tmp1013 = Array{Taylor1{_S}}(undef, size(tmp1012))
    tmp1013 .= Taylor1(zero(_S), order)
    tmp1014 = Array{Taylor1{_S}}(undef, size(tmp1013))
    tmp1014 .= Taylor1(zero(_S), order)
    tmp1016 = Array{Taylor1{_S}}(undef, size(P_n))
    tmp1016 .= Taylor1(zero(_S), order)
    tmp1017 = Array{Taylor1{_S}}(undef, size(tmp1016))
    tmp1017 .= Taylor1(zero(_S), order)
    tmp1019 = Array{Taylor1{_S}}(undef, size(dP_n))
    tmp1019 .= Taylor1(zero(_S), order)
    tmp1020 = Array{Taylor1{_S}}(undef, size(tmp1019))
    tmp1020 .= Taylor1(zero(_S), order)
    tmp1021 = Array{Taylor1{_S}}(undef, size(tmp1020))
    tmp1021 .= Taylor1(zero(_S), order)
    tmp1023 = Array{Taylor1{_S}}(undef, size(P_n))
    tmp1023 .= Taylor1(zero(_S), order)
    tmp1024 = Array{Taylor1{_S}}(undef, size(tmp1023))
    tmp1024 .= Taylor1(zero(_S), order)
    tmp1026 = Array{Taylor1{_S}}(undef, size(dP_n))
    tmp1026 .= Taylor1(zero(_S), order)
    tmp1027 = Array{Taylor1{_S}}(undef, size(tmp1026))
    tmp1027 .= Taylor1(zero(_S), order)
    tmp1028 = Array{Taylor1{_S}}(undef, size(tmp1027))
    tmp1028 .= Taylor1(zero(_S), order)
    tmp1030 = Array{Taylor1{_S}}(undef, size(P_n))
    tmp1030 .= Taylor1(zero(_S), order)
    tmp1031 = Array{Taylor1{_S}}(undef, size(tmp1030))
    tmp1031 .= Taylor1(zero(_S), order)
    tmp1032 = Array{Taylor1{_S}}(undef, size(tmp1031))
    tmp1032 .= Taylor1(zero(_S), order)
    tmp1034 = Array{Taylor1{_S}}(undef, size(dP_n))
    tmp1034 .= Taylor1(zero(_S), order)
    tmp1035 = Array{Taylor1{_S}}(undef, size(tmp1034))
    tmp1035 .= Taylor1(zero(_S), order)
    tmp1036 = Array{Taylor1{_S}}(undef, size(tmp1035))
    tmp1036 .= Taylor1(zero(_S), order)
    tmp1037 = Array{Taylor1{_S}}(undef, size(tmp1036))
    tmp1037 .= Taylor1(zero(_S), order)
    tmp1039 = Array{Taylor1{_S}}(undef, size(sin_mλ))
    tmp1039 .= Taylor1(zero(_S), order)
    tmp1040 = Array{Taylor1{_S}}(undef, size(cos_mλ))
    tmp1040 .= Taylor1(zero(_S), order)
    tmp1042 = Array{Taylor1{_S}}(undef, size(cos_mλ))
    tmp1042 .= Taylor1(zero(_S), order)
    tmp1043 = Array{Taylor1{_S}}(undef, size(sin_mλ))
    tmp1043 .= Taylor1(zero(_S), order)
    tmp1045 = Array{Taylor1{_S}}(undef, size(secϕ_P_nm))
    tmp1045 .= Taylor1(zero(_S), order)
    tmp1048 = Array{Taylor1{_S}}(undef, size(secϕ_P_nm))
    tmp1048 .= Taylor1(zero(_S), order)
    tmp1050 = Array{Taylor1{_S}}(undef, size(secϕ_P_nm))
    tmp1050 .= Taylor1(zero(_S), order)
    tmp1052 = Array{Taylor1{_S}}(undef, size(secϕ_P_nm))
    tmp1052 .= Taylor1(zero(_S), order)
    tmp1053 = Array{Taylor1{_S}}(undef, size(tmp1052))
    tmp1053 .= Taylor1(zero(_S), order)
    tmp1054 = Array{Taylor1{_S}}(undef, size(secϕ_P_nm))
    tmp1054 .= Taylor1(zero(_S), order)
    tmp1057 = Array{Taylor1{_S}}(undef, size(secϕ_P_nm))
    tmp1057 .= Taylor1(zero(_S), order)
    tmp1058 = Array{Taylor1{_S}}(undef, size(tmp1057))
    tmp1058 .= Taylor1(zero(_S), order)
    tmp1059 = Array{Taylor1{_S}}(undef, size(secϕ_P_nm))
    tmp1059 .= Taylor1(zero(_S), order)
    tmp1061 = Array{Taylor1{_S}}(undef, size(P_nm))
    tmp1061 .= Taylor1(zero(_S), order)
    tmp1062 = Array{Taylor1{_S}}(undef, size(cos_mλ))
    tmp1062 .= Taylor1(zero(_S), order)
    tmp1063 = Array{Taylor1{_S}}(undef, size(sin_mλ))
    tmp1063 .= Taylor1(zero(_S), order)
    tmp1064 = Array{Taylor1{_S}}(undef, size(tmp1062))
    tmp1064 .= Taylor1(zero(_S), order)
    tmp1065 = Array{Taylor1{_S}}(undef, size(tmp1061))
    tmp1065 .= Taylor1(zero(_S), order)
    tmp1066 = Array{Taylor1{_S}}(undef, size(P_nm))
    tmp1066 .= Taylor1(zero(_S), order)
    tmp1067 = Array{Taylor1{_S}}(undef, size(cos_mλ))
    tmp1067 .= Taylor1(zero(_S), order)
    tmp1068 = Array{Taylor1{_S}}(undef, size(sin_mλ))
    tmp1068 .= Taylor1(zero(_S), order)
    tmp1069 = Array{Taylor1{_S}}(undef, size(tmp1067))
    tmp1069 .= Taylor1(zero(_S), order)
    tmp1070 = Array{Taylor1{_S}}(undef, size(tmp1066))
    tmp1070 .= Taylor1(zero(_S), order)
    tmp1071 = Array{Taylor1{_S}}(undef, size(tmp1065))
    tmp1071 .= Taylor1(zero(_S), order)
    tmp1073 = Array{Taylor1{_S}}(undef, size(secϕ_P_nm))
    tmp1073 .= Taylor1(zero(_S), order)
    tmp1074 = Array{Taylor1{_S}}(undef, size(cos_mλ))
    tmp1074 .= Taylor1(zero(_S), order)
    tmp1075 = Array{Taylor1{_S}}(undef, size(sin_mλ))
    tmp1075 .= Taylor1(zero(_S), order)
    tmp1076 = Array{Taylor1{_S}}(undef, size(tmp1074))
    tmp1076 .= Taylor1(zero(_S), order)
    tmp1077 = Array{Taylor1{_S}}(undef, size(tmp1073))
    tmp1077 .= Taylor1(zero(_S), order)
    tmp1078 = Array{Taylor1{_S}}(undef, size(secϕ_P_nm))
    tmp1078 .= Taylor1(zero(_S), order)
    tmp1079 = Array{Taylor1{_S}}(undef, size(cos_mλ))
    tmp1079 .= Taylor1(zero(_S), order)
    tmp1080 = Array{Taylor1{_S}}(undef, size(sin_mλ))
    tmp1080 .= Taylor1(zero(_S), order)
    tmp1081 = Array{Taylor1{_S}}(undef, size(tmp1079))
    tmp1081 .= Taylor1(zero(_S), order)
    tmp1082 = Array{Taylor1{_S}}(undef, size(tmp1078))
    tmp1082 .= Taylor1(zero(_S), order)
    tmp1083 = Array{Taylor1{_S}}(undef, size(tmp1077))
    tmp1083 .= Taylor1(zero(_S), order)
    tmp1085 = Array{Taylor1{_S}}(undef, size(cos_mλ))
    tmp1085 .= Taylor1(zero(_S), order)
    tmp1086 = Array{Taylor1{_S}}(undef, size(sin_mλ))
    tmp1086 .= Taylor1(zero(_S), order)
    tmp1087 = Array{Taylor1{_S}}(undef, size(tmp1085))
    tmp1087 .= Taylor1(zero(_S), order)
    tmp1088 = Array{Taylor1{_S}}(undef, size(cosϕ_dP_nm))
    tmp1088 .= Taylor1(zero(_S), order)
    tmp1089 = Array{Taylor1{_S}}(undef, size(cos_mλ))
    tmp1089 .= Taylor1(zero(_S), order)
    tmp1090 = Array{Taylor1{_S}}(undef, size(sin_mλ))
    tmp1090 .= Taylor1(zero(_S), order)
    tmp1091 = Array{Taylor1{_S}}(undef, size(tmp1089))
    tmp1091 .= Taylor1(zero(_S), order)
    tmp1092 = Array{Taylor1{_S}}(undef, size(cosϕ_dP_nm))
    tmp1092 .= Taylor1(zero(_S), order)
    tmp1093 = Array{Taylor1{_S}}(undef, size(tmp1088))
    tmp1093 .= Taylor1(zero(_S), order)
    tmp1095 = Array{Taylor1{_S}}(undef, size(P_nm))
    tmp1095 .= Taylor1(zero(_S), order)
    tmp1096 = Array{Taylor1{_S}}(undef, size(cos_mλ))
    tmp1096 .= Taylor1(zero(_S), order)
    tmp1097 = Array{Taylor1{_S}}(undef, size(sin_mλ))
    tmp1097 .= Taylor1(zero(_S), order)
    tmp1098 = Array{Taylor1{_S}}(undef, size(tmp1096))
    tmp1098 .= Taylor1(zero(_S), order)
    tmp1099 = Array{Taylor1{_S}}(undef, size(tmp1095))
    tmp1099 .= Taylor1(zero(_S), order)
    tmp1100 = Array{Taylor1{_S}}(undef, size(tmp1099))
    tmp1100 .= Taylor1(zero(_S), order)
    temp_CS_ξ = Array{Taylor1{_S}}(undef, size(tmp1100))
    temp_CS_ξ .= Taylor1(zero(_S), order)
    tmp1102 = Array{Taylor1{_S}}(undef, size(secϕ_P_nm))
    tmp1102 .= Taylor1(zero(_S), order)
    tmp1103 = Array{Taylor1{_S}}(undef, size(cos_mλ))
    tmp1103 .= Taylor1(zero(_S), order)
    tmp1104 = Array{Taylor1{_S}}(undef, size(sin_mλ))
    tmp1104 .= Taylor1(zero(_S), order)
    tmp1105 = Array{Taylor1{_S}}(undef, size(tmp1103))
    tmp1105 .= Taylor1(zero(_S), order)
    tmp1106 = Array{Taylor1{_S}}(undef, size(tmp1102))
    tmp1106 .= Taylor1(zero(_S), order)
    tmp1107 = Array{Taylor1{_S}}(undef, size(tmp1106))
    tmp1107 .= Taylor1(zero(_S), order)
    temp_CS_η = Array{Taylor1{_S}}(undef, size(tmp1107))
    temp_CS_η .= Taylor1(zero(_S), order)
    tmp1109 = Array{Taylor1{_S}}(undef, size(cos_mλ))
    tmp1109 .= Taylor1(zero(_S), order)
    tmp1110 = Array{Taylor1{_S}}(undef, size(sin_mλ))
    tmp1110 .= Taylor1(zero(_S), order)
    tmp1111 = Array{Taylor1{_S}}(undef, size(tmp1109))
    tmp1111 .= Taylor1(zero(_S), order)
    tmp1112 = Array{Taylor1{_S}}(undef, size(cosϕ_dP_nm))
    tmp1112 .= Taylor1(zero(_S), order)
    tmp1113 = Array{Taylor1{_S}}(undef, size(tmp1112))
    tmp1113 .= Taylor1(zero(_S), order)
    temp_CS_ζ = Array{Taylor1{_S}}(undef, size(tmp1113))
    temp_CS_ζ .= Taylor1(zero(_S), order)
    tmp1115 = Array{Taylor1{_S}}(undef, size(F_J_ξ))
    tmp1115 .= Taylor1(zero(_S), order)
    tmp1116 = Array{Taylor1{_S}}(undef, size(F_CS_ξ))
    tmp1116 .= Taylor1(zero(_S), order)
    tmp1119 = Array{Taylor1{_S}}(undef, size(F_J_ζ))
    tmp1119 .= Taylor1(zero(_S), order)
    tmp1120 = Array{Taylor1{_S}}(undef, size(F_CS_ζ))
    tmp1120 .= Taylor1(zero(_S), order)
    tmp1126 = Array{Taylor1{_S}}(undef, size(sin_ϕ))
    tmp1126 .= Taylor1(zero(_S), order)
    tmp1129 = Array{Taylor1{_S}}(undef, size(sin_ϕ))
    tmp1129 .= Taylor1(zero(_S), order)
    tmp1131 = Array{Taylor1{_S}}(undef, size(Rb2p))
    tmp1131 .= Taylor1(zero(_S), order)
    tmp1132 = Array{Taylor1{_S}}(undef, size(Rb2p))
    tmp1132 .= Taylor1(zero(_S), order)
    tmp1133 = Array{Taylor1{_S}}(undef, size(tmp1131))
    tmp1133 .= Taylor1(zero(_S), order)
    tmp1134 = Array{Taylor1{_S}}(undef, size(Rb2p))
    tmp1134 .= Taylor1(zero(_S), order)
    tmp1136 = Array{Taylor1{_S}}(undef, size(Rb2p))
    tmp1136 .= Taylor1(zero(_S), order)
    tmp1137 = Array{Taylor1{_S}}(undef, size(Rb2p))
    tmp1137 .= Taylor1(zero(_S), order)
    tmp1138 = Array{Taylor1{_S}}(undef, size(tmp1136))
    tmp1138 .= Taylor1(zero(_S), order)
    tmp1139 = Array{Taylor1{_S}}(undef, size(Rb2p))
    tmp1139 .= Taylor1(zero(_S), order)
    tmp1141 = Array{Taylor1{_S}}(undef, size(Rb2p))
    tmp1141 .= Taylor1(zero(_S), order)
    tmp1142 = Array{Taylor1{_S}}(undef, size(Rb2p))
    tmp1142 .= Taylor1(zero(_S), order)
    tmp1143 = Array{Taylor1{_S}}(undef, size(tmp1141))
    tmp1143 .= Taylor1(zero(_S), order)
    tmp1144 = Array{Taylor1{_S}}(undef, size(Rb2p))
    tmp1144 .= Taylor1(zero(_S), order)
    tmp1146 = Array{Taylor1{_S}}(undef, size(Rb2p))
    tmp1146 .= Taylor1(zero(_S), order)
    tmp1147 = Array{Taylor1{_S}}(undef, size(Rb2p))
    tmp1147 .= Taylor1(zero(_S), order)
    tmp1148 = Array{Taylor1{_S}}(undef, size(tmp1146))
    tmp1148 .= Taylor1(zero(_S), order)
    tmp1149 = Array{Taylor1{_S}}(undef, size(Rb2p))
    tmp1149 .= Taylor1(zero(_S), order)
    tmp1151 = Array{Taylor1{_S}}(undef, size(Rb2p))
    tmp1151 .= Taylor1(zero(_S), order)
    tmp1152 = Array{Taylor1{_S}}(undef, size(Rb2p))
    tmp1152 .= Taylor1(zero(_S), order)
    tmp1153 = Array{Taylor1{_S}}(undef, size(tmp1151))
    tmp1153 .= Taylor1(zero(_S), order)
    tmp1154 = Array{Taylor1{_S}}(undef, size(Rb2p))
    tmp1154 .= Taylor1(zero(_S), order)
    tmp1156 = Array{Taylor1{_S}}(undef, size(Rb2p))
    tmp1156 .= Taylor1(zero(_S), order)
    tmp1157 = Array{Taylor1{_S}}(undef, size(Rb2p))
    tmp1157 .= Taylor1(zero(_S), order)
    tmp1158 = Array{Taylor1{_S}}(undef, size(tmp1156))
    tmp1158 .= Taylor1(zero(_S), order)
    tmp1159 = Array{Taylor1{_S}}(undef, size(Rb2p))
    tmp1159 .= Taylor1(zero(_S), order)
    tmp1161 = Array{Taylor1{_S}}(undef, size(Rb2p))
    tmp1161 .= Taylor1(zero(_S), order)
    tmp1162 = Array{Taylor1{_S}}(undef, size(Rb2p))
    tmp1162 .= Taylor1(zero(_S), order)
    tmp1163 = Array{Taylor1{_S}}(undef, size(tmp1161))
    tmp1163 .= Taylor1(zero(_S), order)
    tmp1164 = Array{Taylor1{_S}}(undef, size(Rb2p))
    tmp1164 .= Taylor1(zero(_S), order)
    tmp1166 = Array{Taylor1{_S}}(undef, size(Rb2p))
    tmp1166 .= Taylor1(zero(_S), order)
    tmp1167 = Array{Taylor1{_S}}(undef, size(Rb2p))
    tmp1167 .= Taylor1(zero(_S), order)
    tmp1168 = Array{Taylor1{_S}}(undef, size(tmp1166))
    tmp1168 .= Taylor1(zero(_S), order)
    tmp1169 = Array{Taylor1{_S}}(undef, size(Rb2p))
    tmp1169 .= Taylor1(zero(_S), order)
    tmp1171 = Array{Taylor1{_S}}(undef, size(Rb2p))
    tmp1171 .= Taylor1(zero(_S), order)
    tmp1172 = Array{Taylor1{_S}}(undef, size(Rb2p))
    tmp1172 .= Taylor1(zero(_S), order)
    tmp1173 = Array{Taylor1{_S}}(undef, size(tmp1171))
    tmp1173 .= Taylor1(zero(_S), order)
    tmp1174 = Array{Taylor1{_S}}(undef, size(Rb2p))
    tmp1174 .= Taylor1(zero(_S), order)
    tmp1176 = Array{Taylor1{_S}}(undef, size(Gc2p))
    tmp1176 .= Taylor1(zero(_S), order)
    tmp1177 = Array{Taylor1{_S}}(undef, size(Gc2p))
    tmp1177 .= Taylor1(zero(_S), order)
    tmp1178 = Array{Taylor1{_S}}(undef, size(tmp1176))
    tmp1178 .= Taylor1(zero(_S), order)
    tmp1179 = Array{Taylor1{_S}}(undef, size(Gc2p))
    tmp1179 .= Taylor1(zero(_S), order)
    tmp1181 = Array{Taylor1{_S}}(undef, size(Gc2p))
    tmp1181 .= Taylor1(zero(_S), order)
    tmp1182 = Array{Taylor1{_S}}(undef, size(Gc2p))
    tmp1182 .= Taylor1(zero(_S), order)
    tmp1183 = Array{Taylor1{_S}}(undef, size(tmp1181))
    tmp1183 .= Taylor1(zero(_S), order)
    tmp1184 = Array{Taylor1{_S}}(undef, size(Gc2p))
    tmp1184 .= Taylor1(zero(_S), order)
    tmp1186 = Array{Taylor1{_S}}(undef, size(Gc2p))
    tmp1186 .= Taylor1(zero(_S), order)
    tmp1187 = Array{Taylor1{_S}}(undef, size(Gc2p))
    tmp1187 .= Taylor1(zero(_S), order)
    tmp1188 = Array{Taylor1{_S}}(undef, size(tmp1186))
    tmp1188 .= Taylor1(zero(_S), order)
    tmp1189 = Array{Taylor1{_S}}(undef, size(Gc2p))
    tmp1189 .= Taylor1(zero(_S), order)
    for j = 1:N
        for i = 1:N
            if i == j
                continue
            else
                if UJ_interaction[i, j]
                    X_bf_1[i, j] = Taylor1(constant_term(X[i, j]) * constant_term(M_[1, 1, j]), order)
                    X_bf_2[i, j] = Taylor1(constant_term(Y[i, j]) * constant_term(M_[1, 2, j]), order)
                    X_bf_3[i, j] = Taylor1(constant_term(Z[i, j]) * constant_term(M_[1, 3, j]), order)
                    Y_bf_1[i, j] = Taylor1(constant_term(X[i, j]) * constant_term(M_[2, 1, j]), order)
                    Y_bf_2[i, j] = Taylor1(constant_term(Y[i, j]) * constant_term(M_[2, 2, j]), order)
                    Y_bf_3[i, j] = Taylor1(constant_term(Z[i, j]) * constant_term(M_[2, 3, j]), order)
                    Z_bf_1[i, j] = Taylor1(constant_term(X[i, j]) * constant_term(M_[3, 1, j]), order)
                    Z_bf_2[i, j] = Taylor1(constant_term(Y[i, j]) * constant_term(M_[3, 2, j]), order)
                    Z_bf_3[i, j] = Taylor1(constant_term(Z[i, j]) * constant_term(M_[3, 3, j]), order)
                    tmp983[i, j] = Taylor1(constant_term(X_bf_1[i, j]) + constant_term(X_bf_2[i, j]), order)
                    X_bf[i, j] = Taylor1(constant_term(tmp983[i, j]) + constant_term(X_bf_3[i, j]), order)
                    tmp985[i, j] = Taylor1(constant_term(Y_bf_1[i, j]) + constant_term(Y_bf_2[i, j]), order)
                    Y_bf[i, j] = Taylor1(constant_term(tmp985[i, j]) + constant_term(Y_bf_3[i, j]), order)
                    tmp987[i, j] = Taylor1(constant_term(Z_bf_1[i, j]) + constant_term(Z_bf_2[i, j]), order)
                    Z_bf[i, j] = Taylor1(constant_term(tmp987[i, j]) + constant_term(Z_bf_3[i, j]), order)
                    sin_ϕ[i, j] = Taylor1(constant_term(Z_bf[i, j]) / constant_term(r_p1d2[i, j]), order)
                    tmp991[i, j] = Taylor1(constant_term(X_bf[i, j]) ^ constant_term(2), order)
                    tmp993[i, j] = Taylor1(constant_term(Y_bf[i, j]) ^ constant_term(2), order)
                    tmp994[i, j] = Taylor1(constant_term(tmp991[i, j]) + constant_term(tmp993[i, j]), order)
                    r_xy[i, j] = Taylor1(sqrt(constant_term(tmp994[i, j])), order)
                    cos_ϕ[i, j] = Taylor1(constant_term(r_xy[i, j]) / constant_term(r_p1d2[i, j]), order)
                    sin_λ[i, j] = Taylor1(constant_term(Y_bf[i, j]) / constant_term(r_xy[i, j]), order)
                    cos_λ[i, j] = Taylor1(constant_term(X_bf[i, j]) / constant_term(r_xy[i, j]), order)
                    P_n[i, j, 1] = Taylor1(identity(constant_term(one_t)), order)
                    P_n[i, j, 2] = Taylor1(identity(constant_term(sin_ϕ[i, j])), order)
                    dP_n[i, j, 1] = Taylor1(identity(constant_term(zero_q_1)), order)
                    dP_n[i, j, 2] = Taylor1(identity(constant_term(one_t)), order)
                    for n = 2:n1SEM[j]
                        tmp999[i, j, n] = Taylor1(constant_term(P_n[i, j, n]) * constant_term(sin_ϕ[i, j]), order)
                        tmp1000[i, j, n] = Taylor1(constant_term(tmp999[i, j, n]) * constant_term(fact1_jsem[n]), order)
                        tmp1001[i, j, n - 1] = Taylor1(constant_term(P_n[i, j, n - 1]) * constant_term(fact2_jsem[n]), order)
                        P_n[i, j, n + 1] = Taylor1(constant_term(tmp1000[i, j, n]) - constant_term(tmp1001[i, j, n - 1]), order)
                        tmp1003[i, j, n] = Taylor1(constant_term(dP_n[i, j, n]) * constant_term(sin_ϕ[i, j]), order)
                        tmp1004[i, j, n] = Taylor1(constant_term(P_n[i, j, n]) * constant_term(fact3_jsem[n]), order)
                        dP_n[i, j, n + 1] = Taylor1(constant_term(tmp1003[i, j, n]) + constant_term(tmp1004[i, j, n]), order)
                        temp_rn[i, j, n] = Taylor1(constant_term(r_p1d2[i, j]) ^ constant_term(fact5_jsem[n]), order)
                    end
                    r_p4[i, j] = Taylor1(constant_term(r_p2[i, j]) ^ constant_term(2), order)
                    if j == mo
                        tmp1009[i, j, 3] = Taylor1(constant_term(P_n[i, j, 3]) * constant_term(fact4_jsem[2]), order)
                        tmp1010[i, j, 3] = Taylor1(constant_term(tmp1009[i, j, 3]) * constant_term(J2M_t), order)
                        F_J_ξ[i, j] = Taylor1(constant_term(tmp1010[i, j, 3]) / constant_term(r_p4[i, j]), order)
                        tmp1012[i, j, 3] = Taylor1(-(constant_term(dP_n[i, j, 3])), order)
                        tmp1013[i, j, 3] = Taylor1(constant_term(tmp1012[i, j, 3]) * constant_term(cos_ϕ[i, j]), order)
                        tmp1014[i, j, 3] = Taylor1(constant_term(tmp1013[i, j, 3]) * constant_term(J2M_t), order)
                        F_J_ζ[i, j] = Taylor1(constant_term(tmp1014[i, j, 3]) / constant_term(r_p4[i, j]), order)
                    else
                        if j == ea
                            tmp1016[i, j, 3] = Taylor1(constant_term(P_n[i, j, 3]) * constant_term(fact4_jsem[2]), order)
                            tmp1017[i, j, 3] = Taylor1(constant_term(tmp1016[i, j, 3]) * constant_term(J2E_t), order)
                            F_J_ξ[i, j] = Taylor1(constant_term(tmp1017[i, j, 3]) / constant_term(r_p4[i, j]), order)
                            tmp1019[i, j, 3] = Taylor1(-(constant_term(dP_n[i, j, 3])), order)
                            tmp1020[i, j, 3] = Taylor1(constant_term(tmp1019[i, j, 3]) * constant_term(cos_ϕ[i, j]), order)
                            tmp1021[i, j, 3] = Taylor1(constant_term(tmp1020[i, j, 3]) * constant_term(J2E_t), order)
                            F_J_ζ[i, j] = Taylor1(constant_term(tmp1021[i, j, 3]) / constant_term(r_p4[i, j]), order)
                        else
                            if j == su
                                tmp1023[i, j, 3] = Taylor1(constant_term(P_n[i, j, 3]) * constant_term(fact4_jsem[2]), order)
                                tmp1024[i, j, 3] = Taylor1(constant_term(tmp1023[i, j, 3]) * constant_term(J2S_t), order)
                                F_J_ξ[i, j] = Taylor1(constant_term(tmp1024[i, j, 3]) / constant_term(r_p4[i, j]), order)
                                tmp1026[i, j, 3] = Taylor1(-(constant_term(dP_n[i, j, 3])), order)
                                tmp1027[i, j, 3] = Taylor1(constant_term(tmp1026[i, j, 3]) * constant_term(cos_ϕ[i, j]), order)
                                tmp1028[i, j, 3] = Taylor1(constant_term(tmp1027[i, j, 3]) * constant_term(J2S_t), order)
                                F_J_ζ[i, j] = Taylor1(constant_term(tmp1028[i, j, 3]) / constant_term(r_p4[i, j]), order)
                            end
                        end
                    end
                    F_J_ξ_36[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                    F_J_ζ_36[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                    for n = 3:n1SEM[j]
                        tmp1030[i, j, n + 1] = Taylor1(constant_term(P_n[i, j, n + 1]) * constant_term(fact4_jsem[n]), order)
                        tmp1031[i, j, n + 1] = Taylor1(constant_term(tmp1030[i, j, n + 1]) * constant_term(JSEM[j, n]), order)
                        tmp1032[i, j, n + 1] = Taylor1(constant_term(tmp1031[i, j, n + 1]) / constant_term(temp_rn[i, j, n]), order)
                        temp_fjξ[i, j, n] = Taylor1(constant_term(tmp1032[i, j, n + 1]) + constant_term(F_J_ξ_36[i, j]), order)
                        tmp1034[i, j, n + 1] = Taylor1(-(constant_term(dP_n[i, j, n + 1])), order)
                        tmp1035[i, j, n + 1] = Taylor1(constant_term(tmp1034[i, j, n + 1]) * constant_term(cos_ϕ[i, j]), order)
                        tmp1036[i, j, n + 1] = Taylor1(constant_term(tmp1035[i, j, n + 1]) * constant_term(JSEM[j, n]), order)
                        tmp1037[i, j, n + 1] = Taylor1(constant_term(tmp1036[i, j, n + 1]) / constant_term(temp_rn[i, j, n]), order)
                        temp_fjζ[i, j, n] = Taylor1(constant_term(tmp1037[i, j, n + 1]) + constant_term(F_J_ζ_36[i, j]), order)
                        F_J_ξ_36[i, j] = Taylor1(identity(constant_term(temp_fjξ[i, j, n])), order)
                        F_J_ζ_36[i, j] = Taylor1(identity(constant_term(temp_fjζ[i, j, n])), order)
                    end
                    if j == mo
                        for m = 1:n1SEM[mo]
                            if m == 1
                                sin_mλ[i, j, 1] = Taylor1(identity(constant_term(sin_λ[i, j])), order)
                                cos_mλ[i, j, 1] = Taylor1(identity(constant_term(cos_λ[i, j])), order)
                                secϕ_P_nm[i, j, 1, 1] = Taylor1(identity(constant_term(one_t)), order)
                            else
                                tmp1039[i, j, 1] = Taylor1(constant_term(sin_mλ[i, j, 1]) * constant_term(cos_mλ[i, j, m - 1]), order)
                                tmp1040[i, j, 1] = Taylor1(constant_term(cos_mλ[i, j, 1]) * constant_term(sin_mλ[i, j, m - 1]), order)
                                sin_mλ[i, j, m] = Taylor1(constant_term(tmp1039[i, j, 1]) + constant_term(tmp1040[i, j, 1]), order)
                                tmp1042[i, j, 1] = Taylor1(constant_term(cos_mλ[i, j, 1]) * constant_term(cos_mλ[i, j, m - 1]), order)
                                tmp1043[i, j, 1] = Taylor1(constant_term(sin_mλ[i, j, 1]) * constant_term(sin_mλ[i, j, m - 1]), order)
                                cos_mλ[i, j, m] = Taylor1(constant_term(tmp1042[i, j, 1]) - constant_term(tmp1043[i, j, 1]), order)
                                tmp1045[i, j, m - 1, m - 1] = Taylor1(constant_term(secϕ_P_nm[i, j, m - 1, m - 1]) * constant_term(cos_ϕ[i, j]), order)
                                secϕ_P_nm[i, j, m, m] = Taylor1(constant_term(tmp1045[i, j, m - 1, m - 1]) * constant_term(lnm5[m]), order)
                                P_nm[i, j, m, m] = Taylor1(constant_term(secϕ_P_nm[i, j, m, m]) * constant_term(cos_ϕ[i, j]), order)
                                tmp1048[i, j, m, m] = Taylor1(constant_term(secϕ_P_nm[i, j, m, m]) * constant_term(sin_ϕ[i, j]), order)
                                cosϕ_dP_nm[i, j, m, m] = Taylor1(constant_term(tmp1048[i, j, m, m]) * constant_term(lnm3[m]), order)
                            end
                            for n = m + 1:n1SEM[mo]
                                if n == m + 1
                                    tmp1050[i, j, n - 1, m] = Taylor1(constant_term(secϕ_P_nm[i, j, n - 1, m]) * constant_term(sin_ϕ[i, j]), order)
                                    secϕ_P_nm[i, j, n, m] = Taylor1(constant_term(tmp1050[i, j, n - 1, m]) * constant_term(lnm1[n, m]), order)
                                else
                                    tmp1052[i, j, n - 1, m] = Taylor1(constant_term(secϕ_P_nm[i, j, n - 1, m]) * constant_term(sin_ϕ[i, j]), order)
                                    tmp1053[i, j, n - 1, m] = Taylor1(constant_term(tmp1052[i, j, n - 1, m]) * constant_term(lnm1[n, m]), order)
                                    tmp1054[i, j, n - 2, m] = Taylor1(constant_term(secϕ_P_nm[i, j, n - 2, m]) * constant_term(lnm2[n, m]), order)
                                    secϕ_P_nm[i, j, n, m] = Taylor1(constant_term(tmp1053[i, j, n - 1, m]) + constant_term(tmp1054[i, j, n - 2, m]), order)
                                end
                                P_nm[i, j, n, m] = Taylor1(constant_term(secϕ_P_nm[i, j, n, m]) * constant_term(cos_ϕ[i, j]), order)
                                tmp1057[i, j, n, m] = Taylor1(constant_term(secϕ_P_nm[i, j, n, m]) * constant_term(sin_ϕ[i, j]), order)
                                tmp1058[i, j, n, m] = Taylor1(constant_term(tmp1057[i, j, n, m]) * constant_term(lnm3[n]), order)
                                tmp1059[i, j, n - 1, m] = Taylor1(constant_term(secϕ_P_nm[i, j, n - 1, m]) * constant_term(lnm4[n, m]), order)
                                cosϕ_dP_nm[i, j, n, m] = Taylor1(constant_term(tmp1058[i, j, n, m]) + constant_term(tmp1059[i, j, n - 1, m]), order)
                            end
                        end
                        tmp1061[i, j, 2, 1] = Taylor1(constant_term(P_nm[i, j, 2, 1]) * constant_term(lnm6[2]), order)
                        tmp1062[i, j, 1] = Taylor1(constant_term(C21M_t) * constant_term(cos_mλ[i, j, 1]), order)
                        tmp1063[i, j, 1] = Taylor1(constant_term(S21M_t) * constant_term(sin_mλ[i, j, 1]), order)
                        tmp1064[i, j, 1] = Taylor1(constant_term(tmp1062[i, j, 1]) + constant_term(tmp1063[i, j, 1]), order)
                        tmp1065[i, j, 2, 1] = Taylor1(constant_term(tmp1061[i, j, 2, 1]) * constant_term(tmp1064[i, j, 1]), order)
                        tmp1066[i, j, 2, 2] = Taylor1(constant_term(P_nm[i, j, 2, 2]) * constant_term(lnm6[2]), order)
                        tmp1067[i, j, 2] = Taylor1(constant_term(C22M_t) * constant_term(cos_mλ[i, j, 2]), order)
                        tmp1068[i, j, 2] = Taylor1(constant_term(S22M_t) * constant_term(sin_mλ[i, j, 2]), order)
                        tmp1069[i, j, 2] = Taylor1(constant_term(tmp1067[i, j, 2]) + constant_term(tmp1068[i, j, 2]), order)
                        tmp1070[i, j, 2, 2] = Taylor1(constant_term(tmp1066[i, j, 2, 2]) * constant_term(tmp1069[i, j, 2]), order)
                        tmp1071[i, j, 2, 1] = Taylor1(constant_term(tmp1065[i, j, 2, 1]) + constant_term(tmp1070[i, j, 2, 2]), order)
                        F_CS_ξ[i, j] = Taylor1(constant_term(tmp1071[i, j, 2, 1]) / constant_term(r_p4[i, j]), order)
                        tmp1073[i, j, 2, 1] = Taylor1(constant_term(secϕ_P_nm[i, j, 2, 1]) * constant_term(lnm7[1]), order)
                        tmp1074[i, j, 1] = Taylor1(constant_term(S21M_t) * constant_term(cos_mλ[i, j, 1]), order)
                        tmp1075[i, j, 1] = Taylor1(constant_term(C21M_t) * constant_term(sin_mλ[i, j, 1]), order)
                        tmp1076[i, j, 1] = Taylor1(constant_term(tmp1074[i, j, 1]) - constant_term(tmp1075[i, j, 1]), order)
                        tmp1077[i, j, 2, 1] = Taylor1(constant_term(tmp1073[i, j, 2, 1]) * constant_term(tmp1076[i, j, 1]), order)
                        tmp1078[i, j, 2, 2] = Taylor1(constant_term(secϕ_P_nm[i, j, 2, 2]) * constant_term(lnm7[2]), order)
                        tmp1079[i, j, 2] = Taylor1(constant_term(S22M_t) * constant_term(cos_mλ[i, j, 2]), order)
                        tmp1080[i, j, 2] = Taylor1(constant_term(C22M_t) * constant_term(sin_mλ[i, j, 2]), order)
                        tmp1081[i, j, 2] = Taylor1(constant_term(tmp1079[i, j, 2]) - constant_term(tmp1080[i, j, 2]), order)
                        tmp1082[i, j, 2, 2] = Taylor1(constant_term(tmp1078[i, j, 2, 2]) * constant_term(tmp1081[i, j, 2]), order)
                        tmp1083[i, j, 2, 1] = Taylor1(constant_term(tmp1077[i, j, 2, 1]) + constant_term(tmp1082[i, j, 2, 2]), order)
                        F_CS_η[i, j] = Taylor1(constant_term(tmp1083[i, j, 2, 1]) / constant_term(r_p4[i, j]), order)
                        tmp1085[i, j, 1] = Taylor1(constant_term(C21M_t) * constant_term(cos_mλ[i, j, 1]), order)
                        tmp1086[i, j, 1] = Taylor1(constant_term(S21M_t) * constant_term(sin_mλ[i, j, 1]), order)
                        tmp1087[i, j, 1] = Taylor1(constant_term(tmp1085[i, j, 1]) + constant_term(tmp1086[i, j, 1]), order)
                        tmp1088[i, j, 2, 1] = Taylor1(constant_term(cosϕ_dP_nm[i, j, 2, 1]) * constant_term(tmp1087[i, j, 1]), order)
                        tmp1089[i, j, 2] = Taylor1(constant_term(C22M_t) * constant_term(cos_mλ[i, j, 2]), order)
                        tmp1090[i, j, 2] = Taylor1(constant_term(S22M_t) * constant_term(sin_mλ[i, j, 2]), order)
                        tmp1091[i, j, 2] = Taylor1(constant_term(tmp1089[i, j, 2]) + constant_term(tmp1090[i, j, 2]), order)
                        tmp1092[i, j, 2, 2] = Taylor1(constant_term(cosϕ_dP_nm[i, j, 2, 2]) * constant_term(tmp1091[i, j, 2]), order)
                        tmp1093[i, j, 2, 1] = Taylor1(constant_term(tmp1088[i, j, 2, 1]) + constant_term(tmp1092[i, j, 2, 2]), order)
                        F_CS_ζ[i, j] = Taylor1(constant_term(tmp1093[i, j, 2, 1]) / constant_term(r_p4[i, j]), order)
                        F_CS_ξ_36[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                        F_CS_η_36[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                        F_CS_ζ_36[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                        for n = 3:n1SEM[mo]
                            for m = 1:n
                                tmp1095[i, j, n, m] = Taylor1(constant_term(P_nm[i, j, n, m]) * constant_term(lnm6[n]), order)
                                tmp1096[i, j, m] = Taylor1(constant_term(cos_mλ[i, j, m]) * constant_term(CM[n, m]), order)
                                tmp1097[i, j, m] = Taylor1(constant_term(sin_mλ[i, j, m]) * constant_term(SM[n, m]), order)
                                tmp1098[i, j, m] = Taylor1(constant_term(tmp1096[i, j, m]) + constant_term(tmp1097[i, j, m]), order)
                                tmp1099[i, j, n, m] = Taylor1(constant_term(tmp1095[i, j, n, m]) * constant_term(tmp1098[i, j, m]), order)
                                tmp1100[i, j, n, m] = Taylor1(constant_term(tmp1099[i, j, n, m]) / constant_term(temp_rn[i, j, n]), order)
                                temp_CS_ξ[i, j, n, m] = Taylor1(constant_term(tmp1100[i, j, n, m]) + constant_term(F_CS_ξ_36[i, j]), order)
                                tmp1102[i, j, n, m] = Taylor1(constant_term(secϕ_P_nm[i, j, n, m]) * constant_term(lnm7[m]), order)
                                tmp1103[i, j, m] = Taylor1(constant_term(cos_mλ[i, j, m]) * constant_term(SM[n, m]), order)
                                tmp1104[i, j, m] = Taylor1(constant_term(sin_mλ[i, j, m]) * constant_term(CM[n, m]), order)
                                tmp1105[i, j, m] = Taylor1(constant_term(tmp1103[i, j, m]) - constant_term(tmp1104[i, j, m]), order)
                                tmp1106[i, j, n, m] = Taylor1(constant_term(tmp1102[i, j, n, m]) * constant_term(tmp1105[i, j, m]), order)
                                tmp1107[i, j, n, m] = Taylor1(constant_term(tmp1106[i, j, n, m]) / constant_term(temp_rn[i, j, n]), order)
                                temp_CS_η[i, j, n, m] = Taylor1(constant_term(tmp1107[i, j, n, m]) + constant_term(F_CS_η_36[i, j]), order)
                                tmp1109[i, j, m] = Taylor1(constant_term(cos_mλ[i, j, m]) * constant_term(CM[n, m]), order)
                                tmp1110[i, j, m] = Taylor1(constant_term(sin_mλ[i, j, m]) * constant_term(SM[n, m]), order)
                                tmp1111[i, j, m] = Taylor1(constant_term(tmp1109[i, j, m]) + constant_term(tmp1110[i, j, m]), order)
                                tmp1112[i, j, n, m] = Taylor1(constant_term(cosϕ_dP_nm[i, j, n, m]) * constant_term(tmp1111[i, j, m]), order)
                                tmp1113[i, j, n, m] = Taylor1(constant_term(tmp1112[i, j, n, m]) / constant_term(temp_rn[i, j, n]), order)
                                temp_CS_ζ[i, j, n, m] = Taylor1(constant_term(tmp1113[i, j, n, m]) + constant_term(F_CS_ζ_36[i, j]), order)
                                F_CS_ξ_36[i, j] = Taylor1(identity(constant_term(temp_CS_ξ[i, j, n, m])), order)
                                F_CS_η_36[i, j] = Taylor1(identity(constant_term(temp_CS_η[i, j, n, m])), order)
                                F_CS_ζ_36[i, j] = Taylor1(identity(constant_term(temp_CS_ζ[i, j, n, m])), order)
                            end
                        end
                        tmp1115[i, j] = Taylor1(constant_term(F_J_ξ[i, j]) + constant_term(F_J_ξ_36[i, j]), order)
                        tmp1116[i, j] = Taylor1(constant_term(F_CS_ξ[i, j]) + constant_term(F_CS_ξ_36[i, j]), order)
                        F_JCS_ξ[i, j] = Taylor1(constant_term(tmp1115[i, j]) + constant_term(tmp1116[i, j]), order)
                        F_JCS_η[i, j] = Taylor1(constant_term(F_CS_η[i, j]) + constant_term(F_CS_η_36[i, j]), order)
                        tmp1119[i, j] = Taylor1(constant_term(F_J_ζ[i, j]) + constant_term(F_J_ζ_36[i, j]), order)
                        tmp1120[i, j] = Taylor1(constant_term(F_CS_ζ[i, j]) + constant_term(F_CS_ζ_36[i, j]), order)
                        F_JCS_ζ[i, j] = Taylor1(constant_term(tmp1119[i, j]) + constant_term(tmp1120[i, j]), order)
                    else
                        F_JCS_ξ[i, j] = Taylor1(constant_term(F_J_ξ[i, j]) + constant_term(F_J_ξ_36[i, j]), order)
                        F_JCS_η[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                        F_JCS_ζ[i, j] = Taylor1(constant_term(F_J_ζ[i, j]) + constant_term(F_J_ζ_36[i, j]), order)
                    end
                    Rb2p[i, j, 1, 1] = Taylor1(constant_term(cos_ϕ[i, j]) * constant_term(cos_λ[i, j]), order)
                    Rb2p[i, j, 2, 1] = Taylor1(-(constant_term(sin_λ[i, j])), order)
                    tmp1126[i, j] = Taylor1(-(constant_term(sin_ϕ[i, j])), order)
                    Rb2p[i, j, 3, 1] = Taylor1(constant_term(tmp1126[i, j]) * constant_term(cos_λ[i, j]), order)
                    Rb2p[i, j, 1, 2] = Taylor1(constant_term(cos_ϕ[i, j]) * constant_term(sin_λ[i, j]), order)
                    Rb2p[i, j, 2, 2] = Taylor1(identity(constant_term(cos_λ[i, j])), order)
                    tmp1129[i, j] = Taylor1(-(constant_term(sin_ϕ[i, j])), order)
                    Rb2p[i, j, 3, 2] = Taylor1(constant_term(tmp1129[i, j]) * constant_term(sin_λ[i, j]), order)
                    Rb2p[i, j, 1, 3] = Taylor1(identity(constant_term(sin_ϕ[i, j])), order)
                    Rb2p[i, j, 2, 3] = Taylor1(identity(constant_term(zero_q_1)), order)
                    Rb2p[i, j, 3, 3] = Taylor1(identity(constant_term(cos_ϕ[i, j])), order)
                    tmp1131[i, j, 1, 1] = Taylor1(constant_term(Rb2p[i, j, 1, 1]) * constant_term(M_[1, 1, j]), order)
                    tmp1132[i, j, 1, 2] = Taylor1(constant_term(Rb2p[i, j, 1, 2]) * constant_term(M_[2, 1, j]), order)
                    tmp1133[i, j, 1, 1] = Taylor1(constant_term(tmp1131[i, j, 1, 1]) + constant_term(tmp1132[i, j, 1, 2]), order)
                    tmp1134[i, j, 1, 3] = Taylor1(constant_term(Rb2p[i, j, 1, 3]) * constant_term(M_[3, 1, j]), order)
                    Gc2p[i, j, 1, 1] = Taylor1(constant_term(tmp1133[i, j, 1, 1]) + constant_term(tmp1134[i, j, 1, 3]), order)
                    tmp1136[i, j, 2, 1] = Taylor1(constant_term(Rb2p[i, j, 2, 1]) * constant_term(M_[1, 1, j]), order)
                    tmp1137[i, j, 2, 2] = Taylor1(constant_term(Rb2p[i, j, 2, 2]) * constant_term(M_[2, 1, j]), order)
                    tmp1138[i, j, 2, 1] = Taylor1(constant_term(tmp1136[i, j, 2, 1]) + constant_term(tmp1137[i, j, 2, 2]), order)
                    tmp1139[i, j, 2, 3] = Taylor1(constant_term(Rb2p[i, j, 2, 3]) * constant_term(M_[3, 1, j]), order)
                    Gc2p[i, j, 2, 1] = Taylor1(constant_term(tmp1138[i, j, 2, 1]) + constant_term(tmp1139[i, j, 2, 3]), order)
                    tmp1141[i, j, 3, 1] = Taylor1(constant_term(Rb2p[i, j, 3, 1]) * constant_term(M_[1, 1, j]), order)
                    tmp1142[i, j, 3, 2] = Taylor1(constant_term(Rb2p[i, j, 3, 2]) * constant_term(M_[2, 1, j]), order)
                    tmp1143[i, j, 3, 1] = Taylor1(constant_term(tmp1141[i, j, 3, 1]) + constant_term(tmp1142[i, j, 3, 2]), order)
                    tmp1144[i, j, 3, 3] = Taylor1(constant_term(Rb2p[i, j, 3, 3]) * constant_term(M_[3, 1, j]), order)
                    Gc2p[i, j, 3, 1] = Taylor1(constant_term(tmp1143[i, j, 3, 1]) + constant_term(tmp1144[i, j, 3, 3]), order)
                    tmp1146[i, j, 1, 1] = Taylor1(constant_term(Rb2p[i, j, 1, 1]) * constant_term(M_[1, 2, j]), order)
                    tmp1147[i, j, 1, 2] = Taylor1(constant_term(Rb2p[i, j, 1, 2]) * constant_term(M_[2, 2, j]), order)
                    tmp1148[i, j, 1, 1] = Taylor1(constant_term(tmp1146[i, j, 1, 1]) + constant_term(tmp1147[i, j, 1, 2]), order)
                    tmp1149[i, j, 1, 3] = Taylor1(constant_term(Rb2p[i, j, 1, 3]) * constant_term(M_[3, 2, j]), order)
                    Gc2p[i, j, 1, 2] = Taylor1(constant_term(tmp1148[i, j, 1, 1]) + constant_term(tmp1149[i, j, 1, 3]), order)
                    tmp1151[i, j, 2, 1] = Taylor1(constant_term(Rb2p[i, j, 2, 1]) * constant_term(M_[1, 2, j]), order)
                    tmp1152[i, j, 2, 2] = Taylor1(constant_term(Rb2p[i, j, 2, 2]) * constant_term(M_[2, 2, j]), order)
                    tmp1153[i, j, 2, 1] = Taylor1(constant_term(tmp1151[i, j, 2, 1]) + constant_term(tmp1152[i, j, 2, 2]), order)
                    tmp1154[i, j, 2, 3] = Taylor1(constant_term(Rb2p[i, j, 2, 3]) * constant_term(M_[3, 2, j]), order)
                    Gc2p[i, j, 2, 2] = Taylor1(constant_term(tmp1153[i, j, 2, 1]) + constant_term(tmp1154[i, j, 2, 3]), order)
                    tmp1156[i, j, 3, 1] = Taylor1(constant_term(Rb2p[i, j, 3, 1]) * constant_term(M_[1, 2, j]), order)
                    tmp1157[i, j, 3, 2] = Taylor1(constant_term(Rb2p[i, j, 3, 2]) * constant_term(M_[2, 2, j]), order)
                    tmp1158[i, j, 3, 1] = Taylor1(constant_term(tmp1156[i, j, 3, 1]) + constant_term(tmp1157[i, j, 3, 2]), order)
                    tmp1159[i, j, 3, 3] = Taylor1(constant_term(Rb2p[i, j, 3, 3]) * constant_term(M_[3, 2, j]), order)
                    Gc2p[i, j, 3, 2] = Taylor1(constant_term(tmp1158[i, j, 3, 1]) + constant_term(tmp1159[i, j, 3, 3]), order)
                    tmp1161[i, j, 1, 1] = Taylor1(constant_term(Rb2p[i, j, 1, 1]) * constant_term(M_[1, 3, j]), order)
                    tmp1162[i, j, 1, 2] = Taylor1(constant_term(Rb2p[i, j, 1, 2]) * constant_term(M_[2, 3, j]), order)
                    tmp1163[i, j, 1, 1] = Taylor1(constant_term(tmp1161[i, j, 1, 1]) + constant_term(tmp1162[i, j, 1, 2]), order)
                    tmp1164[i, j, 1, 3] = Taylor1(constant_term(Rb2p[i, j, 1, 3]) * constant_term(M_[3, 3, j]), order)
                    Gc2p[i, j, 1, 3] = Taylor1(constant_term(tmp1163[i, j, 1, 1]) + constant_term(tmp1164[i, j, 1, 3]), order)
                    tmp1166[i, j, 2, 1] = Taylor1(constant_term(Rb2p[i, j, 2, 1]) * constant_term(M_[1, 3, j]), order)
                    tmp1167[i, j, 2, 2] = Taylor1(constant_term(Rb2p[i, j, 2, 2]) * constant_term(M_[2, 3, j]), order)
                    tmp1168[i, j, 2, 1] = Taylor1(constant_term(tmp1166[i, j, 2, 1]) + constant_term(tmp1167[i, j, 2, 2]), order)
                    tmp1169[i, j, 2, 3] = Taylor1(constant_term(Rb2p[i, j, 2, 3]) * constant_term(M_[3, 3, j]), order)
                    Gc2p[i, j, 2, 3] = Taylor1(constant_term(tmp1168[i, j, 2, 1]) + constant_term(tmp1169[i, j, 2, 3]), order)
                    tmp1171[i, j, 3, 1] = Taylor1(constant_term(Rb2p[i, j, 3, 1]) * constant_term(M_[1, 3, j]), order)
                    tmp1172[i, j, 3, 2] = Taylor1(constant_term(Rb2p[i, j, 3, 2]) * constant_term(M_[2, 3, j]), order)
                    tmp1173[i, j, 3, 1] = Taylor1(constant_term(tmp1171[i, j, 3, 1]) + constant_term(tmp1172[i, j, 3, 2]), order)
                    tmp1174[i, j, 3, 3] = Taylor1(constant_term(Rb2p[i, j, 3, 3]) * constant_term(M_[3, 3, j]), order)
                    Gc2p[i, j, 3, 3] = Taylor1(constant_term(tmp1173[i, j, 3, 1]) + constant_term(tmp1174[i, j, 3, 3]), order)
                    tmp1176[i, j, 1, 1] = Taylor1(constant_term(F_JCS_ξ[i, j]) * constant_term(Gc2p[i, j, 1, 1]), order)
                    tmp1177[i, j, 2, 1] = Taylor1(constant_term(F_JCS_η[i, j]) * constant_term(Gc2p[i, j, 2, 1]), order)
                    tmp1178[i, j, 1, 1] = Taylor1(constant_term(tmp1176[i, j, 1, 1]) + constant_term(tmp1177[i, j, 2, 1]), order)
                    tmp1179[i, j, 3, 1] = Taylor1(constant_term(F_JCS_ζ[i, j]) * constant_term(Gc2p[i, j, 3, 1]), order)
                    F_JCS_x[i, j] = Taylor1(constant_term(tmp1178[i, j, 1, 1]) + constant_term(tmp1179[i, j, 3, 1]), order)
                    tmp1181[i, j, 1, 2] = Taylor1(constant_term(F_JCS_ξ[i, j]) * constant_term(Gc2p[i, j, 1, 2]), order)
                    tmp1182[i, j, 2, 2] = Taylor1(constant_term(F_JCS_η[i, j]) * constant_term(Gc2p[i, j, 2, 2]), order)
                    tmp1183[i, j, 1, 2] = Taylor1(constant_term(tmp1181[i, j, 1, 2]) + constant_term(tmp1182[i, j, 2, 2]), order)
                    tmp1184[i, j, 3, 2] = Taylor1(constant_term(F_JCS_ζ[i, j]) * constant_term(Gc2p[i, j, 3, 2]), order)
                    F_JCS_y[i, j] = Taylor1(constant_term(tmp1183[i, j, 1, 2]) + constant_term(tmp1184[i, j, 3, 2]), order)
                    tmp1186[i, j, 1, 3] = Taylor1(constant_term(F_JCS_ξ[i, j]) * constant_term(Gc2p[i, j, 1, 3]), order)
                    tmp1187[i, j, 2, 3] = Taylor1(constant_term(F_JCS_η[i, j]) * constant_term(Gc2p[i, j, 2, 3]), order)
                    tmp1188[i, j, 1, 3] = Taylor1(constant_term(tmp1186[i, j, 1, 3]) + constant_term(tmp1187[i, j, 2, 3]), order)
                    tmp1189[i, j, 3, 3] = Taylor1(constant_term(F_JCS_ζ[i, j]) * constant_term(Gc2p[i, j, 3, 3]), order)
                    F_JCS_z[i, j] = Taylor1(constant_term(tmp1188[i, j, 1, 3]) + constant_term(tmp1189[i, j, 3, 3]), order)
                end
            end
        end
    end
    tmp1191 = Array{Taylor1{_S}}(undef, size(F_JCS_x))
    tmp1191 .= Taylor1(zero(_S), order)
    tmp1193 = Array{Taylor1{_S}}(undef, size(F_JCS_y))
    tmp1193 .= Taylor1(zero(_S), order)
    tmp1195 = Array{Taylor1{_S}}(undef, size(F_JCS_z))
    tmp1195 .= Taylor1(zero(_S), order)
    tmp1197 = Array{Taylor1{_S}}(undef, size(F_JCS_x))
    tmp1197 .= Taylor1(zero(_S), order)
    tmp1199 = Array{Taylor1{_S}}(undef, size(F_JCS_y))
    tmp1199 .= Taylor1(zero(_S), order)
    tmp1201 = Array{Taylor1{_S}}(undef, size(F_JCS_z))
    tmp1201 .= Taylor1(zero(_S), order)
    for j = 1:N
        for i = 1:N
            if i == j
                continue
            else
                if UJ_interaction[i, j]
                    tmp1191[i, j] = Taylor1(constant_term(μ[i]) * constant_term(F_JCS_x[i, j]), order)
                    temp_accX_j[i, j] = Taylor1(constant_term(accX[j]) - constant_term(tmp1191[i, j]), order)
                    accX[j] = Taylor1(identity(constant_term(temp_accX_j[i, j])), order)
                    tmp1193[i, j] = Taylor1(constant_term(μ[i]) * constant_term(F_JCS_y[i, j]), order)
                    temp_accY_j[i, j] = Taylor1(constant_term(accY[j]) - constant_term(tmp1193[i, j]), order)
                    accY[j] = Taylor1(identity(constant_term(temp_accY_j[i, j])), order)
                    tmp1195[i, j] = Taylor1(constant_term(μ[i]) * constant_term(F_JCS_z[i, j]), order)
                    temp_accZ_j[i, j] = Taylor1(constant_term(accZ[j]) - constant_term(tmp1195[i, j]), order)
                    accZ[j] = Taylor1(identity(constant_term(temp_accZ_j[i, j])), order)
                    tmp1197[i, j] = Taylor1(constant_term(μ[j]) * constant_term(F_JCS_x[i, j]), order)
                    temp_accX_i[i, j] = Taylor1(constant_term(accX[i]) + constant_term(tmp1197[i, j]), order)
                    accX[i] = Taylor1(identity(constant_term(temp_accX_i[i, j])), order)
                    tmp1199[i, j] = Taylor1(constant_term(μ[j]) * constant_term(F_JCS_y[i, j]), order)
                    temp_accY_i[i, j] = Taylor1(constant_term(accY[i]) + constant_term(tmp1199[i, j]), order)
                    accY[i] = Taylor1(identity(constant_term(temp_accY_i[i, j])), order)
                    tmp1201[i, j] = Taylor1(constant_term(μ[j]) * constant_term(F_JCS_z[i, j]), order)
                    temp_accZ_i[i, j] = Taylor1(constant_term(accZ[i]) + constant_term(tmp1201[i, j]), order)
                    accZ[i] = Taylor1(identity(constant_term(temp_accZ_i[i, j])), order)
                end
            end
        end
    end
    tmp1207 = Array{Taylor1{_S}}(undef, size(v2))
    tmp1207 .= Taylor1(zero(_S), order)
    tmp1208 = Array{Taylor1{_S}}(undef, size(v2))
    tmp1208 .= Taylor1(zero(_S), order)
    tmp1210 = Array{Taylor1{_S}}(undef, size(vi_dot_vj))
    tmp1210 .= Taylor1(zero(_S), order)
    Xij_t_Ui = Array{Taylor1{_S}}(undef, size(X))
    Xij_t_Ui .= Taylor1(zero(_S), order)
    Yij_t_Vi = Array{Taylor1{_S}}(undef, size(Y))
    Yij_t_Vi .= Taylor1(zero(_S), order)
    Zij_t_Wi = Array{Taylor1{_S}}(undef, size(Z))
    Zij_t_Wi .= Taylor1(zero(_S), order)
    tmp1216 = Array{Taylor1{_S}}(undef, size(Xij_t_Ui))
    tmp1216 .= Taylor1(zero(_S), order)
    Rij_dot_Vi = Array{Taylor1{_S}}(undef, size(tmp1216))
    Rij_dot_Vi .= Taylor1(zero(_S), order)
    tmp1219 = Array{Taylor1{_S}}(undef, size(Rij_dot_Vi))
    tmp1219 .= Taylor1(zero(_S), order)
    pn1t7 = Array{Taylor1{_S}}(undef, size(tmp1219))
    pn1t7 .= Taylor1(zero(_S), order)
    tmp1222 = Array{Taylor1{_S}}(undef, size(pn1t7))
    tmp1222 .= Taylor1(zero(_S), order)
    pn1t2_7 = Array{Taylor1{_S}}(undef, size(ϕs_and_vs))
    pn1t2_7 .= Taylor1(zero(_S), order)
    for j = 1:N
        for i = 1:N
            if i == j
                continue
            else
                _4ϕj[i, j] = Taylor1(constant_term(4) * constant_term(newtonianNb_Potential[j]), order)
                ϕi_plus_4ϕj[i, j] = Taylor1(constant_term(newtonianNb_Potential[i]) + constant_term(_4ϕj[i, j]), order)
                tmp1207[i] = Taylor1(constant_term(2) * constant_term(v2[i]), order)
                tmp1208[j] = Taylor1(constant_term(v2[j]) + constant_term(tmp1207[i]), order)
                tmp1210[i, j] = Taylor1(constant_term(4) * constant_term(vi_dot_vj[i, j]), order)
                sj2_plus_2si2_minus_4vivj[i, j] = Taylor1(constant_term(tmp1208[j]) - constant_term(tmp1210[i, j]), order)
                ϕs_and_vs[i, j] = Taylor1(constant_term(sj2_plus_2si2_minus_4vivj[i, j]) - constant_term(ϕi_plus_4ϕj[i, j]), order)
                Xij_t_Ui[i, j] = Taylor1(constant_term(X[i, j]) * constant_term(dq[3i - 2]), order)
                Yij_t_Vi[i, j] = Taylor1(constant_term(Y[i, j]) * constant_term(dq[3i - 1]), order)
                Zij_t_Wi[i, j] = Taylor1(constant_term(Z[i, j]) * constant_term(dq[3i]), order)
                tmp1216[i, j] = Taylor1(constant_term(Xij_t_Ui[i, j]) + constant_term(Yij_t_Vi[i, j]), order)
                Rij_dot_Vi[i, j] = Taylor1(constant_term(tmp1216[i, j]) + constant_term(Zij_t_Wi[i, j]), order)
                tmp1219[i, j] = Taylor1(constant_term(Rij_dot_Vi[i, j]) ^ constant_term(2), order)
                pn1t7[i, j] = Taylor1(constant_term(tmp1219[i, j]) / constant_term(r_p2[i, j]), order)
                tmp1222[i, j] = Taylor1(constant_term(1.5) * constant_term(pn1t7[i, j]), order)
                pn1t2_7[i, j] = Taylor1(constant_term(ϕs_and_vs[i, j]) - constant_term(tmp1222[i, j]), order)
                pn1t1_7[i, j] = Taylor1(constant_term(c_p2) + constant_term(pn1t2_7[i, j]), order)
                for k = 1:postnewton_iter
                    pn1[i, j, k] = Taylor1(identity(constant_term(zero_q_1)), order)
                    X_t_pn1[i, j, k] = Taylor1(identity(constant_term(zero_q_1)), order)
                    Y_t_pn1[i, j, k] = Taylor1(identity(constant_term(zero_q_1)), order)
                    Z_t_pn1[i, j, k] = Taylor1(identity(constant_term(zero_q_1)), order)
                    pNX_t_pn3[i, j, k] = Taylor1(identity(constant_term(zero_q_1)), order)
                    pNY_t_pn3[i, j, k] = Taylor1(identity(constant_term(zero_q_1)), order)
                    pNZ_t_pn3[i, j, k] = Taylor1(identity(constant_term(zero_q_1)), order)
                    pNX_t_X[i, j, k] = Taylor1(identity(constant_term(zero_q_1)), order)
                    pNY_t_Y[i, j, k] = Taylor1(identity(constant_term(zero_q_1)), order)
                    pNZ_t_Z[i, j, k] = Taylor1(identity(constant_term(zero_q_1)), order)
                end
            end
        end
        postNewtonX[j, 1] = Taylor1(identity(constant_term(newtonX[j])), order)
        postNewtonY[j, 1] = Taylor1(identity(constant_term(newtonY[j])), order)
        postNewtonZ[j, 1] = Taylor1(identity(constant_term(newtonZ[j])), order)
        for k = 1:postnewton_iter
            pntempX[j, k] = Taylor1(identity(constant_term(zero_q_1)), order)
            pntempY[j, k] = Taylor1(identity(constant_term(zero_q_1)), order)
            pntempZ[j, k] = Taylor1(identity(constant_term(zero_q_1)), order)
        end
    end
    tmp1229 = Array{Taylor1{_S}}(undef, size(pNX_t_X))
    tmp1229 .= Taylor1(zero(_S), order)
    tmp1230 = Array{Taylor1{_S}}(undef, size(tmp1229))
    tmp1230 .= Taylor1(zero(_S), order)
    tmp1231 = Array{Taylor1{_S}}(undef, size(tmp1230))
    tmp1231 .= Taylor1(zero(_S), order)
    tmp1239 = Array{Taylor1{_S}}(undef, size(pNX_t_pn3))
    tmp1239 .= Taylor1(zero(_S), order)
    termpnx = Array{Taylor1{_S}}(undef, size(X_t_pn1))
    termpnx .= Taylor1(zero(_S), order)
    sumpnx = Array{Taylor1{_S}}(undef, size(termpnx))
    sumpnx .= Taylor1(zero(_S), order)
    tmp1242 = Array{Taylor1{_S}}(undef, size(pNY_t_pn3))
    tmp1242 .= Taylor1(zero(_S), order)
    termpny = Array{Taylor1{_S}}(undef, size(Y_t_pn1))
    termpny .= Taylor1(zero(_S), order)
    sumpny = Array{Taylor1{_S}}(undef, size(termpny))
    sumpny .= Taylor1(zero(_S), order)
    tmp1245 = Array{Taylor1{_S}}(undef, size(pNZ_t_pn3))
    tmp1245 .= Taylor1(zero(_S), order)
    termpnz = Array{Taylor1{_S}}(undef, size(Z_t_pn1))
    termpnz .= Taylor1(zero(_S), order)
    sumpnz = Array{Taylor1{_S}}(undef, size(termpnz))
    sumpnz .= Taylor1(zero(_S), order)
    for k = 1:postnewton_iter
        for j = 1:N
            for i = 1:N
                if i == j
                    continue
                else
                    pNX_t_X[i, j, k] = Taylor1(constant_term(postNewtonX[i, k]) * constant_term(X[i, j]), order)
                    pNY_t_Y[i, j, k] = Taylor1(constant_term(postNewtonY[i, k]) * constant_term(Y[i, j]), order)
                    pNZ_t_Z[i, j, k] = Taylor1(constant_term(postNewtonZ[i, k]) * constant_term(Z[i, j]), order)
                    tmp1229[i, j, k] = Taylor1(constant_term(pNX_t_X[i, j, k]) + constant_term(pNY_t_Y[i, j, k]), order)
                    tmp1230[i, j, k] = Taylor1(constant_term(tmp1229[i, j, k]) + constant_term(pNZ_t_Z[i, j, k]), order)
                    tmp1231[i, j, k] = Taylor1(constant_term(0.5) * constant_term(tmp1230[i, j, k]), order)
                    pn1[i, j, k] = Taylor1(constant_term(pn1t1_7[i, j]) + constant_term(tmp1231[i, j, k]), order)
                    X_t_pn1[i, j, k] = Taylor1(constant_term(newton_acc_X[i, j]) * constant_term(pn1[i, j, k]), order)
                    Y_t_pn1[i, j, k] = Taylor1(constant_term(newton_acc_Y[i, j]) * constant_term(pn1[i, j, k]), order)
                    Z_t_pn1[i, j, k] = Taylor1(constant_term(newton_acc_Z[i, j]) * constant_term(pn1[i, j, k]), order)
                    pNX_t_pn3[i, j, k] = Taylor1(constant_term(postNewtonX[i, k]) * constant_term(pn3[i, j]), order)
                    pNY_t_pn3[i, j, k] = Taylor1(constant_term(postNewtonY[i, k]) * constant_term(pn3[i, j]), order)
                    pNZ_t_pn3[i, j, k] = Taylor1(constant_term(postNewtonZ[i, k]) * constant_term(pn3[i, j]), order)
                    tmp1239[i, j, k] = Taylor1(constant_term(U_t_pn2[i, j]) + constant_term(pNX_t_pn3[i, j, k]), order)
                    termpnx[i, j, k] = Taylor1(constant_term(X_t_pn1[i, j, k]) + constant_term(tmp1239[i, j, k]), order)
                    sumpnx[i, j, k] = Taylor1(constant_term(pntempX[j, k]) + constant_term(termpnx[i, j, k]), order)
                    pntempX[j, k] = Taylor1(identity(constant_term(sumpnx[i, j, k])), order)
                    tmp1242[i, j, k] = Taylor1(constant_term(V_t_pn2[i, j]) + constant_term(pNY_t_pn3[i, j, k]), order)
                    termpny[i, j, k] = Taylor1(constant_term(Y_t_pn1[i, j, k]) + constant_term(tmp1242[i, j, k]), order)
                    sumpny[i, j, k] = Taylor1(constant_term(pntempY[j, k]) + constant_term(termpny[i, j, k]), order)
                    pntempY[j, k] = Taylor1(identity(constant_term(sumpny[i, j, k])), order)
                    tmp1245[i, j, k] = Taylor1(constant_term(W_t_pn2[i, j]) + constant_term(pNZ_t_pn3[i, j, k]), order)
                    termpnz[i, j, k] = Taylor1(constant_term(Z_t_pn1[i, j, k]) + constant_term(tmp1245[i, j, k]), order)
                    sumpnz[i, j, k] = Taylor1(constant_term(pntempZ[j, k]) + constant_term(termpnz[i, j, k]), order)
                    pntempZ[j, k] = Taylor1(identity(constant_term(sumpnz[i, j, k])), order)
                end
            end
            postNewtonX[j, k + 1] = Taylor1(constant_term(pntempX[j, k]) * constant_term(c_m2), order)
            postNewtonY[j, k + 1] = Taylor1(constant_term(pntempY[j, k]) * constant_term(c_m2), order)
            postNewtonZ[j, k + 1] = Taylor1(constant_term(pntempZ[j, k]) * constant_term(c_m2), order)
        end
    end
    for i = 1:N
        dq[3 * (N + i) - 2] = Taylor1(constant_term(postNewtonX[i, postnewton_iter + 1]) + constant_term(accX[i]), order)
        dq[3 * (N + i) - 1] = Taylor1(constant_term(postNewtonY[i, postnewton_iter + 1]) + constant_term(accY[i]), order)
        dq[3 * (N + i)] = Taylor1(constant_term(postNewtonZ[i, postnewton_iter + 1]) + constant_term(accZ[i]), order)
    end
    for __idx = eachindex(q)
        (q[__idx]).coeffs[2] = (dq[__idx]).coeffs[1]
    end
    for ord = 1:order - 1
        ordnext = ord + 1
        for j = 1:N
            TaylorSeries.identity!(newtonX[j], zero_q_1, ord)
            TaylorSeries.identity!(newtonY[j], zero_q_1, ord)
            TaylorSeries.identity!(newtonZ[j], zero_q_1, ord)
            TaylorSeries.identity!(newtonianNb_Potential[j], zero_q_1, ord)
            TaylorSeries.identity!(accX[j], zero_q_1, ord)
            TaylorSeries.identity!(accY[j], zero_q_1, ord)
            TaylorSeries.identity!(accZ[j], zero_q_1, ord)
            TaylorSeries.identity!(dq[3j - 2], q[3 * (N + j) - 2], ord)
            TaylorSeries.identity!(dq[3j - 1], q[3 * (N + j) - 1], ord)
            TaylorSeries.identity!(dq[3j], q[3 * (N + j)], ord)
        end
        for j = 1:N
            for i = 1:N
                if i == j
                    continue
                else
                    TaylorSeries.subst!(X[i, j], q[3i - 2], q[3j - 2], ord)
                    TaylorSeries.subst!(Y[i, j], q[3i - 1], q[3j - 1], ord)
                    TaylorSeries.subst!(Z[i, j], q[3i], q[3j], ord)
                    TaylorSeries.subst!(U[i, j], dq[3i - 2], dq[3j - 2], ord)
                    TaylorSeries.subst!(V[i, j], dq[3i - 1], dq[3j - 1], ord)
                    TaylorSeries.subst!(W[i, j], dq[3i], dq[3j], ord)
                    TaylorSeries.mul!(tmp894[3j - 2], 4, dq[3j - 2], ord)
                    TaylorSeries.mul!(tmp896[3i - 2], 3, dq[3i - 2], ord)
                    TaylorSeries.subst!(_4U_m_3X[i, j], tmp894[3j - 2], tmp896[3i - 2], ord)
                    TaylorSeries.mul!(tmp899[3j - 1], 4, dq[3j - 1], ord)
                    TaylorSeries.mul!(tmp901[3i - 1], 3, dq[3i - 1], ord)
                    TaylorSeries.subst!(_4V_m_3Y[i, j], tmp899[3j - 1], tmp901[3i - 1], ord)
                    TaylorSeries.mul!(tmp904[3j], 4, dq[3j], ord)
                    TaylorSeries.mul!(tmp906[3i], 3, dq[3i], ord)
                    TaylorSeries.subst!(_4W_m_3Z[i, j], tmp904[3j], tmp906[3i], ord)
                    TaylorSeries.mul!(pn2x[i, j], X[i, j], _4U_m_3X[i, j], ord)
                    TaylorSeries.mul!(pn2y[i, j], Y[i, j], _4V_m_3Y[i, j], ord)
                    TaylorSeries.mul!(pn2z[i, j], Z[i, j], _4W_m_3Z[i, j], ord)
                    TaylorSeries.mul!(UU[i, j], dq[3i - 2], dq[3j - 2], ord)
                    TaylorSeries.mul!(VV[i, j], dq[3i - 1], dq[3j - 1], ord)
                    TaylorSeries.mul!(WW[i, j], dq[3i], dq[3j], ord)
                    TaylorSeries.add!(tmp914[i, j], UU[i, j], VV[i, j], ord)
                    TaylorSeries.add!(vi_dot_vj[i, j], tmp914[i, j], WW[i, j], ord)
                    TaylorSeries.pow!(tmp917[i, j], X[i, j], 2, ord)
                    TaylorSeries.pow!(tmp919[i, j], Y[i, j], 2, ord)
                    TaylorSeries.add!(tmp920[i, j], tmp917[i, j], tmp919[i, j], ord)
                    TaylorSeries.pow!(tmp922[i, j], Z[i, j], 2, ord)
                    TaylorSeries.add!(r_p2[i, j], tmp920[i, j], tmp922[i, j], ord)
                    TaylorSeries.sqrt!(r_p1d2[i, j], r_p2[i, j], ord)
                    TaylorSeries.pow!(r_p3d2[i, j], r_p2[i, j], 1.5, ord)
                    TaylorSeries.pow!(r_p7d2[i, j], r_p2[i, j], 3.5, ord)
                    TaylorSeries.div!(newtonianCoeff[i, j], μ[i], r_p3d2[i, j], ord)
                    TaylorSeries.add!(tmp930[i, j], pn2x[i, j], pn2y[i, j], ord)
                    TaylorSeries.add!(tmp931[i, j], tmp930[i, j], pn2z[i, j], ord)
                    TaylorSeries.mul!(pn2[i, j], newtonianCoeff[i, j], tmp931[i, j], ord)
                    TaylorSeries.mul!(newton_acc_X[i, j], X[i, j], newtonianCoeff[i, j], ord)
                    TaylorSeries.mul!(newton_acc_Y[i, j], Y[i, j], newtonianCoeff[i, j], ord)
                    TaylorSeries.mul!(newton_acc_Z[i, j], Z[i, j], newtonianCoeff[i, j], ord)
                    TaylorSeries.div!(newtonian1b_Potential[i, j], μ[i], r_p1d2[i, j], ord)
                    TaylorSeries.mul!(pn3[i, j], 3.5, newtonian1b_Potential[i, j], ord)
                    TaylorSeries.mul!(U_t_pn2[i, j], pn2[i, j], U[i, j], ord)
                    TaylorSeries.mul!(V_t_pn2[i, j], pn2[i, j], V[i, j], ord)
                    TaylorSeries.mul!(W_t_pn2[i, j], pn2[i, j], W[i, j], ord)
                    TaylorSeries.mul!(tmp942[i, j], X[i, j], newtonianCoeff[i, j], ord)
                    TaylorSeries.add!(temp_001[i, j], newtonX[j], tmp942[i, j], ord)
                    TaylorSeries.identity!(newtonX[j], temp_001[i, j], ord)
                    TaylorSeries.mul!(tmp944[i, j], Y[i, j], newtonianCoeff[i, j], ord)
                    TaylorSeries.add!(temp_002[i, j], newtonY[j], tmp944[i, j], ord)
                    TaylorSeries.identity!(newtonY[j], temp_002[i, j], ord)
                    TaylorSeries.mul!(tmp946[i, j], Z[i, j], newtonianCoeff[i, j], ord)
                    TaylorSeries.add!(temp_003[i, j], newtonZ[j], tmp946[i, j], ord)
                    TaylorSeries.identity!(newtonZ[j], temp_003[i, j], ord)
                    TaylorSeries.add!(temp_004[i, j], newtonianNb_Potential[j], newtonian1b_Potential[i, j], ord)
                    TaylorSeries.identity!(newtonianNb_Potential[j], temp_004[i, j], ord)
                end
            end
            TaylorSeries.pow!(tmp950[3j - 2], dq[3j - 2], 2, ord)
            TaylorSeries.pow!(tmp952[3j - 1], dq[3j - 1], 2, ord)
            TaylorSeries.add!(tmp953[3j - 2], tmp950[3j - 2], tmp952[3j - 1], ord)
            TaylorSeries.pow!(tmp955[3j], dq[3j], 2, ord)
            TaylorSeries.add!(v2[j], tmp953[3j - 2], tmp955[3j], ord)
        end
        TaylorSeries.add!(tmp957, ITM_t[1, 1], ITM_t[2, 2], ord)
        TaylorSeries.div!(tmp959, tmp957, 2, ord)
        TaylorSeries.subst!(tmp960, ITM_t[3, 3], tmp959, ord)
        TaylorSeries.div!(J2M_t, tmp960, μ[mo], ord)
        TaylorSeries.subst!(tmp962, ITM_t[2, 2], ITM_t[1, 1], ord)
        TaylorSeries.div!(tmp963, tmp962, μ[mo], ord)
        TaylorSeries.div!(C22M_t, tmp963, 4, ord)
        TaylorSeries.subst!(tmp966, ITM_t[1, 3], ord)
        TaylorSeries.div!(C21M_t, tmp966, μ[mo], ord)
        TaylorSeries.subst!(tmp968, ITM_t[3, 2], ord)
        TaylorSeries.div!(S21M_t, tmp968, μ[mo], ord)
        TaylorSeries.subst!(tmp970, ITM_t[2, 1], ord)
        TaylorSeries.div!(tmp971, tmp970, μ[mo], ord)
        TaylorSeries.div!(S22M_t, tmp971, 2, ord)
        for j = 1:N
            for i = 1:N
                if i == j
                    continue
                else
                    if UJ_interaction[i, j]
                        TaylorSeries.mul!(X_bf_1[i, j], X[i, j], M_[1, 1, j], ord)
                        TaylorSeries.mul!(X_bf_2[i, j], Y[i, j], M_[1, 2, j], ord)
                        TaylorSeries.mul!(X_bf_3[i, j], Z[i, j], M_[1, 3, j], ord)
                        TaylorSeries.mul!(Y_bf_1[i, j], X[i, j], M_[2, 1, j], ord)
                        TaylorSeries.mul!(Y_bf_2[i, j], Y[i, j], M_[2, 2, j], ord)
                        TaylorSeries.mul!(Y_bf_3[i, j], Z[i, j], M_[2, 3, j], ord)
                        TaylorSeries.mul!(Z_bf_1[i, j], X[i, j], M_[3, 1, j], ord)
                        TaylorSeries.mul!(Z_bf_2[i, j], Y[i, j], M_[3, 2, j], ord)
                        TaylorSeries.mul!(Z_bf_3[i, j], Z[i, j], M_[3, 3, j], ord)
                        TaylorSeries.add!(tmp983[i, j], X_bf_1[i, j], X_bf_2[i, j], ord)
                        TaylorSeries.add!(X_bf[i, j], tmp983[i, j], X_bf_3[i, j], ord)
                        TaylorSeries.add!(tmp985[i, j], Y_bf_1[i, j], Y_bf_2[i, j], ord)
                        TaylorSeries.add!(Y_bf[i, j], tmp985[i, j], Y_bf_3[i, j], ord)
                        TaylorSeries.add!(tmp987[i, j], Z_bf_1[i, j], Z_bf_2[i, j], ord)
                        TaylorSeries.add!(Z_bf[i, j], tmp987[i, j], Z_bf_3[i, j], ord)
                        TaylorSeries.div!(sin_ϕ[i, j], Z_bf[i, j], r_p1d2[i, j], ord)
                        TaylorSeries.pow!(tmp991[i, j], X_bf[i, j], 2, ord)
                        TaylorSeries.pow!(tmp993[i, j], Y_bf[i, j], 2, ord)
                        TaylorSeries.add!(tmp994[i, j], tmp991[i, j], tmp993[i, j], ord)
                        TaylorSeries.sqrt!(r_xy[i, j], tmp994[i, j], ord)
                        TaylorSeries.div!(cos_ϕ[i, j], r_xy[i, j], r_p1d2[i, j], ord)
                        TaylorSeries.div!(sin_λ[i, j], Y_bf[i, j], r_xy[i, j], ord)
                        TaylorSeries.div!(cos_λ[i, j], X_bf[i, j], r_xy[i, j], ord)
                        TaylorSeries.identity!(P_n[i, j, 1], one_t, ord)
                        TaylorSeries.identity!(P_n[i, j, 2], sin_ϕ[i, j], ord)
                        TaylorSeries.identity!(dP_n[i, j, 1], zero_q_1, ord)
                        TaylorSeries.identity!(dP_n[i, j, 2], one_t, ord)
                        for n = 2:n1SEM[j]
                            TaylorSeries.mul!(tmp999[i, j, n], P_n[i, j, n], sin_ϕ[i, j], ord)
                            TaylorSeries.mul!(tmp1000[i, j, n], tmp999[i, j, n], fact1_jsem[n], ord)
                            TaylorSeries.mul!(tmp1001[i, j, n - 1], P_n[i, j, n - 1], fact2_jsem[n], ord)
                            TaylorSeries.subst!(P_n[i, j, n + 1], tmp1000[i, j, n], tmp1001[i, j, n - 1], ord)
                            TaylorSeries.mul!(tmp1003[i, j, n], dP_n[i, j, n], sin_ϕ[i, j], ord)
                            TaylorSeries.mul!(tmp1004[i, j, n], P_n[i, j, n], fact3_jsem[n], ord)
                            TaylorSeries.add!(dP_n[i, j, n + 1], tmp1003[i, j, n], tmp1004[i, j, n], ord)
                            TaylorSeries.pow!(temp_rn[i, j, n], r_p1d2[i, j], fact5_jsem[n], ord)
                        end
                        TaylorSeries.pow!(r_p4[i, j], r_p2[i, j], 2, ord)
                        if j == mo
                            TaylorSeries.mul!(tmp1009[i, j, 3], P_n[i, j, 3], fact4_jsem[2], ord)
                            TaylorSeries.mul!(tmp1010[i, j, 3], tmp1009[i, j, 3], J2M_t, ord)
                            TaylorSeries.div!(F_J_ξ[i, j], tmp1010[i, j, 3], r_p4[i, j], ord)
                            TaylorSeries.subst!(tmp1012[i, j, 3], dP_n[i, j, 3], ord)
                            TaylorSeries.mul!(tmp1013[i, j, 3], tmp1012[i, j, 3], cos_ϕ[i, j], ord)
                            TaylorSeries.mul!(tmp1014[i, j, 3], tmp1013[i, j, 3], J2M_t, ord)
                            TaylorSeries.div!(F_J_ζ[i, j], tmp1014[i, j, 3], r_p4[i, j], ord)
                        else
                            if j == ea
                                TaylorSeries.mul!(tmp1016[i, j, 3], P_n[i, j, 3], fact4_jsem[2], ord)
                                TaylorSeries.mul!(tmp1017[i, j, 3], tmp1016[i, j, 3], J2E_t, ord)
                                TaylorSeries.div!(F_J_ξ[i, j], tmp1017[i, j, 3], r_p4[i, j], ord)
                                TaylorSeries.subst!(tmp1019[i, j, 3], dP_n[i, j, 3], ord)
                                TaylorSeries.mul!(tmp1020[i, j, 3], tmp1019[i, j, 3], cos_ϕ[i, j], ord)
                                TaylorSeries.mul!(tmp1021[i, j, 3], tmp1020[i, j, 3], J2E_t, ord)
                                TaylorSeries.div!(F_J_ζ[i, j], tmp1021[i, j, 3], r_p4[i, j], ord)
                            else
                                if j == su
                                    TaylorSeries.mul!(tmp1023[i, j, 3], P_n[i, j, 3], fact4_jsem[2], ord)
                                    TaylorSeries.mul!(tmp1024[i, j, 3], tmp1023[i, j, 3], J2S_t, ord)
                                    TaylorSeries.div!(F_J_ξ[i, j], tmp1024[i, j, 3], r_p4[i, j], ord)
                                    TaylorSeries.subst!(tmp1026[i, j, 3], dP_n[i, j, 3], ord)
                                    TaylorSeries.mul!(tmp1027[i, j, 3], tmp1026[i, j, 3], cos_ϕ[i, j], ord)
                                    TaylorSeries.mul!(tmp1028[i, j, 3], tmp1027[i, j, 3], J2S_t, ord)
                                    TaylorSeries.div!(F_J_ζ[i, j], tmp1028[i, j, 3], r_p4[i, j], ord)
                                end
                            end
                        end
                        TaylorSeries.identity!(F_J_ξ_36[i, j], zero_q_1, ord)
                        TaylorSeries.identity!(F_J_ζ_36[i, j], zero_q_1, ord)
                        for n = 3:n1SEM[j]
                            TaylorSeries.mul!(tmp1030[i, j, n + 1], P_n[i, j, n + 1], fact4_jsem[n], ord)
                            TaylorSeries.mul!(tmp1031[i, j, n + 1], tmp1030[i, j, n + 1], JSEM[j, n], ord)
                            TaylorSeries.div!(tmp1032[i, j, n + 1], tmp1031[i, j, n + 1], temp_rn[i, j, n], ord)
                            TaylorSeries.add!(temp_fjξ[i, j, n], tmp1032[i, j, n + 1], F_J_ξ_36[i, j], ord)
                            TaylorSeries.subst!(tmp1034[i, j, n + 1], dP_n[i, j, n + 1], ord)
                            TaylorSeries.mul!(tmp1035[i, j, n + 1], tmp1034[i, j, n + 1], cos_ϕ[i, j], ord)
                            TaylorSeries.mul!(tmp1036[i, j, n + 1], tmp1035[i, j, n + 1], JSEM[j, n], ord)
                            TaylorSeries.div!(tmp1037[i, j, n + 1], tmp1036[i, j, n + 1], temp_rn[i, j, n], ord)
                            TaylorSeries.add!(temp_fjζ[i, j, n], tmp1037[i, j, n + 1], F_J_ζ_36[i, j], ord)
                            TaylorSeries.identity!(F_J_ξ_36[i, j], temp_fjξ[i, j, n], ord)
                            TaylorSeries.identity!(F_J_ζ_36[i, j], temp_fjζ[i, j, n], ord)
                        end
                        if j == mo
                            for m = 1:n1SEM[mo]
                                if m == 1
                                    TaylorSeries.identity!(sin_mλ[i, j, 1], sin_λ[i, j], ord)
                                    TaylorSeries.identity!(cos_mλ[i, j, 1], cos_λ[i, j], ord)
                                    TaylorSeries.identity!(secϕ_P_nm[i, j, 1, 1], one_t, ord)
                                else
                                    TaylorSeries.mul!(tmp1039[i, j, 1], sin_mλ[i, j, 1], cos_mλ[i, j, m - 1], ord)
                                    TaylorSeries.mul!(tmp1040[i, j, 1], cos_mλ[i, j, 1], sin_mλ[i, j, m - 1], ord)
                                    TaylorSeries.add!(sin_mλ[i, j, m], tmp1039[i, j, 1], tmp1040[i, j, 1], ord)
                                    TaylorSeries.mul!(tmp1042[i, j, 1], cos_mλ[i, j, 1], cos_mλ[i, j, m - 1], ord)
                                    TaylorSeries.mul!(tmp1043[i, j, 1], sin_mλ[i, j, 1], sin_mλ[i, j, m - 1], ord)
                                    TaylorSeries.subst!(cos_mλ[i, j, m], tmp1042[i, j, 1], tmp1043[i, j, 1], ord)
                                    TaylorSeries.mul!(tmp1045[i, j, m - 1, m - 1], secϕ_P_nm[i, j, m - 1, m - 1], cos_ϕ[i, j], ord)
                                    TaylorSeries.mul!(secϕ_P_nm[i, j, m, m], tmp1045[i, j, m - 1, m - 1], lnm5[m], ord)
                                    TaylorSeries.mul!(P_nm[i, j, m, m], secϕ_P_nm[i, j, m, m], cos_ϕ[i, j], ord)
                                    TaylorSeries.mul!(tmp1048[i, j, m, m], secϕ_P_nm[i, j, m, m], sin_ϕ[i, j], ord)
                                    TaylorSeries.mul!(cosϕ_dP_nm[i, j, m, m], tmp1048[i, j, m, m], lnm3[m], ord)
                                end
                                for n = m + 1:n1SEM[mo]
                                    if n == m + 1
                                        TaylorSeries.mul!(tmp1050[i, j, n - 1, m], secϕ_P_nm[i, j, n - 1, m], sin_ϕ[i, j], ord)
                                        TaylorSeries.mul!(secϕ_P_nm[i, j, n, m], tmp1050[i, j, n - 1, m], lnm1[n, m], ord)
                                    else
                                        TaylorSeries.mul!(tmp1052[i, j, n - 1, m], secϕ_P_nm[i, j, n - 1, m], sin_ϕ[i, j], ord)
                                        TaylorSeries.mul!(tmp1053[i, j, n - 1, m], tmp1052[i, j, n - 1, m], lnm1[n, m], ord)
                                        TaylorSeries.mul!(tmp1054[i, j, n - 2, m], secϕ_P_nm[i, j, n - 2, m], lnm2[n, m], ord)
                                        TaylorSeries.add!(secϕ_P_nm[i, j, n, m], tmp1053[i, j, n - 1, m], tmp1054[i, j, n - 2, m], ord)
                                    end
                                    TaylorSeries.mul!(P_nm[i, j, n, m], secϕ_P_nm[i, j, n, m], cos_ϕ[i, j], ord)
                                    TaylorSeries.mul!(tmp1057[i, j, n, m], secϕ_P_nm[i, j, n, m], sin_ϕ[i, j], ord)
                                    TaylorSeries.mul!(tmp1058[i, j, n, m], tmp1057[i, j, n, m], lnm3[n], ord)
                                    TaylorSeries.mul!(tmp1059[i, j, n - 1, m], secϕ_P_nm[i, j, n - 1, m], lnm4[n, m], ord)
                                    TaylorSeries.add!(cosϕ_dP_nm[i, j, n, m], tmp1058[i, j, n, m], tmp1059[i, j, n - 1, m], ord)
                                end
                            end
                            TaylorSeries.mul!(tmp1061[i, j, 2, 1], P_nm[i, j, 2, 1], lnm6[2], ord)
                            TaylorSeries.mul!(tmp1062[i, j, 1], C21M_t, cos_mλ[i, j, 1], ord)
                            TaylorSeries.mul!(tmp1063[i, j, 1], S21M_t, sin_mλ[i, j, 1], ord)
                            TaylorSeries.add!(tmp1064[i, j, 1], tmp1062[i, j, 1], tmp1063[i, j, 1], ord)
                            TaylorSeries.mul!(tmp1065[i, j, 2, 1], tmp1061[i, j, 2, 1], tmp1064[i, j, 1], ord)
                            TaylorSeries.mul!(tmp1066[i, j, 2, 2], P_nm[i, j, 2, 2], lnm6[2], ord)
                            TaylorSeries.mul!(tmp1067[i, j, 2], C22M_t, cos_mλ[i, j, 2], ord)
                            TaylorSeries.mul!(tmp1068[i, j, 2], S22M_t, sin_mλ[i, j, 2], ord)
                            TaylorSeries.add!(tmp1069[i, j, 2], tmp1067[i, j, 2], tmp1068[i, j, 2], ord)
                            TaylorSeries.mul!(tmp1070[i, j, 2, 2], tmp1066[i, j, 2, 2], tmp1069[i, j, 2], ord)
                            TaylorSeries.add!(tmp1071[i, j, 2, 1], tmp1065[i, j, 2, 1], tmp1070[i, j, 2, 2], ord)
                            TaylorSeries.div!(F_CS_ξ[i, j], tmp1071[i, j, 2, 1], r_p4[i, j], ord)
                            TaylorSeries.mul!(tmp1073[i, j, 2, 1], secϕ_P_nm[i, j, 2, 1], lnm7[1], ord)
                            TaylorSeries.mul!(tmp1074[i, j, 1], S21M_t, cos_mλ[i, j, 1], ord)
                            TaylorSeries.mul!(tmp1075[i, j, 1], C21M_t, sin_mλ[i, j, 1], ord)
                            TaylorSeries.subst!(tmp1076[i, j, 1], tmp1074[i, j, 1], tmp1075[i, j, 1], ord)
                            TaylorSeries.mul!(tmp1077[i, j, 2, 1], tmp1073[i, j, 2, 1], tmp1076[i, j, 1], ord)
                            TaylorSeries.mul!(tmp1078[i, j, 2, 2], secϕ_P_nm[i, j, 2, 2], lnm7[2], ord)
                            TaylorSeries.mul!(tmp1079[i, j, 2], S22M_t, cos_mλ[i, j, 2], ord)
                            TaylorSeries.mul!(tmp1080[i, j, 2], C22M_t, sin_mλ[i, j, 2], ord)
                            TaylorSeries.subst!(tmp1081[i, j, 2], tmp1079[i, j, 2], tmp1080[i, j, 2], ord)
                            TaylorSeries.mul!(tmp1082[i, j, 2, 2], tmp1078[i, j, 2, 2], tmp1081[i, j, 2], ord)
                            TaylorSeries.add!(tmp1083[i, j, 2, 1], tmp1077[i, j, 2, 1], tmp1082[i, j, 2, 2], ord)
                            TaylorSeries.div!(F_CS_η[i, j], tmp1083[i, j, 2, 1], r_p4[i, j], ord)
                            TaylorSeries.mul!(tmp1085[i, j, 1], C21M_t, cos_mλ[i, j, 1], ord)
                            TaylorSeries.mul!(tmp1086[i, j, 1], S21M_t, sin_mλ[i, j, 1], ord)
                            TaylorSeries.add!(tmp1087[i, j, 1], tmp1085[i, j, 1], tmp1086[i, j, 1], ord)
                            TaylorSeries.mul!(tmp1088[i, j, 2, 1], cosϕ_dP_nm[i, j, 2, 1], tmp1087[i, j, 1], ord)
                            TaylorSeries.mul!(tmp1089[i, j, 2], C22M_t, cos_mλ[i, j, 2], ord)
                            TaylorSeries.mul!(tmp1090[i, j, 2], S22M_t, sin_mλ[i, j, 2], ord)
                            TaylorSeries.add!(tmp1091[i, j, 2], tmp1089[i, j, 2], tmp1090[i, j, 2], ord)
                            TaylorSeries.mul!(tmp1092[i, j, 2, 2], cosϕ_dP_nm[i, j, 2, 2], tmp1091[i, j, 2], ord)
                            TaylorSeries.add!(tmp1093[i, j, 2, 1], tmp1088[i, j, 2, 1], tmp1092[i, j, 2, 2], ord)
                            TaylorSeries.div!(F_CS_ζ[i, j], tmp1093[i, j, 2, 1], r_p4[i, j], ord)
                            TaylorSeries.identity!(F_CS_ξ_36[i, j], zero_q_1, ord)
                            TaylorSeries.identity!(F_CS_η_36[i, j], zero_q_1, ord)
                            TaylorSeries.identity!(F_CS_ζ_36[i, j], zero_q_1, ord)
                            for n = 3:n1SEM[mo]
                                for m = 1:n
                                    TaylorSeries.mul!(tmp1095[i, j, n, m], P_nm[i, j, n, m], lnm6[n], ord)
                                    TaylorSeries.mul!(tmp1096[i, j, m], cos_mλ[i, j, m], CM[n, m], ord)
                                    TaylorSeries.mul!(tmp1097[i, j, m], sin_mλ[i, j, m], SM[n, m], ord)
                                    TaylorSeries.add!(tmp1098[i, j, m], tmp1096[i, j, m], tmp1097[i, j, m], ord)
                                    TaylorSeries.mul!(tmp1099[i, j, n, m], tmp1095[i, j, n, m], tmp1098[i, j, m], ord)
                                    TaylorSeries.div!(tmp1100[i, j, n, m], tmp1099[i, j, n, m], temp_rn[i, j, n], ord)
                                    TaylorSeries.add!(temp_CS_ξ[i, j, n, m], tmp1100[i, j, n, m], F_CS_ξ_36[i, j], ord)
                                    TaylorSeries.mul!(tmp1102[i, j, n, m], secϕ_P_nm[i, j, n, m], lnm7[m], ord)
                                    TaylorSeries.mul!(tmp1103[i, j, m], cos_mλ[i, j, m], SM[n, m], ord)
                                    TaylorSeries.mul!(tmp1104[i, j, m], sin_mλ[i, j, m], CM[n, m], ord)
                                    TaylorSeries.subst!(tmp1105[i, j, m], tmp1103[i, j, m], tmp1104[i, j, m], ord)
                                    TaylorSeries.mul!(tmp1106[i, j, n, m], tmp1102[i, j, n, m], tmp1105[i, j, m], ord)
                                    TaylorSeries.div!(tmp1107[i, j, n, m], tmp1106[i, j, n, m], temp_rn[i, j, n], ord)
                                    TaylorSeries.add!(temp_CS_η[i, j, n, m], tmp1107[i, j, n, m], F_CS_η_36[i, j], ord)
                                    TaylorSeries.mul!(tmp1109[i, j, m], cos_mλ[i, j, m], CM[n, m], ord)
                                    TaylorSeries.mul!(tmp1110[i, j, m], sin_mλ[i, j, m], SM[n, m], ord)
                                    TaylorSeries.add!(tmp1111[i, j, m], tmp1109[i, j, m], tmp1110[i, j, m], ord)
                                    TaylorSeries.mul!(tmp1112[i, j, n, m], cosϕ_dP_nm[i, j, n, m], tmp1111[i, j, m], ord)
                                    TaylorSeries.div!(tmp1113[i, j, n, m], tmp1112[i, j, n, m], temp_rn[i, j, n], ord)
                                    TaylorSeries.add!(temp_CS_ζ[i, j, n, m], tmp1113[i, j, n, m], F_CS_ζ_36[i, j], ord)
                                    TaylorSeries.identity!(F_CS_ξ_36[i, j], temp_CS_ξ[i, j, n, m], ord)
                                    TaylorSeries.identity!(F_CS_η_36[i, j], temp_CS_η[i, j, n, m], ord)
                                    TaylorSeries.identity!(F_CS_ζ_36[i, j], temp_CS_ζ[i, j, n, m], ord)
                                end
                            end
                            TaylorSeries.add!(tmp1115[i, j], F_J_ξ[i, j], F_J_ξ_36[i, j], ord)
                            TaylorSeries.add!(tmp1116[i, j], F_CS_ξ[i, j], F_CS_ξ_36[i, j], ord)
                            TaylorSeries.add!(F_JCS_ξ[i, j], tmp1115[i, j], tmp1116[i, j], ord)
                            TaylorSeries.add!(F_JCS_η[i, j], F_CS_η[i, j], F_CS_η_36[i, j], ord)
                            TaylorSeries.add!(tmp1119[i, j], F_J_ζ[i, j], F_J_ζ_36[i, j], ord)
                            TaylorSeries.add!(tmp1120[i, j], F_CS_ζ[i, j], F_CS_ζ_36[i, j], ord)
                            TaylorSeries.add!(F_JCS_ζ[i, j], tmp1119[i, j], tmp1120[i, j], ord)
                        else
                            TaylorSeries.add!(F_JCS_ξ[i, j], F_J_ξ[i, j], F_J_ξ_36[i, j], ord)
                            TaylorSeries.identity!(F_JCS_η[i, j], zero_q_1, ord)
                            TaylorSeries.add!(F_JCS_ζ[i, j], F_J_ζ[i, j], F_J_ζ_36[i, j], ord)
                        end
                        TaylorSeries.mul!(Rb2p[i, j, 1, 1], cos_ϕ[i, j], cos_λ[i, j], ord)
                        TaylorSeries.subst!(Rb2p[i, j, 2, 1], sin_λ[i, j], ord)
                        TaylorSeries.subst!(tmp1126[i, j], sin_ϕ[i, j], ord)
                        TaylorSeries.mul!(Rb2p[i, j, 3, 1], tmp1126[i, j], cos_λ[i, j], ord)
                        TaylorSeries.mul!(Rb2p[i, j, 1, 2], cos_ϕ[i, j], sin_λ[i, j], ord)
                        TaylorSeries.identity!(Rb2p[i, j, 2, 2], cos_λ[i, j], ord)
                        TaylorSeries.subst!(tmp1129[i, j], sin_ϕ[i, j], ord)
                        TaylorSeries.mul!(Rb2p[i, j, 3, 2], tmp1129[i, j], sin_λ[i, j], ord)
                        TaylorSeries.identity!(Rb2p[i, j, 1, 3], sin_ϕ[i, j], ord)
                        TaylorSeries.identity!(Rb2p[i, j, 2, 3], zero_q_1, ord)
                        TaylorSeries.identity!(Rb2p[i, j, 3, 3], cos_ϕ[i, j], ord)
                        TaylorSeries.mul!(tmp1131[i, j, 1, 1], Rb2p[i, j, 1, 1], M_[1, 1, j], ord)
                        TaylorSeries.mul!(tmp1132[i, j, 1, 2], Rb2p[i, j, 1, 2], M_[2, 1, j], ord)
                        TaylorSeries.add!(tmp1133[i, j, 1, 1], tmp1131[i, j, 1, 1], tmp1132[i, j, 1, 2], ord)
                        TaylorSeries.mul!(tmp1134[i, j, 1, 3], Rb2p[i, j, 1, 3], M_[3, 1, j], ord)
                        TaylorSeries.add!(Gc2p[i, j, 1, 1], tmp1133[i, j, 1, 1], tmp1134[i, j, 1, 3], ord)
                        TaylorSeries.mul!(tmp1136[i, j, 2, 1], Rb2p[i, j, 2, 1], M_[1, 1, j], ord)
                        TaylorSeries.mul!(tmp1137[i, j, 2, 2], Rb2p[i, j, 2, 2], M_[2, 1, j], ord)
                        TaylorSeries.add!(tmp1138[i, j, 2, 1], tmp1136[i, j, 2, 1], tmp1137[i, j, 2, 2], ord)
                        TaylorSeries.mul!(tmp1139[i, j, 2, 3], Rb2p[i, j, 2, 3], M_[3, 1, j], ord)
                        TaylorSeries.add!(Gc2p[i, j, 2, 1], tmp1138[i, j, 2, 1], tmp1139[i, j, 2, 3], ord)
                        TaylorSeries.mul!(tmp1141[i, j, 3, 1], Rb2p[i, j, 3, 1], M_[1, 1, j], ord)
                        TaylorSeries.mul!(tmp1142[i, j, 3, 2], Rb2p[i, j, 3, 2], M_[2, 1, j], ord)
                        TaylorSeries.add!(tmp1143[i, j, 3, 1], tmp1141[i, j, 3, 1], tmp1142[i, j, 3, 2], ord)
                        TaylorSeries.mul!(tmp1144[i, j, 3, 3], Rb2p[i, j, 3, 3], M_[3, 1, j], ord)
                        TaylorSeries.add!(Gc2p[i, j, 3, 1], tmp1143[i, j, 3, 1], tmp1144[i, j, 3, 3], ord)
                        TaylorSeries.mul!(tmp1146[i, j, 1, 1], Rb2p[i, j, 1, 1], M_[1, 2, j], ord)
                        TaylorSeries.mul!(tmp1147[i, j, 1, 2], Rb2p[i, j, 1, 2], M_[2, 2, j], ord)
                        TaylorSeries.add!(tmp1148[i, j, 1, 1], tmp1146[i, j, 1, 1], tmp1147[i, j, 1, 2], ord)
                        TaylorSeries.mul!(tmp1149[i, j, 1, 3], Rb2p[i, j, 1, 3], M_[3, 2, j], ord)
                        TaylorSeries.add!(Gc2p[i, j, 1, 2], tmp1148[i, j, 1, 1], tmp1149[i, j, 1, 3], ord)
                        TaylorSeries.mul!(tmp1151[i, j, 2, 1], Rb2p[i, j, 2, 1], M_[1, 2, j], ord)
                        TaylorSeries.mul!(tmp1152[i, j, 2, 2], Rb2p[i, j, 2, 2], M_[2, 2, j], ord)
                        TaylorSeries.add!(tmp1153[i, j, 2, 1], tmp1151[i, j, 2, 1], tmp1152[i, j, 2, 2], ord)
                        TaylorSeries.mul!(tmp1154[i, j, 2, 3], Rb2p[i, j, 2, 3], M_[3, 2, j], ord)
                        TaylorSeries.add!(Gc2p[i, j, 2, 2], tmp1153[i, j, 2, 1], tmp1154[i, j, 2, 3], ord)
                        TaylorSeries.mul!(tmp1156[i, j, 3, 1], Rb2p[i, j, 3, 1], M_[1, 2, j], ord)
                        TaylorSeries.mul!(tmp1157[i, j, 3, 2], Rb2p[i, j, 3, 2], M_[2, 2, j], ord)
                        TaylorSeries.add!(tmp1158[i, j, 3, 1], tmp1156[i, j, 3, 1], tmp1157[i, j, 3, 2], ord)
                        TaylorSeries.mul!(tmp1159[i, j, 3, 3], Rb2p[i, j, 3, 3], M_[3, 2, j], ord)
                        TaylorSeries.add!(Gc2p[i, j, 3, 2], tmp1158[i, j, 3, 1], tmp1159[i, j, 3, 3], ord)
                        TaylorSeries.mul!(tmp1161[i, j, 1, 1], Rb2p[i, j, 1, 1], M_[1, 3, j], ord)
                        TaylorSeries.mul!(tmp1162[i, j, 1, 2], Rb2p[i, j, 1, 2], M_[2, 3, j], ord)
                        TaylorSeries.add!(tmp1163[i, j, 1, 1], tmp1161[i, j, 1, 1], tmp1162[i, j, 1, 2], ord)
                        TaylorSeries.mul!(tmp1164[i, j, 1, 3], Rb2p[i, j, 1, 3], M_[3, 3, j], ord)
                        TaylorSeries.add!(Gc2p[i, j, 1, 3], tmp1163[i, j, 1, 1], tmp1164[i, j, 1, 3], ord)
                        TaylorSeries.mul!(tmp1166[i, j, 2, 1], Rb2p[i, j, 2, 1], M_[1, 3, j], ord)
                        TaylorSeries.mul!(tmp1167[i, j, 2, 2], Rb2p[i, j, 2, 2], M_[2, 3, j], ord)
                        TaylorSeries.add!(tmp1168[i, j, 2, 1], tmp1166[i, j, 2, 1], tmp1167[i, j, 2, 2], ord)
                        TaylorSeries.mul!(tmp1169[i, j, 2, 3], Rb2p[i, j, 2, 3], M_[3, 3, j], ord)
                        TaylorSeries.add!(Gc2p[i, j, 2, 3], tmp1168[i, j, 2, 1], tmp1169[i, j, 2, 3], ord)
                        TaylorSeries.mul!(tmp1171[i, j, 3, 1], Rb2p[i, j, 3, 1], M_[1, 3, j], ord)
                        TaylorSeries.mul!(tmp1172[i, j, 3, 2], Rb2p[i, j, 3, 2], M_[2, 3, j], ord)
                        TaylorSeries.add!(tmp1173[i, j, 3, 1], tmp1171[i, j, 3, 1], tmp1172[i, j, 3, 2], ord)
                        TaylorSeries.mul!(tmp1174[i, j, 3, 3], Rb2p[i, j, 3, 3], M_[3, 3, j], ord)
                        TaylorSeries.add!(Gc2p[i, j, 3, 3], tmp1173[i, j, 3, 1], tmp1174[i, j, 3, 3], ord)
                        TaylorSeries.mul!(tmp1176[i, j, 1, 1], F_JCS_ξ[i, j], Gc2p[i, j, 1, 1], ord)
                        TaylorSeries.mul!(tmp1177[i, j, 2, 1], F_JCS_η[i, j], Gc2p[i, j, 2, 1], ord)
                        TaylorSeries.add!(tmp1178[i, j, 1, 1], tmp1176[i, j, 1, 1], tmp1177[i, j, 2, 1], ord)
                        TaylorSeries.mul!(tmp1179[i, j, 3, 1], F_JCS_ζ[i, j], Gc2p[i, j, 3, 1], ord)
                        TaylorSeries.add!(F_JCS_x[i, j], tmp1178[i, j, 1, 1], tmp1179[i, j, 3, 1], ord)
                        TaylorSeries.mul!(tmp1181[i, j, 1, 2], F_JCS_ξ[i, j], Gc2p[i, j, 1, 2], ord)
                        TaylorSeries.mul!(tmp1182[i, j, 2, 2], F_JCS_η[i, j], Gc2p[i, j, 2, 2], ord)
                        TaylorSeries.add!(tmp1183[i, j, 1, 2], tmp1181[i, j, 1, 2], tmp1182[i, j, 2, 2], ord)
                        TaylorSeries.mul!(tmp1184[i, j, 3, 2], F_JCS_ζ[i, j], Gc2p[i, j, 3, 2], ord)
                        TaylorSeries.add!(F_JCS_y[i, j], tmp1183[i, j, 1, 2], tmp1184[i, j, 3, 2], ord)
                        TaylorSeries.mul!(tmp1186[i, j, 1, 3], F_JCS_ξ[i, j], Gc2p[i, j, 1, 3], ord)
                        TaylorSeries.mul!(tmp1187[i, j, 2, 3], F_JCS_η[i, j], Gc2p[i, j, 2, 3], ord)
                        TaylorSeries.add!(tmp1188[i, j, 1, 3], tmp1186[i, j, 1, 3], tmp1187[i, j, 2, 3], ord)
                        TaylorSeries.mul!(tmp1189[i, j, 3, 3], F_JCS_ζ[i, j], Gc2p[i, j, 3, 3], ord)
                        TaylorSeries.add!(F_JCS_z[i, j], tmp1188[i, j, 1, 3], tmp1189[i, j, 3, 3], ord)
                    end
                end
            end
        end
        for j = 1:N
            for i = 1:N
                if i == j
                    continue
                else
                    if UJ_interaction[i, j]
                        TaylorSeries.mul!(tmp1191[i, j], μ[i], F_JCS_x[i, j], ord)
                        TaylorSeries.subst!(temp_accX_j[i, j], accX[j], tmp1191[i, j], ord)
                        TaylorSeries.identity!(accX[j], temp_accX_j[i, j], ord)
                        TaylorSeries.mul!(tmp1193[i, j], μ[i], F_JCS_y[i, j], ord)
                        TaylorSeries.subst!(temp_accY_j[i, j], accY[j], tmp1193[i, j], ord)
                        TaylorSeries.identity!(accY[j], temp_accY_j[i, j], ord)
                        TaylorSeries.mul!(tmp1195[i, j], μ[i], F_JCS_z[i, j], ord)
                        TaylorSeries.subst!(temp_accZ_j[i, j], accZ[j], tmp1195[i, j], ord)
                        TaylorSeries.identity!(accZ[j], temp_accZ_j[i, j], ord)
                        TaylorSeries.mul!(tmp1197[i, j], μ[j], F_JCS_x[i, j], ord)
                        TaylorSeries.add!(temp_accX_i[i, j], accX[i], tmp1197[i, j], ord)
                        TaylorSeries.identity!(accX[i], temp_accX_i[i, j], ord)
                        TaylorSeries.mul!(tmp1199[i, j], μ[j], F_JCS_y[i, j], ord)
                        TaylorSeries.add!(temp_accY_i[i, j], accY[i], tmp1199[i, j], ord)
                        TaylorSeries.identity!(accY[i], temp_accY_i[i, j], ord)
                        TaylorSeries.mul!(tmp1201[i, j], μ[j], F_JCS_z[i, j], ord)
                        TaylorSeries.add!(temp_accZ_i[i, j], accZ[i], tmp1201[i, j], ord)
                        TaylorSeries.identity!(accZ[i], temp_accZ_i[i, j], ord)
                    end
                end
            end
        end
        for j = 1:N
            for i = 1:N
                if i == j
                    continue
                else
                    TaylorSeries.mul!(_4ϕj[i, j], 4, newtonianNb_Potential[j], ord)
                    TaylorSeries.add!(ϕi_plus_4ϕj[i, j], newtonianNb_Potential[i], _4ϕj[i, j], ord)
                    TaylorSeries.mul!(tmp1207[i], 2, v2[i], ord)
                    TaylorSeries.add!(tmp1208[j], v2[j], tmp1207[i], ord)
                    TaylorSeries.mul!(tmp1210[i, j], 4, vi_dot_vj[i, j], ord)
                    TaylorSeries.subst!(sj2_plus_2si2_minus_4vivj[i, j], tmp1208[j], tmp1210[i, j], ord)
                    TaylorSeries.subst!(ϕs_and_vs[i, j], sj2_plus_2si2_minus_4vivj[i, j], ϕi_plus_4ϕj[i, j], ord)
                    TaylorSeries.mul!(Xij_t_Ui[i, j], X[i, j], dq[3i - 2], ord)
                    TaylorSeries.mul!(Yij_t_Vi[i, j], Y[i, j], dq[3i - 1], ord)
                    TaylorSeries.mul!(Zij_t_Wi[i, j], Z[i, j], dq[3i], ord)
                    TaylorSeries.add!(tmp1216[i, j], Xij_t_Ui[i, j], Yij_t_Vi[i, j], ord)
                    TaylorSeries.add!(Rij_dot_Vi[i, j], tmp1216[i, j], Zij_t_Wi[i, j], ord)
                    TaylorSeries.pow!(tmp1219[i, j], Rij_dot_Vi[i, j], 2, ord)
                    TaylorSeries.div!(pn1t7[i, j], tmp1219[i, j], r_p2[i, j], ord)
                    TaylorSeries.mul!(tmp1222[i, j], 1.5, pn1t7[i, j], ord)
                    TaylorSeries.subst!(pn1t2_7[i, j], ϕs_and_vs[i, j], tmp1222[i, j], ord)
                    TaylorSeries.add!(pn1t1_7[i, j], c_p2, pn1t2_7[i, j], ord)
                    for k = 1:postnewton_iter
                        TaylorSeries.identity!(pn1[i, j, k], zero_q_1, ord)
                        TaylorSeries.identity!(X_t_pn1[i, j, k], zero_q_1, ord)
                        TaylorSeries.identity!(Y_t_pn1[i, j, k], zero_q_1, ord)
                        TaylorSeries.identity!(Z_t_pn1[i, j, k], zero_q_1, ord)
                        TaylorSeries.identity!(pNX_t_pn3[i, j, k], zero_q_1, ord)
                        TaylorSeries.identity!(pNY_t_pn3[i, j, k], zero_q_1, ord)
                        TaylorSeries.identity!(pNZ_t_pn3[i, j, k], zero_q_1, ord)
                        TaylorSeries.identity!(pNX_t_X[i, j, k], zero_q_1, ord)
                        TaylorSeries.identity!(pNY_t_Y[i, j, k], zero_q_1, ord)
                        TaylorSeries.identity!(pNZ_t_Z[i, j, k], zero_q_1, ord)
                    end
                end
            end
            TaylorSeries.identity!(postNewtonX[j, 1], newtonX[j], ord)
            TaylorSeries.identity!(postNewtonY[j, 1], newtonY[j], ord)
            TaylorSeries.identity!(postNewtonZ[j, 1], newtonZ[j], ord)
            for k = 1:postnewton_iter
                TaylorSeries.identity!(pntempX[j, k], zero_q_1, ord)
                TaylorSeries.identity!(pntempY[j, k], zero_q_1, ord)
                TaylorSeries.identity!(pntempZ[j, k], zero_q_1, ord)
            end
        end
        for k = 1:postnewton_iter
            for j = 1:N
                for i = 1:N
                    if i == j
                        continue
                    else
                        TaylorSeries.mul!(pNX_t_X[i, j, k], postNewtonX[i, k], X[i, j], ord)
                        TaylorSeries.mul!(pNY_t_Y[i, j, k], postNewtonY[i, k], Y[i, j], ord)
                        TaylorSeries.mul!(pNZ_t_Z[i, j, k], postNewtonZ[i, k], Z[i, j], ord)
                        TaylorSeries.add!(tmp1229[i, j, k], pNX_t_X[i, j, k], pNY_t_Y[i, j, k], ord)
                        TaylorSeries.add!(tmp1230[i, j, k], tmp1229[i, j, k], pNZ_t_Z[i, j, k], ord)
                        TaylorSeries.mul!(tmp1231[i, j, k], 0.5, tmp1230[i, j, k], ord)
                        TaylorSeries.add!(pn1[i, j, k], pn1t1_7[i, j], tmp1231[i, j, k], ord)
                        TaylorSeries.mul!(X_t_pn1[i, j, k], newton_acc_X[i, j], pn1[i, j, k], ord)
                        TaylorSeries.mul!(Y_t_pn1[i, j, k], newton_acc_Y[i, j], pn1[i, j, k], ord)
                        TaylorSeries.mul!(Z_t_pn1[i, j, k], newton_acc_Z[i, j], pn1[i, j, k], ord)
                        TaylorSeries.mul!(pNX_t_pn3[i, j, k], postNewtonX[i, k], pn3[i, j], ord)
                        TaylorSeries.mul!(pNY_t_pn3[i, j, k], postNewtonY[i, k], pn3[i, j], ord)
                        TaylorSeries.mul!(pNZ_t_pn3[i, j, k], postNewtonZ[i, k], pn3[i, j], ord)
                        TaylorSeries.add!(tmp1239[i, j, k], U_t_pn2[i, j], pNX_t_pn3[i, j, k], ord)
                        TaylorSeries.add!(termpnx[i, j, k], X_t_pn1[i, j, k], tmp1239[i, j, k], ord)
                        TaylorSeries.add!(sumpnx[i, j, k], pntempX[j, k], termpnx[i, j, k], ord)
                        TaylorSeries.identity!(pntempX[j, k], sumpnx[i, j, k], ord)
                        TaylorSeries.add!(tmp1242[i, j, k], V_t_pn2[i, j], pNY_t_pn3[i, j, k], ord)
                        TaylorSeries.add!(termpny[i, j, k], Y_t_pn1[i, j, k], tmp1242[i, j, k], ord)
                        TaylorSeries.add!(sumpny[i, j, k], pntempY[j, k], termpny[i, j, k], ord)
                        TaylorSeries.identity!(pntempY[j, k], sumpny[i, j, k], ord)
                        TaylorSeries.add!(tmp1245[i, j, k], W_t_pn2[i, j], pNZ_t_pn3[i, j, k], ord)
                        TaylorSeries.add!(termpnz[i, j, k], Z_t_pn1[i, j, k], tmp1245[i, j, k], ord)
                        TaylorSeries.add!(sumpnz[i, j, k], pntempZ[j, k], termpnz[i, j, k], ord)
                        TaylorSeries.identity!(pntempZ[j, k], sumpnz[i, j, k], ord)
                    end
                end
                TaylorSeries.mul!(postNewtonX[j, k + 1], pntempX[j, k], c_m2, ord)
                TaylorSeries.mul!(postNewtonY[j, k + 1], pntempY[j, k], c_m2, ord)
                TaylorSeries.mul!(postNewtonZ[j, k + 1], pntempZ[j, k], c_m2, ord)
            end
        end
        for i = 1:N
            TaylorSeries.add!(dq[3 * (N + i) - 2], postNewtonX[i, postnewton_iter + 1], accX[i], ord)
            TaylorSeries.add!(dq[3 * (N + i) - 1], postNewtonY[i, postnewton_iter + 1], accY[i], ord)
            TaylorSeries.add!(dq[3 * (N + i)], postNewtonZ[i, postnewton_iter + 1], accZ[i], ord)
        end
        for __idx = eachindex(q)
            (q[__idx]).coeffs[ordnext + 1] = (dq[__idx]).coeffs[ordnext] / ordnext
        end
    end
    return nothing
end

function TaylorIntegration.jetcoeffs!(::Val{NBP_pN_A_J23E_J23M_J2S_threads!}, t::Taylor1{_T}, q::AbstractVector{Taylor1{_S}}, dq::AbstractVector{Taylor1{_S}}, params) where {_T <: Real, _S <: Number}
    order = t.order
    local (N, S, eulang_de430_, jd0) = params
    local N_ext = 11
    local eulang_t = eulang_de430_((t + (jd0 - J2000)) * daysec)
    local postnewton_iter = 1
    local c_p2 = 29979.063823897606
    local c_m2 = 3.3356611996764786e-5
    local zero_q_1 = zero(q[1])
    local one_t = one(t)
    X = Array{Taylor1{S}}(undef, N, N)
    Y = Array{Taylor1{S}}(undef, N, N)
    Z = Array{Taylor1{S}}(undef, N, N)
    r_p2 = Array{Taylor1{S}}(undef, N, N)
    r_p3d2 = Array{Taylor1{S}}(undef, N, N)
    r_p7d2 = Array{Taylor1{S}}(undef, N, N)
    newtonX = Array{Taylor1{S}}(undef, N)
    newtonY = Array{Taylor1{S}}(undef, N)
    newtonZ = Array{Taylor1{S}}(undef, N)
    newtonianCoeff = Array{Taylor1{S}}(undef, N, N)
    U = Array{Taylor1{S}}(undef, N, N)
    V = Array{Taylor1{S}}(undef, N, N)
    W = Array{Taylor1{S}}(undef, N, N)
    _4U_m_3X = Array{Taylor1{S}}(undef, N, N)
    _4V_m_3Y = Array{Taylor1{S}}(undef, N, N)
    _4W_m_3Z = Array{Taylor1{S}}(undef, N, N)
    UU = Array{Taylor1{S}}(undef, N, N)
    VV = Array{Taylor1{S}}(undef, N, N)
    WW = Array{Taylor1{S}}(undef, N, N)
    r_p1d2 = Array{Taylor1{S}}(undef, N, N)
    newtonianNb_Potential = Array{Taylor1{S}}(undef, N)
    newtonian1b_Potential = Array{Taylor1{S}}(undef, N, N)
    newton_acc_X = Array{Taylor1{S}}(undef, N, N)
    newton_acc_Y = Array{Taylor1{S}}(undef, N, N)
    newton_acc_Z = Array{Taylor1{S}}(undef, N, N)
    v2 = Array{Taylor1{S}}(undef, N)
    vi_dot_vj = Array{Taylor1{S}}(undef, N, N)
    pn2 = Array{Taylor1{S}}(undef, N, N)
    pn3 = Array{Taylor1{S}}(undef, N, N)
    _4ϕj = Array{Taylor1{S}}(undef, N, N)
    ϕi_plus_4ϕj = Array{Taylor1{S}}(undef, N, N)
    sj2_plus_2si2_minus_4vivj = Array{Taylor1{S}}(undef, N, N)
    ϕs_and_vs = Array{Taylor1{S}}(undef, N, N)
    U_t_pn2 = Array{Taylor1{S}}(undef, N, N)
    V_t_pn2 = Array{Taylor1{S}}(undef, N, N)
    W_t_pn2 = Array{Taylor1{S}}(undef, N, N)
    pn1t1_7 = Array{Taylor1{S}}(undef, N, N)
    pntempX = Array{Taylor1{S}}(undef, N, postnewton_iter)
    pntempY = Array{Taylor1{S}}(undef, N, postnewton_iter)
    pntempZ = Array{Taylor1{S}}(undef, N, postnewton_iter)
    pn1 = Array{Taylor1{S}}(undef, N, N, postnewton_iter)
    X_t_pn1 = Array{Taylor1{S}}(undef, N, N, postnewton_iter)
    Y_t_pn1 = Array{Taylor1{S}}(undef, N, N, postnewton_iter)
    Z_t_pn1 = Array{Taylor1{S}}(undef, N, N, postnewton_iter)
    pNX_t_pn3 = Array{Taylor1{S}}(undef, N, N, postnewton_iter)
    pNY_t_pn3 = Array{Taylor1{S}}(undef, N, N, postnewton_iter)
    pNZ_t_pn3 = Array{Taylor1{S}}(undef, N, N, postnewton_iter)
    pNX_t_X = Array{Taylor1{S}}(undef, N, N, postnewton_iter)
    pNY_t_Y = Array{Taylor1{S}}(undef, N, N, postnewton_iter)
    pNZ_t_Z = Array{Taylor1{S}}(undef, N, N, postnewton_iter)
    postNewtonX = Array{Taylor1{S}}(undef, N, postnewton_iter + 1)
    postNewtonY = Array{Taylor1{S}}(undef, N, postnewton_iter + 1)
    postNewtonZ = Array{Taylor1{S}}(undef, N, postnewton_iter + 1)
    X_bf_1 = Array{Taylor1{S}}(undef, N_ext, N_ext)
    Y_bf_1 = Array{Taylor1{S}}(undef, N_ext, N_ext)
    Z_bf_1 = Array{Taylor1{S}}(undef, N_ext, N_ext)
    X_bf_2 = Array{Taylor1{S}}(undef, N_ext, N_ext)
    Y_bf_2 = Array{Taylor1{S}}(undef, N_ext, N_ext)
    Z_bf_2 = Array{Taylor1{S}}(undef, N_ext, N_ext)
    X_bf_3 = Array{Taylor1{S}}(undef, N_ext, N_ext)
    Y_bf_3 = Array{Taylor1{S}}(undef, N_ext, N_ext)
    Z_bf_3 = Array{Taylor1{S}}(undef, N_ext, N_ext)
    X_bf = Array{Taylor1{S}}(undef, N_ext, N_ext)
    Y_bf = Array{Taylor1{S}}(undef, N_ext, N_ext)
    Z_bf = Array{Taylor1{S}}(undef, N_ext, N_ext)
    F_JCS_x = Array{Taylor1{S}}(undef, N_ext, N_ext)
    F_JCS_y = Array{Taylor1{S}}(undef, N_ext, N_ext)
    F_JCS_z = Array{Taylor1{S}}(undef, N_ext, N_ext)
    temp_accX_j = Array{Taylor1{S}}(undef, N_ext, N_ext)
    temp_accY_j = Array{Taylor1{S}}(undef, N_ext, N_ext)
    temp_accZ_j = Array{Taylor1{S}}(undef, N_ext, N_ext)
    temp_accX_i = Array{Taylor1{S}}(undef, N_ext, N_ext)
    temp_accY_i = Array{Taylor1{S}}(undef, N_ext, N_ext)
    temp_accZ_i = Array{Taylor1{S}}(undef, N_ext, N_ext)
    sin_ϕ = Array{Taylor1{S}}(undef, N_ext, N_ext)
    cos_ϕ = Array{Taylor1{S}}(undef, N_ext, N_ext)
    sin_λ = Array{Taylor1{S}}(undef, N_ext, N_ext)
    cos_λ = Array{Taylor1{S}}(undef, N_ext, N_ext)
    r_xy = Array{Taylor1{S}}(undef, N_ext, N_ext)
    r_p4 = Array{Taylor1{S}}(undef, N_ext, N_ext)
    P_n = Array{Taylor1{S}}(undef, N_ext, N_ext, maximum(n1SEM) + 1)
    dP_n = Array{Taylor1{S}}(undef, N_ext, N_ext, maximum(n1SEM) + 1)
    temp_fjξ = Array{Taylor1{S}}(undef, N_ext, N_ext, maximum(n1SEM) + 1)
    temp_fjζ = Array{Taylor1{S}}(undef, N_ext, N_ext, maximum(n1SEM) + 1)
    temp_rn = Array{Taylor1{S}}(undef, N_ext, N_ext, maximum(n1SEM) + 1)
    F_CS_ξ_36 = Array{Taylor1{S}}(undef, N_ext, N_ext)
    F_CS_η_36 = Array{Taylor1{S}}(undef, N_ext, N_ext)
    F_CS_ζ_36 = Array{Taylor1{S}}(undef, N_ext, N_ext)
    F_J_ξ_36 = Array{Taylor1{S}}(undef, N_ext, N_ext)
    F_J_ζ_36 = Array{Taylor1{S}}(undef, N_ext, N_ext)
    sin_mλ = Array{Taylor1{S}}(undef, N_ext, N_ext, n1SEM[mo])
    cos_mλ = Array{Taylor1{S}}(undef, N_ext, N_ext, n1SEM[mo])
    secϕ_P_nm = Array{Taylor1{S}}(undef, N_ext, N_ext, n1SEM[mo] + 1, n1SEM[mo] + 1)
    P_nm = Array{Taylor1{S}}(undef, N_ext, N_ext, n1SEM[mo] + 1, n1SEM[mo] + 1)
    cosϕ_dP_nm = Array{Taylor1{S}}(undef, N_ext, N_ext, n1SEM[mo] + 1, n1SEM[mo] + 1)
    F_J_ξ = Array{Taylor1{S}}(undef, N_ext, N_ext)
    F_J_η = Array{Taylor1{S}}(undef, N_ext, N_ext)
    F_J_ζ = Array{Taylor1{S}}(undef, N_ext, N_ext)
    F_CS_ξ = Array{Taylor1{S}}(undef, N_ext, N_ext)
    F_CS_η = Array{Taylor1{S}}(undef, N_ext, N_ext)
    F_CS_ζ = Array{Taylor1{S}}(undef, N_ext, N_ext)
    F_JCS_ξ = Array{Taylor1{S}}(undef, N_ext, N_ext)
    F_JCS_η = Array{Taylor1{S}}(undef, N_ext, N_ext)
    F_JCS_ζ = Array{Taylor1{S}}(undef, N_ext, N_ext)
    Rb2p = Array{Taylor1{S}}(undef, N_ext, N_ext, 3, 3)
    Gc2p = Array{Taylor1{S}}(undef, N_ext, N_ext, 3, 3)
    accX = Array{Taylor1{S}}(undef, N_ext)
    accY = Array{Taylor1{S}}(undef, N_ext)
    accZ = Array{Taylor1{S}}(undef, N_ext)
    local dsj2k = t + (jd0 - 2.451545e6)
    local αs = deg2rad(α_p_sun * one_t)
    local δs = deg2rad(δ_p_sun * one_t)
    local αm = eulang_t[1] - pi / 2
    local δm = pi / 2 - eulang_t[2]
    local Wm = eulang_t[3]
    local M_ = Array{Taylor1{S}}(undef, 3, 3, 5)
    local M_[:, :, ea] = c2t_jpl_de430(dsj2k)
    local M_[:, :, su] = pole_rotation(αs, δs)
    local M_[:, :, mo] = pole_rotation(αm, δm, Wm)
    local ITM_t = ITM_und .* one_t
    local fact_num = -4.5257273867882326e-36
    local fact1_jsem = [(2n - 1) / n for n = 1:maximum(n1SEM)]
    local fact2_jsem = [(n - 1) / n for n = 1:maximum(n1SEM)]
    local fact3_jsem = [n for n = 1:maximum(n1SEM)]
    local fact4_jsem = [n + 1 for n = 1:maximum(n1SEM)]
    local fact5_jsem = [n + 2 for n = 1:maximum(n1SEM)]
    local lnm1 = [(2n - 1) / (n - m) for n = 1:6, m = 1:6]
    local lnm2 = [-(((n + m) - 1)) / (n - m) for n = 1:6, m = 1:6]
    local lnm3 = [-n for n = 1:6]
    local lnm4 = [n + m for n = 1:6, m = 1:6]
    local lnm5 = [2n - 1 for n = 1:6]
    local lnm6 = [-((n + 1)) for n = 1:6]
    local lnm7 = [m for m = 1:6]
    local RE_au = RE / au
    local J2E_t = (J2E + J2EDOT * (dsj2k / yr)) * RE_au ^ 2
    local J2S_t = JSEM[su, 2] * one_t
    local J2_t = Array{Taylor1{S}}(undef, 5)
    J2_t[su] = Taylor1(identity(constant_term(J2S_t)), order)
    J2_t[ea] = Taylor1(identity(constant_term(J2E_t)), order)
    #= REPL[4]:180 =# Threads.@threads for j = 1:N
            newtonX[j] = Taylor1(identity(constant_term(zero_q_1)), order)
            newtonY[j] = Taylor1(identity(constant_term(zero_q_1)), order)
            newtonZ[j] = Taylor1(identity(constant_term(zero_q_1)), order)
            newtonianNb_Potential[j] = Taylor1(identity(constant_term(zero_q_1)), order)
            dq[3j - 2] = Taylor1(identity(constant_term(q[3 * (N + j) - 2])), order)
            dq[3j - 1] = Taylor1(identity(constant_term(q[3 * (N + j) - 1])), order)
            dq[3j] = Taylor1(identity(constant_term(q[3 * (N + j)])), order)
        end
    #= REPL[4]:190 =# Threads.@threads for j = 1:N_ext
            accX[j] = Taylor1(identity(constant_term(zero_q_1)), order)
            accY[j] = Taylor1(identity(constant_term(zero_q_1)), order)
            accZ[j] = Taylor1(identity(constant_term(zero_q_1)), order)
        end
    tmp1890 = Array{Taylor1{_S}}(undef, size(dq))
    tmp1890 .= Taylor1(zero(_S), order)
    tmp1892 = Array{Taylor1{_S}}(undef, size(dq))
    tmp1892 .= Taylor1(zero(_S), order)
    tmp1895 = Array{Taylor1{_S}}(undef, size(dq))
    tmp1895 .= Taylor1(zero(_S), order)
    tmp1897 = Array{Taylor1{_S}}(undef, size(dq))
    tmp1897 .= Taylor1(zero(_S), order)
    tmp1900 = Array{Taylor1{_S}}(undef, size(dq))
    tmp1900 .= Taylor1(zero(_S), order)
    tmp1902 = Array{Taylor1{_S}}(undef, size(dq))
    tmp1902 .= Taylor1(zero(_S), order)
    pn2x = Array{Taylor1{_S}}(undef, size(X))
    pn2x .= Taylor1(zero(_S), order)
    pn2y = Array{Taylor1{_S}}(undef, size(Y))
    pn2y .= Taylor1(zero(_S), order)
    pn2z = Array{Taylor1{_S}}(undef, size(Z))
    pn2z .= Taylor1(zero(_S), order)
    tmp1910 = Array{Taylor1{_S}}(undef, size(UU))
    tmp1910 .= Taylor1(zero(_S), order)
    tmp1913 = Array{Taylor1{_S}}(undef, size(X))
    tmp1913 .= Taylor1(zero(_S), order)
    tmp1915 = Array{Taylor1{_S}}(undef, size(Y))
    tmp1915 .= Taylor1(zero(_S), order)
    tmp1916 = Array{Taylor1{_S}}(undef, size(tmp1913))
    tmp1916 .= Taylor1(zero(_S), order)
    tmp1918 = Array{Taylor1{_S}}(undef, size(Z))
    tmp1918 .= Taylor1(zero(_S), order)
    tmp1926 = Array{Taylor1{_S}}(undef, size(pn2x))
    tmp1926 .= Taylor1(zero(_S), order)
    tmp1927 = Array{Taylor1{_S}}(undef, size(tmp1926))
    tmp1927 .= Taylor1(zero(_S), order)
    tmp1938 = Array{Taylor1{_S}}(undef, size(X))
    tmp1938 .= Taylor1(zero(_S), order)
    temp_001 = Array{Taylor1{_S}}(undef, size(tmp1938))
    temp_001 .= Taylor1(zero(_S), order)
    tmp1940 = Array{Taylor1{_S}}(undef, size(Y))
    tmp1940 .= Taylor1(zero(_S), order)
    temp_002 = Array{Taylor1{_S}}(undef, size(tmp1940))
    temp_002 .= Taylor1(zero(_S), order)
    tmp1942 = Array{Taylor1{_S}}(undef, size(Z))
    tmp1942 .= Taylor1(zero(_S), order)
    temp_003 = Array{Taylor1{_S}}(undef, size(tmp1942))
    temp_003 .= Taylor1(zero(_S), order)
    temp_004 = Array{Taylor1{_S}}(undef, size(newtonian1b_Potential))
    temp_004 .= Taylor1(zero(_S), order)
    tmp1946 = Array{Taylor1{_S}}(undef, size(dq))
    tmp1946 .= Taylor1(zero(_S), order)
    tmp1948 = Array{Taylor1{_S}}(undef, size(dq))
    tmp1948 .= Taylor1(zero(_S), order)
    tmp1949 = Array{Taylor1{_S}}(undef, size(tmp1946))
    tmp1949 .= Taylor1(zero(_S), order)
    tmp1951 = Array{Taylor1{_S}}(undef, size(dq))
    tmp1951 .= Taylor1(zero(_S), order)
    #= REPL[4]:197 =# Threads.@threads for j = 1:N
            for i = 1:N
                if i == j
                    continue
                else
                    X[i, j] = Taylor1(constant_term(q[3i - 2]) - constant_term(q[3j - 2]), order)
                    Y[i, j] = Taylor1(constant_term(q[3i - 1]) - constant_term(q[3j - 1]), order)
                    Z[i, j] = Taylor1(constant_term(q[3i]) - constant_term(q[3j]), order)
                    U[i, j] = Taylor1(constant_term(dq[3i - 2]) - constant_term(dq[3j - 2]), order)
                    V[i, j] = Taylor1(constant_term(dq[3i - 1]) - constant_term(dq[3j - 1]), order)
                    W[i, j] = Taylor1(constant_term(dq[3i]) - constant_term(dq[3j]), order)
                    tmp1890[3j - 2] = Taylor1(constant_term(4) * constant_term(dq[3j - 2]), order)
                    tmp1892[3i - 2] = Taylor1(constant_term(3) * constant_term(dq[3i - 2]), order)
                    _4U_m_3X[i, j] = Taylor1(constant_term(tmp1890[3j - 2]) - constant_term(tmp1892[3i - 2]), order)
                    tmp1895[3j - 1] = Taylor1(constant_term(4) * constant_term(dq[3j - 1]), order)
                    tmp1897[3i - 1] = Taylor1(constant_term(3) * constant_term(dq[3i - 1]), order)
                    _4V_m_3Y[i, j] = Taylor1(constant_term(tmp1895[3j - 1]) - constant_term(tmp1897[3i - 1]), order)
                    tmp1900[3j] = Taylor1(constant_term(4) * constant_term(dq[3j]), order)
                    tmp1902[3i] = Taylor1(constant_term(3) * constant_term(dq[3i]), order)
                    _4W_m_3Z[i, j] = Taylor1(constant_term(tmp1900[3j]) - constant_term(tmp1902[3i]), order)
                    pn2x[i, j] = Taylor1(constant_term(X[i, j]) * constant_term(_4U_m_3X[i, j]), order)
                    pn2y[i, j] = Taylor1(constant_term(Y[i, j]) * constant_term(_4V_m_3Y[i, j]), order)
                    pn2z[i, j] = Taylor1(constant_term(Z[i, j]) * constant_term(_4W_m_3Z[i, j]), order)
                    UU[i, j] = Taylor1(constant_term(dq[3i - 2]) * constant_term(dq[3j - 2]), order)
                    VV[i, j] = Taylor1(constant_term(dq[3i - 1]) * constant_term(dq[3j - 1]), order)
                    WW[i, j] = Taylor1(constant_term(dq[3i]) * constant_term(dq[3j]), order)
                    tmp1910[i, j] = Taylor1(constant_term(UU[i, j]) + constant_term(VV[i, j]), order)
                    vi_dot_vj[i, j] = Taylor1(constant_term(tmp1910[i, j]) + constant_term(WW[i, j]), order)
                    tmp1913[i, j] = Taylor1(constant_term(X[i, j]) ^ constant_term(2), order)
                    tmp1915[i, j] = Taylor1(constant_term(Y[i, j]) ^ constant_term(2), order)
                    tmp1916[i, j] = Taylor1(constant_term(tmp1913[i, j]) + constant_term(tmp1915[i, j]), order)
                    tmp1918[i, j] = Taylor1(constant_term(Z[i, j]) ^ constant_term(2), order)
                    r_p2[i, j] = Taylor1(constant_term(tmp1916[i, j]) + constant_term(tmp1918[i, j]), order)
                    r_p1d2[i, j] = Taylor1(sqrt(constant_term(r_p2[i, j])), order)
                    r_p3d2[i, j] = Taylor1(constant_term(r_p2[i, j]) ^ constant_term(1.5), order)
                    r_p7d2[i, j] = Taylor1(constant_term(r_p2[i, j]) ^ constant_term(3.5), order)
                    newtonianCoeff[i, j] = Taylor1(constant_term(μ[i]) / constant_term(r_p3d2[i, j]), order)
                    tmp1926[i, j] = Taylor1(constant_term(pn2x[i, j]) + constant_term(pn2y[i, j]), order)
                    tmp1927[i, j] = Taylor1(constant_term(tmp1926[i, j]) + constant_term(pn2z[i, j]), order)
                    pn2[i, j] = Taylor1(constant_term(newtonianCoeff[i, j]) * constant_term(tmp1927[i, j]), order)
                    newton_acc_X[i, j] = Taylor1(constant_term(X[i, j]) * constant_term(newtonianCoeff[i, j]), order)
                    newton_acc_Y[i, j] = Taylor1(constant_term(Y[i, j]) * constant_term(newtonianCoeff[i, j]), order)
                    newton_acc_Z[i, j] = Taylor1(constant_term(Z[i, j]) * constant_term(newtonianCoeff[i, j]), order)
                    newtonian1b_Potential[i, j] = Taylor1(constant_term(μ[i]) / constant_term(r_p1d2[i, j]), order)
                    pn3[i, j] = Taylor1(constant_term(3.5) * constant_term(newtonian1b_Potential[i, j]), order)
                    U_t_pn2[i, j] = Taylor1(constant_term(pn2[i, j]) * constant_term(U[i, j]), order)
                    V_t_pn2[i, j] = Taylor1(constant_term(pn2[i, j]) * constant_term(V[i, j]), order)
                    W_t_pn2[i, j] = Taylor1(constant_term(pn2[i, j]) * constant_term(W[i, j]), order)
                    tmp1938[i, j] = Taylor1(constant_term(X[i, j]) * constant_term(newtonianCoeff[i, j]), order)
                    temp_001[i, j] = Taylor1(constant_term(newtonX[j]) + constant_term(tmp1938[i, j]), order)
                    newtonX[j] = Taylor1(identity(constant_term(temp_001[i, j])), order)
                    tmp1940[i, j] = Taylor1(constant_term(Y[i, j]) * constant_term(newtonianCoeff[i, j]), order)
                    temp_002[i, j] = Taylor1(constant_term(newtonY[j]) + constant_term(tmp1940[i, j]), order)
                    newtonY[j] = Taylor1(identity(constant_term(temp_002[i, j])), order)
                    tmp1942[i, j] = Taylor1(constant_term(Z[i, j]) * constant_term(newtonianCoeff[i, j]), order)
                    temp_003[i, j] = Taylor1(constant_term(newtonZ[j]) + constant_term(tmp1942[i, j]), order)
                    newtonZ[j] = Taylor1(identity(constant_term(temp_003[i, j])), order)
                    temp_004[i, j] = Taylor1(constant_term(newtonianNb_Potential[j]) + constant_term(newtonian1b_Potential[i, j]), order)
                    newtonianNb_Potential[j] = Taylor1(identity(constant_term(temp_004[i, j])), order)
                end
            end
            tmp1946[3j - 2] = Taylor1(constant_term(dq[3j - 2]) ^ constant_term(2), order)
            tmp1948[3j - 1] = Taylor1(constant_term(dq[3j - 1]) ^ constant_term(2), order)
            tmp1949[3j - 2] = Taylor1(constant_term(tmp1946[3j - 2]) + constant_term(tmp1948[3j - 1]), order)
            tmp1951[3j] = Taylor1(constant_term(dq[3j]) ^ constant_term(2), order)
            v2[j] = Taylor1(constant_term(tmp1949[3j - 2]) + constant_term(tmp1951[3j]), order)
        end
    tmp1953 = Taylor1(constant_term(ITM_t[1, 1]) + constant_term(ITM_t[2, 2]), order)
    tmp1955 = Taylor1(constant_term(tmp1953) / constant_term(2), order)
    tmp1956 = Taylor1(constant_term(ITM_t[3, 3]) - constant_term(tmp1955), order)
    J2M_t = Taylor1(constant_term(tmp1956) / constant_term(μ[mo]), order)
    tmp1958 = Taylor1(constant_term(ITM_t[2, 2]) - constant_term(ITM_t[1, 1]), order)
    tmp1959 = Taylor1(constant_term(tmp1958) / constant_term(μ[mo]), order)
    C22M_t = Taylor1(constant_term(tmp1959) / constant_term(4), order)
    tmp1962 = Taylor1(-(constant_term(ITM_t[1, 3])), order)
    C21M_t = Taylor1(constant_term(tmp1962) / constant_term(μ[mo]), order)
    tmp1964 = Taylor1(-(constant_term(ITM_t[3, 2])), order)
    S21M_t = Taylor1(constant_term(tmp1964) / constant_term(μ[mo]), order)
    tmp1966 = Taylor1(-(constant_term(ITM_t[2, 1])), order)
    tmp1967 = Taylor1(constant_term(tmp1966) / constant_term(μ[mo]), order)
    S22M_t = Taylor1(constant_term(tmp1967) / constant_term(2), order)
    J2_t[mo] = Taylor1(identity(constant_term(J2M_t)), order)
    tmp1979 = Array{Taylor1{_S}}(undef, size(X_bf_1))
    tmp1979 .= Taylor1(zero(_S), order)
    tmp1981 = Array{Taylor1{_S}}(undef, size(Y_bf_1))
    tmp1981 .= Taylor1(zero(_S), order)
    tmp1983 = Array{Taylor1{_S}}(undef, size(Z_bf_1))
    tmp1983 .= Taylor1(zero(_S), order)
    tmp1987 = Array{Taylor1{_S}}(undef, size(X_bf))
    tmp1987 .= Taylor1(zero(_S), order)
    tmp1989 = Array{Taylor1{_S}}(undef, size(Y_bf))
    tmp1989 .= Taylor1(zero(_S), order)
    tmp1990 = Array{Taylor1{_S}}(undef, size(tmp1987))
    tmp1990 .= Taylor1(zero(_S), order)
    tmp1995 = Array{Taylor1{_S}}(undef, size(P_n))
    tmp1995 .= Taylor1(zero(_S), order)
    tmp1996 = Array{Taylor1{_S}}(undef, size(tmp1995))
    tmp1996 .= Taylor1(zero(_S), order)
    tmp1997 = Array{Taylor1{_S}}(undef, size(P_n))
    tmp1997 .= Taylor1(zero(_S), order)
    tmp1999 = Array{Taylor1{_S}}(undef, size(dP_n))
    tmp1999 .= Taylor1(zero(_S), order)
    tmp2000 = Array{Taylor1{_S}}(undef, size(P_n))
    tmp2000 .= Taylor1(zero(_S), order)
    tmp2005 = Array{Taylor1{_S}}(undef, size(P_n))
    tmp2005 .= Taylor1(zero(_S), order)
    tmp2006 = Array{Taylor1{_S}}(undef, size(tmp2005))
    tmp2006 .= Taylor1(zero(_S), order)
    tmp2008 = Array{Taylor1{_S}}(undef, size(dP_n))
    tmp2008 .= Taylor1(zero(_S), order)
    tmp2009 = Array{Taylor1{_S}}(undef, size(tmp2008))
    tmp2009 .= Taylor1(zero(_S), order)
    tmp2010 = Array{Taylor1{_S}}(undef, size(tmp2009))
    tmp2010 .= Taylor1(zero(_S), order)
    tmp2012 = Array{Taylor1{_S}}(undef, size(P_n))
    tmp2012 .= Taylor1(zero(_S), order)
    tmp2013 = Array{Taylor1{_S}}(undef, size(tmp2012))
    tmp2013 .= Taylor1(zero(_S), order)
    tmp2014 = Array{Taylor1{_S}}(undef, size(tmp2013))
    tmp2014 .= Taylor1(zero(_S), order)
    tmp2016 = Array{Taylor1{_S}}(undef, size(dP_n))
    tmp2016 .= Taylor1(zero(_S), order)
    tmp2017 = Array{Taylor1{_S}}(undef, size(tmp2016))
    tmp2017 .= Taylor1(zero(_S), order)
    tmp2018 = Array{Taylor1{_S}}(undef, size(tmp2017))
    tmp2018 .= Taylor1(zero(_S), order)
    tmp2019 = Array{Taylor1{_S}}(undef, size(tmp2018))
    tmp2019 .= Taylor1(zero(_S), order)
    tmp2021 = Array{Taylor1{_S}}(undef, size(sin_mλ))
    tmp2021 .= Taylor1(zero(_S), order)
    tmp2022 = Array{Taylor1{_S}}(undef, size(cos_mλ))
    tmp2022 .= Taylor1(zero(_S), order)
    tmp2024 = Array{Taylor1{_S}}(undef, size(cos_mλ))
    tmp2024 .= Taylor1(zero(_S), order)
    tmp2025 = Array{Taylor1{_S}}(undef, size(sin_mλ))
    tmp2025 .= Taylor1(zero(_S), order)
    tmp2027 = Array{Taylor1{_S}}(undef, size(secϕ_P_nm))
    tmp2027 .= Taylor1(zero(_S), order)
    tmp2030 = Array{Taylor1{_S}}(undef, size(secϕ_P_nm))
    tmp2030 .= Taylor1(zero(_S), order)
    tmp2032 = Array{Taylor1{_S}}(undef, size(secϕ_P_nm))
    tmp2032 .= Taylor1(zero(_S), order)
    tmp2034 = Array{Taylor1{_S}}(undef, size(secϕ_P_nm))
    tmp2034 .= Taylor1(zero(_S), order)
    tmp2035 = Array{Taylor1{_S}}(undef, size(tmp2034))
    tmp2035 .= Taylor1(zero(_S), order)
    tmp2036 = Array{Taylor1{_S}}(undef, size(secϕ_P_nm))
    tmp2036 .= Taylor1(zero(_S), order)
    tmp2039 = Array{Taylor1{_S}}(undef, size(secϕ_P_nm))
    tmp2039 .= Taylor1(zero(_S), order)
    tmp2040 = Array{Taylor1{_S}}(undef, size(tmp2039))
    tmp2040 .= Taylor1(zero(_S), order)
    tmp2041 = Array{Taylor1{_S}}(undef, size(secϕ_P_nm))
    tmp2041 .= Taylor1(zero(_S), order)
    tmp2043 = Array{Taylor1{_S}}(undef, size(P_nm))
    tmp2043 .= Taylor1(zero(_S), order)
    tmp2044 = Array{Taylor1{_S}}(undef, size(cos_mλ))
    tmp2044 .= Taylor1(zero(_S), order)
    tmp2045 = Array{Taylor1{_S}}(undef, size(sin_mλ))
    tmp2045 .= Taylor1(zero(_S), order)
    tmp2046 = Array{Taylor1{_S}}(undef, size(tmp2044))
    tmp2046 .= Taylor1(zero(_S), order)
    tmp2047 = Array{Taylor1{_S}}(undef, size(tmp2043))
    tmp2047 .= Taylor1(zero(_S), order)
    tmp2048 = Array{Taylor1{_S}}(undef, size(P_nm))
    tmp2048 .= Taylor1(zero(_S), order)
    tmp2049 = Array{Taylor1{_S}}(undef, size(cos_mλ))
    tmp2049 .= Taylor1(zero(_S), order)
    tmp2050 = Array{Taylor1{_S}}(undef, size(sin_mλ))
    tmp2050 .= Taylor1(zero(_S), order)
    tmp2051 = Array{Taylor1{_S}}(undef, size(tmp2049))
    tmp2051 .= Taylor1(zero(_S), order)
    tmp2052 = Array{Taylor1{_S}}(undef, size(tmp2048))
    tmp2052 .= Taylor1(zero(_S), order)
    tmp2053 = Array{Taylor1{_S}}(undef, size(tmp2047))
    tmp2053 .= Taylor1(zero(_S), order)
    tmp2055 = Array{Taylor1{_S}}(undef, size(secϕ_P_nm))
    tmp2055 .= Taylor1(zero(_S), order)
    tmp2056 = Array{Taylor1{_S}}(undef, size(cos_mλ))
    tmp2056 .= Taylor1(zero(_S), order)
    tmp2057 = Array{Taylor1{_S}}(undef, size(sin_mλ))
    tmp2057 .= Taylor1(zero(_S), order)
    tmp2058 = Array{Taylor1{_S}}(undef, size(tmp2056))
    tmp2058 .= Taylor1(zero(_S), order)
    tmp2059 = Array{Taylor1{_S}}(undef, size(tmp2055))
    tmp2059 .= Taylor1(zero(_S), order)
    tmp2060 = Array{Taylor1{_S}}(undef, size(secϕ_P_nm))
    tmp2060 .= Taylor1(zero(_S), order)
    tmp2061 = Array{Taylor1{_S}}(undef, size(cos_mλ))
    tmp2061 .= Taylor1(zero(_S), order)
    tmp2062 = Array{Taylor1{_S}}(undef, size(sin_mλ))
    tmp2062 .= Taylor1(zero(_S), order)
    tmp2063 = Array{Taylor1{_S}}(undef, size(tmp2061))
    tmp2063 .= Taylor1(zero(_S), order)
    tmp2064 = Array{Taylor1{_S}}(undef, size(tmp2060))
    tmp2064 .= Taylor1(zero(_S), order)
    tmp2065 = Array{Taylor1{_S}}(undef, size(tmp2059))
    tmp2065 .= Taylor1(zero(_S), order)
    tmp2067 = Array{Taylor1{_S}}(undef, size(cos_mλ))
    tmp2067 .= Taylor1(zero(_S), order)
    tmp2068 = Array{Taylor1{_S}}(undef, size(sin_mλ))
    tmp2068 .= Taylor1(zero(_S), order)
    tmp2069 = Array{Taylor1{_S}}(undef, size(tmp2067))
    tmp2069 .= Taylor1(zero(_S), order)
    tmp2070 = Array{Taylor1{_S}}(undef, size(cosϕ_dP_nm))
    tmp2070 .= Taylor1(zero(_S), order)
    tmp2071 = Array{Taylor1{_S}}(undef, size(cos_mλ))
    tmp2071 .= Taylor1(zero(_S), order)
    tmp2072 = Array{Taylor1{_S}}(undef, size(sin_mλ))
    tmp2072 .= Taylor1(zero(_S), order)
    tmp2073 = Array{Taylor1{_S}}(undef, size(tmp2071))
    tmp2073 .= Taylor1(zero(_S), order)
    tmp2074 = Array{Taylor1{_S}}(undef, size(cosϕ_dP_nm))
    tmp2074 .= Taylor1(zero(_S), order)
    tmp2075 = Array{Taylor1{_S}}(undef, size(tmp2070))
    tmp2075 .= Taylor1(zero(_S), order)
    tmp2077 = Array{Taylor1{_S}}(undef, size(P_nm))
    tmp2077 .= Taylor1(zero(_S), order)
    tmp2078 = Array{Taylor1{_S}}(undef, size(cos_mλ))
    tmp2078 .= Taylor1(zero(_S), order)
    tmp2079 = Array{Taylor1{_S}}(undef, size(sin_mλ))
    tmp2079 .= Taylor1(zero(_S), order)
    tmp2080 = Array{Taylor1{_S}}(undef, size(tmp2078))
    tmp2080 .= Taylor1(zero(_S), order)
    tmp2081 = Array{Taylor1{_S}}(undef, size(tmp2077))
    tmp2081 .= Taylor1(zero(_S), order)
    tmp2082 = Array{Taylor1{_S}}(undef, size(tmp2081))
    tmp2082 .= Taylor1(zero(_S), order)
    temp_CS_ξ = Array{Taylor1{_S}}(undef, size(tmp2082))
    temp_CS_ξ .= Taylor1(zero(_S), order)
    tmp2084 = Array{Taylor1{_S}}(undef, size(secϕ_P_nm))
    tmp2084 .= Taylor1(zero(_S), order)
    tmp2085 = Array{Taylor1{_S}}(undef, size(cos_mλ))
    tmp2085 .= Taylor1(zero(_S), order)
    tmp2086 = Array{Taylor1{_S}}(undef, size(sin_mλ))
    tmp2086 .= Taylor1(zero(_S), order)
    tmp2087 = Array{Taylor1{_S}}(undef, size(tmp2085))
    tmp2087 .= Taylor1(zero(_S), order)
    tmp2088 = Array{Taylor1{_S}}(undef, size(tmp2084))
    tmp2088 .= Taylor1(zero(_S), order)
    tmp2089 = Array{Taylor1{_S}}(undef, size(tmp2088))
    tmp2089 .= Taylor1(zero(_S), order)
    temp_CS_η = Array{Taylor1{_S}}(undef, size(tmp2089))
    temp_CS_η .= Taylor1(zero(_S), order)
    tmp2091 = Array{Taylor1{_S}}(undef, size(cos_mλ))
    tmp2091 .= Taylor1(zero(_S), order)
    tmp2092 = Array{Taylor1{_S}}(undef, size(sin_mλ))
    tmp2092 .= Taylor1(zero(_S), order)
    tmp2093 = Array{Taylor1{_S}}(undef, size(tmp2091))
    tmp2093 .= Taylor1(zero(_S), order)
    tmp2094 = Array{Taylor1{_S}}(undef, size(cosϕ_dP_nm))
    tmp2094 .= Taylor1(zero(_S), order)
    tmp2095 = Array{Taylor1{_S}}(undef, size(tmp2094))
    tmp2095 .= Taylor1(zero(_S), order)
    temp_CS_ζ = Array{Taylor1{_S}}(undef, size(tmp2095))
    temp_CS_ζ .= Taylor1(zero(_S), order)
    tmp2097 = Array{Taylor1{_S}}(undef, size(F_J_ξ))
    tmp2097 .= Taylor1(zero(_S), order)
    tmp2098 = Array{Taylor1{_S}}(undef, size(F_CS_ξ))
    tmp2098 .= Taylor1(zero(_S), order)
    tmp2101 = Array{Taylor1{_S}}(undef, size(F_J_ζ))
    tmp2101 .= Taylor1(zero(_S), order)
    tmp2102 = Array{Taylor1{_S}}(undef, size(F_CS_ζ))
    tmp2102 .= Taylor1(zero(_S), order)
    tmp2108 = Array{Taylor1{_S}}(undef, size(sin_ϕ))
    tmp2108 .= Taylor1(zero(_S), order)
    tmp2111 = Array{Taylor1{_S}}(undef, size(sin_ϕ))
    tmp2111 .= Taylor1(zero(_S), order)
    tmp2113 = Array{Taylor1{_S}}(undef, size(Rb2p))
    tmp2113 .= Taylor1(zero(_S), order)
    tmp2114 = Array{Taylor1{_S}}(undef, size(Rb2p))
    tmp2114 .= Taylor1(zero(_S), order)
    tmp2115 = Array{Taylor1{_S}}(undef, size(tmp2113))
    tmp2115 .= Taylor1(zero(_S), order)
    tmp2116 = Array{Taylor1{_S}}(undef, size(Rb2p))
    tmp2116 .= Taylor1(zero(_S), order)
    tmp2118 = Array{Taylor1{_S}}(undef, size(Rb2p))
    tmp2118 .= Taylor1(zero(_S), order)
    tmp2119 = Array{Taylor1{_S}}(undef, size(Rb2p))
    tmp2119 .= Taylor1(zero(_S), order)
    tmp2120 = Array{Taylor1{_S}}(undef, size(tmp2118))
    tmp2120 .= Taylor1(zero(_S), order)
    tmp2121 = Array{Taylor1{_S}}(undef, size(Rb2p))
    tmp2121 .= Taylor1(zero(_S), order)
    tmp2123 = Array{Taylor1{_S}}(undef, size(Rb2p))
    tmp2123 .= Taylor1(zero(_S), order)
    tmp2124 = Array{Taylor1{_S}}(undef, size(Rb2p))
    tmp2124 .= Taylor1(zero(_S), order)
    tmp2125 = Array{Taylor1{_S}}(undef, size(tmp2123))
    tmp2125 .= Taylor1(zero(_S), order)
    tmp2126 = Array{Taylor1{_S}}(undef, size(Rb2p))
    tmp2126 .= Taylor1(zero(_S), order)
    tmp2128 = Array{Taylor1{_S}}(undef, size(Rb2p))
    tmp2128 .= Taylor1(zero(_S), order)
    tmp2129 = Array{Taylor1{_S}}(undef, size(Rb2p))
    tmp2129 .= Taylor1(zero(_S), order)
    tmp2130 = Array{Taylor1{_S}}(undef, size(tmp2128))
    tmp2130 .= Taylor1(zero(_S), order)
    tmp2131 = Array{Taylor1{_S}}(undef, size(Rb2p))
    tmp2131 .= Taylor1(zero(_S), order)
    tmp2133 = Array{Taylor1{_S}}(undef, size(Rb2p))
    tmp2133 .= Taylor1(zero(_S), order)
    tmp2134 = Array{Taylor1{_S}}(undef, size(Rb2p))
    tmp2134 .= Taylor1(zero(_S), order)
    tmp2135 = Array{Taylor1{_S}}(undef, size(tmp2133))
    tmp2135 .= Taylor1(zero(_S), order)
    tmp2136 = Array{Taylor1{_S}}(undef, size(Rb2p))
    tmp2136 .= Taylor1(zero(_S), order)
    tmp2138 = Array{Taylor1{_S}}(undef, size(Rb2p))
    tmp2138 .= Taylor1(zero(_S), order)
    tmp2139 = Array{Taylor1{_S}}(undef, size(Rb2p))
    tmp2139 .= Taylor1(zero(_S), order)
    tmp2140 = Array{Taylor1{_S}}(undef, size(tmp2138))
    tmp2140 .= Taylor1(zero(_S), order)
    tmp2141 = Array{Taylor1{_S}}(undef, size(Rb2p))
    tmp2141 .= Taylor1(zero(_S), order)
    tmp2143 = Array{Taylor1{_S}}(undef, size(Rb2p))
    tmp2143 .= Taylor1(zero(_S), order)
    tmp2144 = Array{Taylor1{_S}}(undef, size(Rb2p))
    tmp2144 .= Taylor1(zero(_S), order)
    tmp2145 = Array{Taylor1{_S}}(undef, size(tmp2143))
    tmp2145 .= Taylor1(zero(_S), order)
    tmp2146 = Array{Taylor1{_S}}(undef, size(Rb2p))
    tmp2146 .= Taylor1(zero(_S), order)
    tmp2148 = Array{Taylor1{_S}}(undef, size(Rb2p))
    tmp2148 .= Taylor1(zero(_S), order)
    tmp2149 = Array{Taylor1{_S}}(undef, size(Rb2p))
    tmp2149 .= Taylor1(zero(_S), order)
    tmp2150 = Array{Taylor1{_S}}(undef, size(tmp2148))
    tmp2150 .= Taylor1(zero(_S), order)
    tmp2151 = Array{Taylor1{_S}}(undef, size(Rb2p))
    tmp2151 .= Taylor1(zero(_S), order)
    tmp2153 = Array{Taylor1{_S}}(undef, size(Rb2p))
    tmp2153 .= Taylor1(zero(_S), order)
    tmp2154 = Array{Taylor1{_S}}(undef, size(Rb2p))
    tmp2154 .= Taylor1(zero(_S), order)
    tmp2155 = Array{Taylor1{_S}}(undef, size(tmp2153))
    tmp2155 .= Taylor1(zero(_S), order)
    tmp2156 = Array{Taylor1{_S}}(undef, size(Rb2p))
    tmp2156 .= Taylor1(zero(_S), order)
    tmp2158 = Array{Taylor1{_S}}(undef, size(Gc2p))
    tmp2158 .= Taylor1(zero(_S), order)
    tmp2159 = Array{Taylor1{_S}}(undef, size(Gc2p))
    tmp2159 .= Taylor1(zero(_S), order)
    tmp2160 = Array{Taylor1{_S}}(undef, size(tmp2158))
    tmp2160 .= Taylor1(zero(_S), order)
    tmp2161 = Array{Taylor1{_S}}(undef, size(Gc2p))
    tmp2161 .= Taylor1(zero(_S), order)
    tmp2163 = Array{Taylor1{_S}}(undef, size(Gc2p))
    tmp2163 .= Taylor1(zero(_S), order)
    tmp2164 = Array{Taylor1{_S}}(undef, size(Gc2p))
    tmp2164 .= Taylor1(zero(_S), order)
    tmp2165 = Array{Taylor1{_S}}(undef, size(tmp2163))
    tmp2165 .= Taylor1(zero(_S), order)
    tmp2166 = Array{Taylor1{_S}}(undef, size(Gc2p))
    tmp2166 .= Taylor1(zero(_S), order)
    tmp2168 = Array{Taylor1{_S}}(undef, size(Gc2p))
    tmp2168 .= Taylor1(zero(_S), order)
    tmp2169 = Array{Taylor1{_S}}(undef, size(Gc2p))
    tmp2169 .= Taylor1(zero(_S), order)
    tmp2170 = Array{Taylor1{_S}}(undef, size(tmp2168))
    tmp2170 .= Taylor1(zero(_S), order)
    tmp2171 = Array{Taylor1{_S}}(undef, size(Gc2p))
    tmp2171 .= Taylor1(zero(_S), order)
    #= REPL[4]:265 =# Threads.@threads for j = 1:N_ext
            for i = 1:N_ext
                if i == j
                    continue
                else
                    if UJ_interaction[i, j]
                        X_bf_1[i, j] = Taylor1(constant_term(X[i, j]) * constant_term(M_[1, 1, j]), order)
                        X_bf_2[i, j] = Taylor1(constant_term(Y[i, j]) * constant_term(M_[1, 2, j]), order)
                        X_bf_3[i, j] = Taylor1(constant_term(Z[i, j]) * constant_term(M_[1, 3, j]), order)
                        Y_bf_1[i, j] = Taylor1(constant_term(X[i, j]) * constant_term(M_[2, 1, j]), order)
                        Y_bf_2[i, j] = Taylor1(constant_term(Y[i, j]) * constant_term(M_[2, 2, j]), order)
                        Y_bf_3[i, j] = Taylor1(constant_term(Z[i, j]) * constant_term(M_[2, 3, j]), order)
                        Z_bf_1[i, j] = Taylor1(constant_term(X[i, j]) * constant_term(M_[3, 1, j]), order)
                        Z_bf_2[i, j] = Taylor1(constant_term(Y[i, j]) * constant_term(M_[3, 2, j]), order)
                        Z_bf_3[i, j] = Taylor1(constant_term(Z[i, j]) * constant_term(M_[3, 3, j]), order)
                        tmp1979[i, j] = Taylor1(constant_term(X_bf_1[i, j]) + constant_term(X_bf_2[i, j]), order)
                        X_bf[i, j] = Taylor1(constant_term(tmp1979[i, j]) + constant_term(X_bf_3[i, j]), order)
                        tmp1981[i, j] = Taylor1(constant_term(Y_bf_1[i, j]) + constant_term(Y_bf_2[i, j]), order)
                        Y_bf[i, j] = Taylor1(constant_term(tmp1981[i, j]) + constant_term(Y_bf_3[i, j]), order)
                        tmp1983[i, j] = Taylor1(constant_term(Z_bf_1[i, j]) + constant_term(Z_bf_2[i, j]), order)
                        Z_bf[i, j] = Taylor1(constant_term(tmp1983[i, j]) + constant_term(Z_bf_3[i, j]), order)
                        sin_ϕ[i, j] = Taylor1(constant_term(Z_bf[i, j]) / constant_term(r_p1d2[i, j]), order)
                        tmp1987[i, j] = Taylor1(constant_term(X_bf[i, j]) ^ constant_term(2), order)
                        tmp1989[i, j] = Taylor1(constant_term(Y_bf[i, j]) ^ constant_term(2), order)
                        tmp1990[i, j] = Taylor1(constant_term(tmp1987[i, j]) + constant_term(tmp1989[i, j]), order)
                        r_xy[i, j] = Taylor1(sqrt(constant_term(tmp1990[i, j])), order)
                        cos_ϕ[i, j] = Taylor1(constant_term(r_xy[i, j]) / constant_term(r_p1d2[i, j]), order)
                        sin_λ[i, j] = Taylor1(constant_term(Y_bf[i, j]) / constant_term(r_xy[i, j]), order)
                        cos_λ[i, j] = Taylor1(constant_term(X_bf[i, j]) / constant_term(r_xy[i, j]), order)
                        P_n[i, j, 1] = Taylor1(identity(constant_term(one_t)), order)
                        P_n[i, j, 2] = Taylor1(identity(constant_term(sin_ϕ[i, j])), order)
                        dP_n[i, j, 1] = Taylor1(identity(constant_term(zero_q_1)), order)
                        dP_n[i, j, 2] = Taylor1(identity(constant_term(one_t)), order)
                        for n = 2:n1SEM[j]
                            tmp1995[i, j, n] = Taylor1(constant_term(P_n[i, j, n]) * constant_term(sin_ϕ[i, j]), order)
                            tmp1996[i, j, n] = Taylor1(constant_term(tmp1995[i, j, n]) * constant_term(fact1_jsem[n]), order)
                            tmp1997[i, j, n - 1] = Taylor1(constant_term(P_n[i, j, n - 1]) * constant_term(fact2_jsem[n]), order)
                            P_n[i, j, n + 1] = Taylor1(constant_term(tmp1996[i, j, n]) - constant_term(tmp1997[i, j, n - 1]), order)
                            tmp1999[i, j, n] = Taylor1(constant_term(dP_n[i, j, n]) * constant_term(sin_ϕ[i, j]), order)
                            tmp2000[i, j, n] = Taylor1(constant_term(P_n[i, j, n]) * constant_term(fact3_jsem[n]), order)
                            dP_n[i, j, n + 1] = Taylor1(constant_term(tmp1999[i, j, n]) + constant_term(tmp2000[i, j, n]), order)
                            temp_rn[i, j, n] = Taylor1(constant_term(r_p1d2[i, j]) ^ constant_term(fact5_jsem[n]), order)
                        end
                        r_p4[i, j] = Taylor1(constant_term(r_p2[i, j]) ^ constant_term(2), order)
                        tmp2005[i, j, 3] = Taylor1(constant_term(P_n[i, j, 3]) * constant_term(fact4_jsem[2]), order)
                        tmp2006[i, j, 3] = Taylor1(constant_term(tmp2005[i, j, 3]) * constant_term(J2_t[j]), order)
                        F_J_ξ[i, j] = Taylor1(constant_term(tmp2006[i, j, 3]) / constant_term(r_p4[i, j]), order)
                        tmp2008[i, j, 3] = Taylor1(-(constant_term(dP_n[i, j, 3])), order)
                        tmp2009[i, j, 3] = Taylor1(constant_term(tmp2008[i, j, 3]) * constant_term(cos_ϕ[i, j]), order)
                        tmp2010[i, j, 3] = Taylor1(constant_term(tmp2009[i, j, 3]) * constant_term(J2_t[j]), order)
                        F_J_ζ[i, j] = Taylor1(constant_term(tmp2010[i, j, 3]) / constant_term(r_p4[i, j]), order)
                        F_J_ξ_36[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                        F_J_ζ_36[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                        for n = 3:n1SEM[j]
                            tmp2012[i, j, n + 1] = Taylor1(constant_term(P_n[i, j, n + 1]) * constant_term(fact4_jsem[n]), order)
                            tmp2013[i, j, n + 1] = Taylor1(constant_term(tmp2012[i, j, n + 1]) * constant_term(JSEM[j, n]), order)
                            tmp2014[i, j, n + 1] = Taylor1(constant_term(tmp2013[i, j, n + 1]) / constant_term(temp_rn[i, j, n]), order)
                            temp_fjξ[i, j, n] = Taylor1(constant_term(tmp2014[i, j, n + 1]) + constant_term(F_J_ξ_36[i, j]), order)
                            tmp2016[i, j, n + 1] = Taylor1(-(constant_term(dP_n[i, j, n + 1])), order)
                            tmp2017[i, j, n + 1] = Taylor1(constant_term(tmp2016[i, j, n + 1]) * constant_term(cos_ϕ[i, j]), order)
                            tmp2018[i, j, n + 1] = Taylor1(constant_term(tmp2017[i, j, n + 1]) * constant_term(JSEM[j, n]), order)
                            tmp2019[i, j, n + 1] = Taylor1(constant_term(tmp2018[i, j, n + 1]) / constant_term(temp_rn[i, j, n]), order)
                            temp_fjζ[i, j, n] = Taylor1(constant_term(tmp2019[i, j, n + 1]) + constant_term(F_J_ζ_36[i, j]), order)
                            F_J_ξ_36[i, j] = Taylor1(identity(constant_term(temp_fjξ[i, j, n])), order)
                            F_J_ζ_36[i, j] = Taylor1(identity(constant_term(temp_fjζ[i, j, n])), order)
                        end
                        if j == mo
                            for m = 1:n1SEM[mo]
                                if m == 1
                                    sin_mλ[i, j, 1] = Taylor1(identity(constant_term(sin_λ[i, j])), order)
                                    cos_mλ[i, j, 1] = Taylor1(identity(constant_term(cos_λ[i, j])), order)
                                    secϕ_P_nm[i, j, 1, 1] = Taylor1(identity(constant_term(one_t)), order)
                                else
                                    tmp2021[i, j, 1] = Taylor1(constant_term(sin_mλ[i, j, 1]) * constant_term(cos_mλ[i, j, m - 1]), order)
                                    tmp2022[i, j, 1] = Taylor1(constant_term(cos_mλ[i, j, 1]) * constant_term(sin_mλ[i, j, m - 1]), order)
                                    sin_mλ[i, j, m] = Taylor1(constant_term(tmp2021[i, j, 1]) + constant_term(tmp2022[i, j, 1]), order)
                                    tmp2024[i, j, 1] = Taylor1(constant_term(cos_mλ[i, j, 1]) * constant_term(cos_mλ[i, j, m - 1]), order)
                                    tmp2025[i, j, 1] = Taylor1(constant_term(sin_mλ[i, j, 1]) * constant_term(sin_mλ[i, j, m - 1]), order)
                                    cos_mλ[i, j, m] = Taylor1(constant_term(tmp2024[i, j, 1]) - constant_term(tmp2025[i, j, 1]), order)
                                    tmp2027[i, j, m - 1, m - 1] = Taylor1(constant_term(secϕ_P_nm[i, j, m - 1, m - 1]) * constant_term(cos_ϕ[i, j]), order)
                                    secϕ_P_nm[i, j, m, m] = Taylor1(constant_term(tmp2027[i, j, m - 1, m - 1]) * constant_term(lnm5[m]), order)
                                    P_nm[i, j, m, m] = Taylor1(constant_term(secϕ_P_nm[i, j, m, m]) * constant_term(cos_ϕ[i, j]), order)
                                    tmp2030[i, j, m, m] = Taylor1(constant_term(secϕ_P_nm[i, j, m, m]) * constant_term(sin_ϕ[i, j]), order)
                                    cosϕ_dP_nm[i, j, m, m] = Taylor1(constant_term(tmp2030[i, j, m, m]) * constant_term(lnm3[m]), order)
                                end
                                for n = m + 1:n1SEM[mo]
                                    if n == m + 1
                                        tmp2032[i, j, n - 1, m] = Taylor1(constant_term(secϕ_P_nm[i, j, n - 1, m]) * constant_term(sin_ϕ[i, j]), order)
                                        secϕ_P_nm[i, j, n, m] = Taylor1(constant_term(tmp2032[i, j, n - 1, m]) * constant_term(lnm1[n, m]), order)
                                    else
                                        tmp2034[i, j, n - 1, m] = Taylor1(constant_term(secϕ_P_nm[i, j, n - 1, m]) * constant_term(sin_ϕ[i, j]), order)
                                        tmp2035[i, j, n - 1, m] = Taylor1(constant_term(tmp2034[i, j, n - 1, m]) * constant_term(lnm1[n, m]), order)
                                        tmp2036[i, j, n - 2, m] = Taylor1(constant_term(secϕ_P_nm[i, j, n - 2, m]) * constant_term(lnm2[n, m]), order)
                                        secϕ_P_nm[i, j, n, m] = Taylor1(constant_term(tmp2035[i, j, n - 1, m]) + constant_term(tmp2036[i, j, n - 2, m]), order)
                                    end
                                    P_nm[i, j, n, m] = Taylor1(constant_term(secϕ_P_nm[i, j, n, m]) * constant_term(cos_ϕ[i, j]), order)
                                    tmp2039[i, j, n, m] = Taylor1(constant_term(secϕ_P_nm[i, j, n, m]) * constant_term(sin_ϕ[i, j]), order)
                                    tmp2040[i, j, n, m] = Taylor1(constant_term(tmp2039[i, j, n, m]) * constant_term(lnm3[n]), order)
                                    tmp2041[i, j, n - 1, m] = Taylor1(constant_term(secϕ_P_nm[i, j, n - 1, m]) * constant_term(lnm4[n, m]), order)
                                    cosϕ_dP_nm[i, j, n, m] = Taylor1(constant_term(tmp2040[i, j, n, m]) + constant_term(tmp2041[i, j, n - 1, m]), order)
                                end
                            end
                            tmp2043[i, j, 2, 1] = Taylor1(constant_term(P_nm[i, j, 2, 1]) * constant_term(lnm6[2]), order)
                            tmp2044[i, j, 1] = Taylor1(constant_term(C21M_t) * constant_term(cos_mλ[i, j, 1]), order)
                            tmp2045[i, j, 1] = Taylor1(constant_term(S21M_t) * constant_term(sin_mλ[i, j, 1]), order)
                            tmp2046[i, j, 1] = Taylor1(constant_term(tmp2044[i, j, 1]) + constant_term(tmp2045[i, j, 1]), order)
                            tmp2047[i, j, 2, 1] = Taylor1(constant_term(tmp2043[i, j, 2, 1]) * constant_term(tmp2046[i, j, 1]), order)
                            tmp2048[i, j, 2, 2] = Taylor1(constant_term(P_nm[i, j, 2, 2]) * constant_term(lnm6[2]), order)
                            tmp2049[i, j, 2] = Taylor1(constant_term(C22M_t) * constant_term(cos_mλ[i, j, 2]), order)
                            tmp2050[i, j, 2] = Taylor1(constant_term(S22M_t) * constant_term(sin_mλ[i, j, 2]), order)
                            tmp2051[i, j, 2] = Taylor1(constant_term(tmp2049[i, j, 2]) + constant_term(tmp2050[i, j, 2]), order)
                            tmp2052[i, j, 2, 2] = Taylor1(constant_term(tmp2048[i, j, 2, 2]) * constant_term(tmp2051[i, j, 2]), order)
                            tmp2053[i, j, 2, 1] = Taylor1(constant_term(tmp2047[i, j, 2, 1]) + constant_term(tmp2052[i, j, 2, 2]), order)
                            F_CS_ξ[i, j] = Taylor1(constant_term(tmp2053[i, j, 2, 1]) / constant_term(r_p4[i, j]), order)
                            tmp2055[i, j, 2, 1] = Taylor1(constant_term(secϕ_P_nm[i, j, 2, 1]) * constant_term(lnm7[1]), order)
                            tmp2056[i, j, 1] = Taylor1(constant_term(S21M_t) * constant_term(cos_mλ[i, j, 1]), order)
                            tmp2057[i, j, 1] = Taylor1(constant_term(C21M_t) * constant_term(sin_mλ[i, j, 1]), order)
                            tmp2058[i, j, 1] = Taylor1(constant_term(tmp2056[i, j, 1]) - constant_term(tmp2057[i, j, 1]), order)
                            tmp2059[i, j, 2, 1] = Taylor1(constant_term(tmp2055[i, j, 2, 1]) * constant_term(tmp2058[i, j, 1]), order)
                            tmp2060[i, j, 2, 2] = Taylor1(constant_term(secϕ_P_nm[i, j, 2, 2]) * constant_term(lnm7[2]), order)
                            tmp2061[i, j, 2] = Taylor1(constant_term(S22M_t) * constant_term(cos_mλ[i, j, 2]), order)
                            tmp2062[i, j, 2] = Taylor1(constant_term(C22M_t) * constant_term(sin_mλ[i, j, 2]), order)
                            tmp2063[i, j, 2] = Taylor1(constant_term(tmp2061[i, j, 2]) - constant_term(tmp2062[i, j, 2]), order)
                            tmp2064[i, j, 2, 2] = Taylor1(constant_term(tmp2060[i, j, 2, 2]) * constant_term(tmp2063[i, j, 2]), order)
                            tmp2065[i, j, 2, 1] = Taylor1(constant_term(tmp2059[i, j, 2, 1]) + constant_term(tmp2064[i, j, 2, 2]), order)
                            F_CS_η[i, j] = Taylor1(constant_term(tmp2065[i, j, 2, 1]) / constant_term(r_p4[i, j]), order)
                            tmp2067[i, j, 1] = Taylor1(constant_term(C21M_t) * constant_term(cos_mλ[i, j, 1]), order)
                            tmp2068[i, j, 1] = Taylor1(constant_term(S21M_t) * constant_term(sin_mλ[i, j, 1]), order)
                            tmp2069[i, j, 1] = Taylor1(constant_term(tmp2067[i, j, 1]) + constant_term(tmp2068[i, j, 1]), order)
                            tmp2070[i, j, 2, 1] = Taylor1(constant_term(cosϕ_dP_nm[i, j, 2, 1]) * constant_term(tmp2069[i, j, 1]), order)
                            tmp2071[i, j, 2] = Taylor1(constant_term(C22M_t) * constant_term(cos_mλ[i, j, 2]), order)
                            tmp2072[i, j, 2] = Taylor1(constant_term(S22M_t) * constant_term(sin_mλ[i, j, 2]), order)
                            tmp2073[i, j, 2] = Taylor1(constant_term(tmp2071[i, j, 2]) + constant_term(tmp2072[i, j, 2]), order)
                            tmp2074[i, j, 2, 2] = Taylor1(constant_term(cosϕ_dP_nm[i, j, 2, 2]) * constant_term(tmp2073[i, j, 2]), order)
                            tmp2075[i, j, 2, 1] = Taylor1(constant_term(tmp2070[i, j, 2, 1]) + constant_term(tmp2074[i, j, 2, 2]), order)
                            F_CS_ζ[i, j] = Taylor1(constant_term(tmp2075[i, j, 2, 1]) / constant_term(r_p4[i, j]), order)
                            F_CS_ξ_36[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                            F_CS_η_36[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                            F_CS_ζ_36[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                            for n = 3:n1SEM[mo]
                                for m = 1:n
                                    tmp2077[i, j, n, m] = Taylor1(constant_term(P_nm[i, j, n, m]) * constant_term(lnm6[n]), order)
                                    tmp2078[i, j, m] = Taylor1(constant_term(cos_mλ[i, j, m]) * constant_term(CM[n, m]), order)
                                    tmp2079[i, j, m] = Taylor1(constant_term(sin_mλ[i, j, m]) * constant_term(SM[n, m]), order)
                                    tmp2080[i, j, m] = Taylor1(constant_term(tmp2078[i, j, m]) + constant_term(tmp2079[i, j, m]), order)
                                    tmp2081[i, j, n, m] = Taylor1(constant_term(tmp2077[i, j, n, m]) * constant_term(tmp2080[i, j, m]), order)
                                    tmp2082[i, j, n, m] = Taylor1(constant_term(tmp2081[i, j, n, m]) / constant_term(temp_rn[i, j, n]), order)
                                    temp_CS_ξ[i, j, n, m] = Taylor1(constant_term(tmp2082[i, j, n, m]) + constant_term(F_CS_ξ_36[i, j]), order)
                                    tmp2084[i, j, n, m] = Taylor1(constant_term(secϕ_P_nm[i, j, n, m]) * constant_term(lnm7[m]), order)
                                    tmp2085[i, j, m] = Taylor1(constant_term(cos_mλ[i, j, m]) * constant_term(SM[n, m]), order)
                                    tmp2086[i, j, m] = Taylor1(constant_term(sin_mλ[i, j, m]) * constant_term(CM[n, m]), order)
                                    tmp2087[i, j, m] = Taylor1(constant_term(tmp2085[i, j, m]) - constant_term(tmp2086[i, j, m]), order)
                                    tmp2088[i, j, n, m] = Taylor1(constant_term(tmp2084[i, j, n, m]) * constant_term(tmp2087[i, j, m]), order)
                                    tmp2089[i, j, n, m] = Taylor1(constant_term(tmp2088[i, j, n, m]) / constant_term(temp_rn[i, j, n]), order)
                                    temp_CS_η[i, j, n, m] = Taylor1(constant_term(tmp2089[i, j, n, m]) + constant_term(F_CS_η_36[i, j]), order)
                                    tmp2091[i, j, m] = Taylor1(constant_term(cos_mλ[i, j, m]) * constant_term(CM[n, m]), order)
                                    tmp2092[i, j, m] = Taylor1(constant_term(sin_mλ[i, j, m]) * constant_term(SM[n, m]), order)
                                    tmp2093[i, j, m] = Taylor1(constant_term(tmp2091[i, j, m]) + constant_term(tmp2092[i, j, m]), order)
                                    tmp2094[i, j, n, m] = Taylor1(constant_term(cosϕ_dP_nm[i, j, n, m]) * constant_term(tmp2093[i, j, m]), order)
                                    tmp2095[i, j, n, m] = Taylor1(constant_term(tmp2094[i, j, n, m]) / constant_term(temp_rn[i, j, n]), order)
                                    temp_CS_ζ[i, j, n, m] = Taylor1(constant_term(tmp2095[i, j, n, m]) + constant_term(F_CS_ζ_36[i, j]), order)
                                    F_CS_ξ_36[i, j] = Taylor1(identity(constant_term(temp_CS_ξ[i, j, n, m])), order)
                                    F_CS_η_36[i, j] = Taylor1(identity(constant_term(temp_CS_η[i, j, n, m])), order)
                                    F_CS_ζ_36[i, j] = Taylor1(identity(constant_term(temp_CS_ζ[i, j, n, m])), order)
                                end
                            end
                            tmp2097[i, j] = Taylor1(constant_term(F_J_ξ[i, j]) + constant_term(F_J_ξ_36[i, j]), order)
                            tmp2098[i, j] = Taylor1(constant_term(F_CS_ξ[i, j]) + constant_term(F_CS_ξ_36[i, j]), order)
                            F_JCS_ξ[i, j] = Taylor1(constant_term(tmp2097[i, j]) + constant_term(tmp2098[i, j]), order)
                            F_JCS_η[i, j] = Taylor1(constant_term(F_CS_η[i, j]) + constant_term(F_CS_η_36[i, j]), order)
                            tmp2101[i, j] = Taylor1(constant_term(F_J_ζ[i, j]) + constant_term(F_J_ζ_36[i, j]), order)
                            tmp2102[i, j] = Taylor1(constant_term(F_CS_ζ[i, j]) + constant_term(F_CS_ζ_36[i, j]), order)
                            F_JCS_ζ[i, j] = Taylor1(constant_term(tmp2101[i, j]) + constant_term(tmp2102[i, j]), order)
                        else
                            F_JCS_ξ[i, j] = Taylor1(constant_term(F_J_ξ[i, j]) + constant_term(F_J_ξ_36[i, j]), order)
                            F_JCS_η[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                            F_JCS_ζ[i, j] = Taylor1(constant_term(F_J_ζ[i, j]) + constant_term(F_J_ζ_36[i, j]), order)
                        end
                        Rb2p[i, j, 1, 1] = Taylor1(constant_term(cos_ϕ[i, j]) * constant_term(cos_λ[i, j]), order)
                        Rb2p[i, j, 2, 1] = Taylor1(-(constant_term(sin_λ[i, j])), order)
                        tmp2108[i, j] = Taylor1(-(constant_term(sin_ϕ[i, j])), order)
                        Rb2p[i, j, 3, 1] = Taylor1(constant_term(tmp2108[i, j]) * constant_term(cos_λ[i, j]), order)
                        Rb2p[i, j, 1, 2] = Taylor1(constant_term(cos_ϕ[i, j]) * constant_term(sin_λ[i, j]), order)
                        Rb2p[i, j, 2, 2] = Taylor1(identity(constant_term(cos_λ[i, j])), order)
                        tmp2111[i, j] = Taylor1(-(constant_term(sin_ϕ[i, j])), order)
                        Rb2p[i, j, 3, 2] = Taylor1(constant_term(tmp2111[i, j]) * constant_term(sin_λ[i, j]), order)
                        Rb2p[i, j, 1, 3] = Taylor1(identity(constant_term(sin_ϕ[i, j])), order)
                        Rb2p[i, j, 2, 3] = Taylor1(identity(constant_term(zero_q_1)), order)
                        Rb2p[i, j, 3, 3] = Taylor1(identity(constant_term(cos_ϕ[i, j])), order)
                        tmp2113[i, j, 1, 1] = Taylor1(constant_term(Rb2p[i, j, 1, 1]) * constant_term(M_[1, 1, j]), order)
                        tmp2114[i, j, 1, 2] = Taylor1(constant_term(Rb2p[i, j, 1, 2]) * constant_term(M_[2, 1, j]), order)
                        tmp2115[i, j, 1, 1] = Taylor1(constant_term(tmp2113[i, j, 1, 1]) + constant_term(tmp2114[i, j, 1, 2]), order)
                        tmp2116[i, j, 1, 3] = Taylor1(constant_term(Rb2p[i, j, 1, 3]) * constant_term(M_[3, 1, j]), order)
                        Gc2p[i, j, 1, 1] = Taylor1(constant_term(tmp2115[i, j, 1, 1]) + constant_term(tmp2116[i, j, 1, 3]), order)
                        tmp2118[i, j, 2, 1] = Taylor1(constant_term(Rb2p[i, j, 2, 1]) * constant_term(M_[1, 1, j]), order)
                        tmp2119[i, j, 2, 2] = Taylor1(constant_term(Rb2p[i, j, 2, 2]) * constant_term(M_[2, 1, j]), order)
                        tmp2120[i, j, 2, 1] = Taylor1(constant_term(tmp2118[i, j, 2, 1]) + constant_term(tmp2119[i, j, 2, 2]), order)
                        tmp2121[i, j, 2, 3] = Taylor1(constant_term(Rb2p[i, j, 2, 3]) * constant_term(M_[3, 1, j]), order)
                        Gc2p[i, j, 2, 1] = Taylor1(constant_term(tmp2120[i, j, 2, 1]) + constant_term(tmp2121[i, j, 2, 3]), order)
                        tmp2123[i, j, 3, 1] = Taylor1(constant_term(Rb2p[i, j, 3, 1]) * constant_term(M_[1, 1, j]), order)
                        tmp2124[i, j, 3, 2] = Taylor1(constant_term(Rb2p[i, j, 3, 2]) * constant_term(M_[2, 1, j]), order)
                        tmp2125[i, j, 3, 1] = Taylor1(constant_term(tmp2123[i, j, 3, 1]) + constant_term(tmp2124[i, j, 3, 2]), order)
                        tmp2126[i, j, 3, 3] = Taylor1(constant_term(Rb2p[i, j, 3, 3]) * constant_term(M_[3, 1, j]), order)
                        Gc2p[i, j, 3, 1] = Taylor1(constant_term(tmp2125[i, j, 3, 1]) + constant_term(tmp2126[i, j, 3, 3]), order)
                        tmp2128[i, j, 1, 1] = Taylor1(constant_term(Rb2p[i, j, 1, 1]) * constant_term(M_[1, 2, j]), order)
                        tmp2129[i, j, 1, 2] = Taylor1(constant_term(Rb2p[i, j, 1, 2]) * constant_term(M_[2, 2, j]), order)
                        tmp2130[i, j, 1, 1] = Taylor1(constant_term(tmp2128[i, j, 1, 1]) + constant_term(tmp2129[i, j, 1, 2]), order)
                        tmp2131[i, j, 1, 3] = Taylor1(constant_term(Rb2p[i, j, 1, 3]) * constant_term(M_[3, 2, j]), order)
                        Gc2p[i, j, 1, 2] = Taylor1(constant_term(tmp2130[i, j, 1, 1]) + constant_term(tmp2131[i, j, 1, 3]), order)
                        tmp2133[i, j, 2, 1] = Taylor1(constant_term(Rb2p[i, j, 2, 1]) * constant_term(M_[1, 2, j]), order)
                        tmp2134[i, j, 2, 2] = Taylor1(constant_term(Rb2p[i, j, 2, 2]) * constant_term(M_[2, 2, j]), order)
                        tmp2135[i, j, 2, 1] = Taylor1(constant_term(tmp2133[i, j, 2, 1]) + constant_term(tmp2134[i, j, 2, 2]), order)
                        tmp2136[i, j, 2, 3] = Taylor1(constant_term(Rb2p[i, j, 2, 3]) * constant_term(M_[3, 2, j]), order)
                        Gc2p[i, j, 2, 2] = Taylor1(constant_term(tmp2135[i, j, 2, 1]) + constant_term(tmp2136[i, j, 2, 3]), order)
                        tmp2138[i, j, 3, 1] = Taylor1(constant_term(Rb2p[i, j, 3, 1]) * constant_term(M_[1, 2, j]), order)
                        tmp2139[i, j, 3, 2] = Taylor1(constant_term(Rb2p[i, j, 3, 2]) * constant_term(M_[2, 2, j]), order)
                        tmp2140[i, j, 3, 1] = Taylor1(constant_term(tmp2138[i, j, 3, 1]) + constant_term(tmp2139[i, j, 3, 2]), order)
                        tmp2141[i, j, 3, 3] = Taylor1(constant_term(Rb2p[i, j, 3, 3]) * constant_term(M_[3, 2, j]), order)
                        Gc2p[i, j, 3, 2] = Taylor1(constant_term(tmp2140[i, j, 3, 1]) + constant_term(tmp2141[i, j, 3, 3]), order)
                        tmp2143[i, j, 1, 1] = Taylor1(constant_term(Rb2p[i, j, 1, 1]) * constant_term(M_[1, 3, j]), order)
                        tmp2144[i, j, 1, 2] = Taylor1(constant_term(Rb2p[i, j, 1, 2]) * constant_term(M_[2, 3, j]), order)
                        tmp2145[i, j, 1, 1] = Taylor1(constant_term(tmp2143[i, j, 1, 1]) + constant_term(tmp2144[i, j, 1, 2]), order)
                        tmp2146[i, j, 1, 3] = Taylor1(constant_term(Rb2p[i, j, 1, 3]) * constant_term(M_[3, 3, j]), order)
                        Gc2p[i, j, 1, 3] = Taylor1(constant_term(tmp2145[i, j, 1, 1]) + constant_term(tmp2146[i, j, 1, 3]), order)
                        tmp2148[i, j, 2, 1] = Taylor1(constant_term(Rb2p[i, j, 2, 1]) * constant_term(M_[1, 3, j]), order)
                        tmp2149[i, j, 2, 2] = Taylor1(constant_term(Rb2p[i, j, 2, 2]) * constant_term(M_[2, 3, j]), order)
                        tmp2150[i, j, 2, 1] = Taylor1(constant_term(tmp2148[i, j, 2, 1]) + constant_term(tmp2149[i, j, 2, 2]), order)
                        tmp2151[i, j, 2, 3] = Taylor1(constant_term(Rb2p[i, j, 2, 3]) * constant_term(M_[3, 3, j]), order)
                        Gc2p[i, j, 2, 3] = Taylor1(constant_term(tmp2150[i, j, 2, 1]) + constant_term(tmp2151[i, j, 2, 3]), order)
                        tmp2153[i, j, 3, 1] = Taylor1(constant_term(Rb2p[i, j, 3, 1]) * constant_term(M_[1, 3, j]), order)
                        tmp2154[i, j, 3, 2] = Taylor1(constant_term(Rb2p[i, j, 3, 2]) * constant_term(M_[2, 3, j]), order)
                        tmp2155[i, j, 3, 1] = Taylor1(constant_term(tmp2153[i, j, 3, 1]) + constant_term(tmp2154[i, j, 3, 2]), order)
                        tmp2156[i, j, 3, 3] = Taylor1(constant_term(Rb2p[i, j, 3, 3]) * constant_term(M_[3, 3, j]), order)
                        Gc2p[i, j, 3, 3] = Taylor1(constant_term(tmp2155[i, j, 3, 1]) + constant_term(tmp2156[i, j, 3, 3]), order)
                        tmp2158[i, j, 1, 1] = Taylor1(constant_term(F_JCS_ξ[i, j]) * constant_term(Gc2p[i, j, 1, 1]), order)
                        tmp2159[i, j, 2, 1] = Taylor1(constant_term(F_JCS_η[i, j]) * constant_term(Gc2p[i, j, 2, 1]), order)
                        tmp2160[i, j, 1, 1] = Taylor1(constant_term(tmp2158[i, j, 1, 1]) + constant_term(tmp2159[i, j, 2, 1]), order)
                        tmp2161[i, j, 3, 1] = Taylor1(constant_term(F_JCS_ζ[i, j]) * constant_term(Gc2p[i, j, 3, 1]), order)
                        F_JCS_x[i, j] = Taylor1(constant_term(tmp2160[i, j, 1, 1]) + constant_term(tmp2161[i, j, 3, 1]), order)
                        tmp2163[i, j, 1, 2] = Taylor1(constant_term(F_JCS_ξ[i, j]) * constant_term(Gc2p[i, j, 1, 2]), order)
                        tmp2164[i, j, 2, 2] = Taylor1(constant_term(F_JCS_η[i, j]) * constant_term(Gc2p[i, j, 2, 2]), order)
                        tmp2165[i, j, 1, 2] = Taylor1(constant_term(tmp2163[i, j, 1, 2]) + constant_term(tmp2164[i, j, 2, 2]), order)
                        tmp2166[i, j, 3, 2] = Taylor1(constant_term(F_JCS_ζ[i, j]) * constant_term(Gc2p[i, j, 3, 2]), order)
                        F_JCS_y[i, j] = Taylor1(constant_term(tmp2165[i, j, 1, 2]) + constant_term(tmp2166[i, j, 3, 2]), order)
                        tmp2168[i, j, 1, 3] = Taylor1(constant_term(F_JCS_ξ[i, j]) * constant_term(Gc2p[i, j, 1, 3]), order)
                        tmp2169[i, j, 2, 3] = Taylor1(constant_term(F_JCS_η[i, j]) * constant_term(Gc2p[i, j, 2, 3]), order)
                        tmp2170[i, j, 1, 3] = Taylor1(constant_term(tmp2168[i, j, 1, 3]) + constant_term(tmp2169[i, j, 2, 3]), order)
                        tmp2171[i, j, 3, 3] = Taylor1(constant_term(F_JCS_ζ[i, j]) * constant_term(Gc2p[i, j, 3, 3]), order)
                        F_JCS_z[i, j] = Taylor1(constant_term(tmp2170[i, j, 1, 3]) + constant_term(tmp2171[i, j, 3, 3]), order)
                    end
                end
            end
        end
    tmp2173 = Array{Taylor1{_S}}(undef, size(F_JCS_x))
    tmp2173 .= Taylor1(zero(_S), order)
    tmp2175 = Array{Taylor1{_S}}(undef, size(F_JCS_y))
    tmp2175 .= Taylor1(zero(_S), order)
    tmp2177 = Array{Taylor1{_S}}(undef, size(F_JCS_z))
    tmp2177 .= Taylor1(zero(_S), order)
    tmp2179 = Array{Taylor1{_S}}(undef, size(F_JCS_x))
    tmp2179 .= Taylor1(zero(_S), order)
    tmp2181 = Array{Taylor1{_S}}(undef, size(F_JCS_y))
    tmp2181 .= Taylor1(zero(_S), order)
    tmp2183 = Array{Taylor1{_S}}(undef, size(F_JCS_z))
    tmp2183 .= Taylor1(zero(_S), order)
    for j = 1:N_ext
        for i = 1:N_ext
            if i == j
                continue
            else
                if UJ_interaction[i, j]
                    tmp2173[i, j] = Taylor1(constant_term(μ[i]) * constant_term(F_JCS_x[i, j]), order)
                    temp_accX_j[i, j] = Taylor1(constant_term(accX[j]) - constant_term(tmp2173[i, j]), order)
                    accX[j] = Taylor1(identity(constant_term(temp_accX_j[i, j])), order)
                    tmp2175[i, j] = Taylor1(constant_term(μ[i]) * constant_term(F_JCS_y[i, j]), order)
                    temp_accY_j[i, j] = Taylor1(constant_term(accY[j]) - constant_term(tmp2175[i, j]), order)
                    accY[j] = Taylor1(identity(constant_term(temp_accY_j[i, j])), order)
                    tmp2177[i, j] = Taylor1(constant_term(μ[i]) * constant_term(F_JCS_z[i, j]), order)
                    temp_accZ_j[i, j] = Taylor1(constant_term(accZ[j]) - constant_term(tmp2177[i, j]), order)
                    accZ[j] = Taylor1(identity(constant_term(temp_accZ_j[i, j])), order)
                    tmp2179[i, j] = Taylor1(constant_term(μ[j]) * constant_term(F_JCS_x[i, j]), order)
                    temp_accX_i[i, j] = Taylor1(constant_term(accX[i]) + constant_term(tmp2179[i, j]), order)
                    accX[i] = Taylor1(identity(constant_term(temp_accX_i[i, j])), order)
                    tmp2181[i, j] = Taylor1(constant_term(μ[j]) * constant_term(F_JCS_y[i, j]), order)
                    temp_accY_i[i, j] = Taylor1(constant_term(accY[i]) + constant_term(tmp2181[i, j]), order)
                    accY[i] = Taylor1(identity(constant_term(temp_accY_i[i, j])), order)
                    tmp2183[i, j] = Taylor1(constant_term(μ[j]) * constant_term(F_JCS_z[i, j]), order)
                    temp_accZ_i[i, j] = Taylor1(constant_term(accZ[i]) + constant_term(tmp2183[i, j]), order)
                    accZ[i] = Taylor1(identity(constant_term(temp_accZ_i[i, j])), order)
                end
            end
        end
    end
    tmp2189 = Array{Taylor1{_S}}(undef, size(v2))
    tmp2189 .= Taylor1(zero(_S), order)
    tmp2190 = Array{Taylor1{_S}}(undef, size(v2))
    tmp2190 .= Taylor1(zero(_S), order)
    tmp2192 = Array{Taylor1{_S}}(undef, size(vi_dot_vj))
    tmp2192 .= Taylor1(zero(_S), order)
    Xij_t_Ui = Array{Taylor1{_S}}(undef, size(X))
    Xij_t_Ui .= Taylor1(zero(_S), order)
    Yij_t_Vi = Array{Taylor1{_S}}(undef, size(Y))
    Yij_t_Vi .= Taylor1(zero(_S), order)
    Zij_t_Wi = Array{Taylor1{_S}}(undef, size(Z))
    Zij_t_Wi .= Taylor1(zero(_S), order)
    tmp2198 = Array{Taylor1{_S}}(undef, size(Xij_t_Ui))
    tmp2198 .= Taylor1(zero(_S), order)
    Rij_dot_Vi = Array{Taylor1{_S}}(undef, size(tmp2198))
    Rij_dot_Vi .= Taylor1(zero(_S), order)
    tmp2201 = Array{Taylor1{_S}}(undef, size(Rij_dot_Vi))
    tmp2201 .= Taylor1(zero(_S), order)
    pn1t7 = Array{Taylor1{_S}}(undef, size(tmp2201))
    pn1t7 .= Taylor1(zero(_S), order)
    tmp2204 = Array{Taylor1{_S}}(undef, size(pn1t7))
    tmp2204 .= Taylor1(zero(_S), order)
    pn1t2_7 = Array{Taylor1{_S}}(undef, size(ϕs_and_vs))
    pn1t2_7 .= Taylor1(zero(_S), order)
    #= REPL[4]:439 =# Threads.@threads for j = 1:N
            for i = 1:N
                if i == j
                    continue
                else
                    _4ϕj[i, j] = Taylor1(constant_term(4) * constant_term(newtonianNb_Potential[j]), order)
                    ϕi_plus_4ϕj[i, j] = Taylor1(constant_term(newtonianNb_Potential[i]) + constant_term(_4ϕj[i, j]), order)
                    tmp2189[i] = Taylor1(constant_term(2) * constant_term(v2[i]), order)
                    tmp2190[j] = Taylor1(constant_term(v2[j]) + constant_term(tmp2189[i]), order)
                    tmp2192[i, j] = Taylor1(constant_term(4) * constant_term(vi_dot_vj[i, j]), order)
                    sj2_plus_2si2_minus_4vivj[i, j] = Taylor1(constant_term(tmp2190[j]) - constant_term(tmp2192[i, j]), order)
                    ϕs_and_vs[i, j] = Taylor1(constant_term(sj2_plus_2si2_minus_4vivj[i, j]) - constant_term(ϕi_plus_4ϕj[i, j]), order)
                    Xij_t_Ui[i, j] = Taylor1(constant_term(X[i, j]) * constant_term(dq[3i - 2]), order)
                    Yij_t_Vi[i, j] = Taylor1(constant_term(Y[i, j]) * constant_term(dq[3i - 1]), order)
                    Zij_t_Wi[i, j] = Taylor1(constant_term(Z[i, j]) * constant_term(dq[3i]), order)
                    tmp2198[i, j] = Taylor1(constant_term(Xij_t_Ui[i, j]) + constant_term(Yij_t_Vi[i, j]), order)
                    Rij_dot_Vi[i, j] = Taylor1(constant_term(tmp2198[i, j]) + constant_term(Zij_t_Wi[i, j]), order)
                    tmp2201[i, j] = Taylor1(constant_term(Rij_dot_Vi[i, j]) ^ constant_term(2), order)
                    pn1t7[i, j] = Taylor1(constant_term(tmp2201[i, j]) / constant_term(r_p2[i, j]), order)
                    tmp2204[i, j] = Taylor1(constant_term(1.5) * constant_term(pn1t7[i, j]), order)
                    pn1t2_7[i, j] = Taylor1(constant_term(ϕs_and_vs[i, j]) - constant_term(tmp2204[i, j]), order)
                    pn1t1_7[i, j] = Taylor1(constant_term(c_p2) + constant_term(pn1t2_7[i, j]), order)
                    for k = 1:postnewton_iter
                        pn1[i, j, k] = Taylor1(identity(constant_term(zero_q_1)), order)
                        X_t_pn1[i, j, k] = Taylor1(identity(constant_term(zero_q_1)), order)
                        Y_t_pn1[i, j, k] = Taylor1(identity(constant_term(zero_q_1)), order)
                        Z_t_pn1[i, j, k] = Taylor1(identity(constant_term(zero_q_1)), order)
                        pNX_t_pn3[i, j, k] = Taylor1(identity(constant_term(zero_q_1)), order)
                        pNY_t_pn3[i, j, k] = Taylor1(identity(constant_term(zero_q_1)), order)
                        pNZ_t_pn3[i, j, k] = Taylor1(identity(constant_term(zero_q_1)), order)
                        pNX_t_X[i, j, k] = Taylor1(identity(constant_term(zero_q_1)), order)
                        pNY_t_Y[i, j, k] = Taylor1(identity(constant_term(zero_q_1)), order)
                        pNZ_t_Z[i, j, k] = Taylor1(identity(constant_term(zero_q_1)), order)
                    end
                end
            end
            postNewtonX[j, 1] = Taylor1(identity(constant_term(newtonX[j])), order)
            postNewtonY[j, 1] = Taylor1(identity(constant_term(newtonY[j])), order)
            postNewtonZ[j, 1] = Taylor1(identity(constant_term(newtonZ[j])), order)
            for k = 1:postnewton_iter
                pntempX[j, k] = Taylor1(identity(constant_term(zero_q_1)), order)
                pntempY[j, k] = Taylor1(identity(constant_term(zero_q_1)), order)
                pntempZ[j, k] = Taylor1(identity(constant_term(zero_q_1)), order)
            end
        end
    tmp2211 = Array{Taylor1{_S}}(undef, size(pNX_t_X))
    tmp2211 .= Taylor1(zero(_S), order)
    tmp2212 = Array{Taylor1{_S}}(undef, size(tmp2211))
    tmp2212 .= Taylor1(zero(_S), order)
    tmp2213 = Array{Taylor1{_S}}(undef, size(tmp2212))
    tmp2213 .= Taylor1(zero(_S), order)
    tmp2221 = Array{Taylor1{_S}}(undef, size(pNX_t_pn3))
    tmp2221 .= Taylor1(zero(_S), order)
    termpnx = Array{Taylor1{_S}}(undef, size(X_t_pn1))
    termpnx .= Taylor1(zero(_S), order)
    sumpnx = Array{Taylor1{_S}}(undef, size(termpnx))
    sumpnx .= Taylor1(zero(_S), order)
    tmp2224 = Array{Taylor1{_S}}(undef, size(pNY_t_pn3))
    tmp2224 .= Taylor1(zero(_S), order)
    termpny = Array{Taylor1{_S}}(undef, size(Y_t_pn1))
    termpny .= Taylor1(zero(_S), order)
    sumpny = Array{Taylor1{_S}}(undef, size(termpny))
    sumpny .= Taylor1(zero(_S), order)
    tmp2227 = Array{Taylor1{_S}}(undef, size(pNZ_t_pn3))
    tmp2227 .= Taylor1(zero(_S), order)
    termpnz = Array{Taylor1{_S}}(undef, size(Z_t_pn1))
    termpnz .= Taylor1(zero(_S), order)
    sumpnz = Array{Taylor1{_S}}(undef, size(termpnz))
    sumpnz .= Taylor1(zero(_S), order)
    for k = 1:postnewton_iter
        #= REPL[4]:484 =# Threads.@threads for j = 1:N
                for i = 1:N
                    if i == j
                        continue
                    else
                        pNX_t_X[i, j, k] = Taylor1(constant_term(postNewtonX[i, k]) * constant_term(X[i, j]), order)
                        pNY_t_Y[i, j, k] = Taylor1(constant_term(postNewtonY[i, k]) * constant_term(Y[i, j]), order)
                        pNZ_t_Z[i, j, k] = Taylor1(constant_term(postNewtonZ[i, k]) * constant_term(Z[i, j]), order)
                        tmp2211[i, j, k] = Taylor1(constant_term(pNX_t_X[i, j, k]) + constant_term(pNY_t_Y[i, j, k]), order)
                        tmp2212[i, j, k] = Taylor1(constant_term(tmp2211[i, j, k]) + constant_term(pNZ_t_Z[i, j, k]), order)
                        tmp2213[i, j, k] = Taylor1(constant_term(0.5) * constant_term(tmp2212[i, j, k]), order)
                        pn1[i, j, k] = Taylor1(constant_term(pn1t1_7[i, j]) + constant_term(tmp2213[i, j, k]), order)
                        X_t_pn1[i, j, k] = Taylor1(constant_term(newton_acc_X[i, j]) * constant_term(pn1[i, j, k]), order)
                        Y_t_pn1[i, j, k] = Taylor1(constant_term(newton_acc_Y[i, j]) * constant_term(pn1[i, j, k]), order)
                        Z_t_pn1[i, j, k] = Taylor1(constant_term(newton_acc_Z[i, j]) * constant_term(pn1[i, j, k]), order)
                        pNX_t_pn3[i, j, k] = Taylor1(constant_term(postNewtonX[i, k]) * constant_term(pn3[i, j]), order)
                        pNY_t_pn3[i, j, k] = Taylor1(constant_term(postNewtonY[i, k]) * constant_term(pn3[i, j]), order)
                        pNZ_t_pn3[i, j, k] = Taylor1(constant_term(postNewtonZ[i, k]) * constant_term(pn3[i, j]), order)
                        tmp2221[i, j, k] = Taylor1(constant_term(U_t_pn2[i, j]) + constant_term(pNX_t_pn3[i, j, k]), order)
                        termpnx[i, j, k] = Taylor1(constant_term(X_t_pn1[i, j, k]) + constant_term(tmp2221[i, j, k]), order)
                        sumpnx[i, j, k] = Taylor1(constant_term(pntempX[j, k]) + constant_term(termpnx[i, j, k]), order)
                        pntempX[j, k] = Taylor1(identity(constant_term(sumpnx[i, j, k])), order)
                        tmp2224[i, j, k] = Taylor1(constant_term(V_t_pn2[i, j]) + constant_term(pNY_t_pn3[i, j, k]), order)
                        termpny[i, j, k] = Taylor1(constant_term(Y_t_pn1[i, j, k]) + constant_term(tmp2224[i, j, k]), order)
                        sumpny[i, j, k] = Taylor1(constant_term(pntempY[j, k]) + constant_term(termpny[i, j, k]), order)
                        pntempY[j, k] = Taylor1(identity(constant_term(sumpny[i, j, k])), order)
                        tmp2227[i, j, k] = Taylor1(constant_term(W_t_pn2[i, j]) + constant_term(pNZ_t_pn3[i, j, k]), order)
                        termpnz[i, j, k] = Taylor1(constant_term(Z_t_pn1[i, j, k]) + constant_term(tmp2227[i, j, k]), order)
                        sumpnz[i, j, k] = Taylor1(constant_term(pntempZ[j, k]) + constant_term(termpnz[i, j, k]), order)
                        pntempZ[j, k] = Taylor1(identity(constant_term(sumpnz[i, j, k])), order)
                    end
                end
                postNewtonX[j, k + 1] = Taylor1(constant_term(pntempX[j, k]) * constant_term(c_m2), order)
                postNewtonY[j, k + 1] = Taylor1(constant_term(pntempY[j, k]) * constant_term(c_m2), order)
                postNewtonZ[j, k + 1] = Taylor1(constant_term(pntempZ[j, k]) * constant_term(c_m2), order)
            end
    end
    #= REPL[4]:521 =# Threads.@threads for i = 1:N_ext
            dq[3 * (N + i) - 2] = Taylor1(constant_term(postNewtonX[i, postnewton_iter + 1]) + constant_term(accX[i]), order)
            dq[3 * (N + i) - 1] = Taylor1(constant_term(postNewtonY[i, postnewton_iter + 1]) + constant_term(accY[i]), order)
            dq[3 * (N + i)] = Taylor1(constant_term(postNewtonZ[i, postnewton_iter + 1]) + constant_term(accZ[i]), order)
        end
    #= REPL[4]:526 =# Threads.@threads for i = N_ext + 1:N
            dq[3 * (N + i) - 2] = Taylor1(identity(constant_term(postNewtonX[i, postnewton_iter + 1])), order)
            dq[3 * (N + i) - 1] = Taylor1(identity(constant_term(postNewtonY[i, postnewton_iter + 1])), order)
            dq[3 * (N + i)] = Taylor1(identity(constant_term(postNewtonZ[i, postnewton_iter + 1])), order)
        end
    for __idx = eachindex(q)
        (q[__idx]).coeffs[2] = (dq[__idx]).coeffs[1]
    end
    for ord = 1:order - 1
        ordnext = ord + 1
        TaylorSeries.identity!(J2_t[su], J2S_t, ord)
        TaylorSeries.identity!(J2_t[ea], J2E_t, ord)
        #= REPL[4]:180 =# Threads.@threads for j = 1:N
                TaylorSeries.identity!(newtonX[j], zero_q_1, ord)
                TaylorSeries.identity!(newtonY[j], zero_q_1, ord)
                TaylorSeries.identity!(newtonZ[j], zero_q_1, ord)
                TaylorSeries.identity!(newtonianNb_Potential[j], zero_q_1, ord)
                TaylorSeries.identity!(dq[3j - 2], q[3 * (N + j) - 2], ord)
                TaylorSeries.identity!(dq[3j - 1], q[3 * (N + j) - 1], ord)
                TaylorSeries.identity!(dq[3j], q[3 * (N + j)], ord)
            end
        #= REPL[4]:190 =# Threads.@threads for j = 1:N_ext
                TaylorSeries.identity!(accX[j], zero_q_1, ord)
                TaylorSeries.identity!(accY[j], zero_q_1, ord)
                TaylorSeries.identity!(accZ[j], zero_q_1, ord)
            end
        #= REPL[4]:197 =# Threads.@threads for j = 1:N
                for i = 1:N
                    if i == j
                        continue
                    else
                        TaylorSeries.subst!(X[i, j], q[3i - 2], q[3j - 2], ord)
                        TaylorSeries.subst!(Y[i, j], q[3i - 1], q[3j - 1], ord)
                        TaylorSeries.subst!(Z[i, j], q[3i], q[3j], ord)
                        TaylorSeries.subst!(U[i, j], dq[3i - 2], dq[3j - 2], ord)
                        TaylorSeries.subst!(V[i, j], dq[3i - 1], dq[3j - 1], ord)
                        TaylorSeries.subst!(W[i, j], dq[3i], dq[3j], ord)
                        TaylorSeries.mul!(tmp1890[3j - 2], 4, dq[3j - 2], ord)
                        TaylorSeries.mul!(tmp1892[3i - 2], 3, dq[3i - 2], ord)
                        TaylorSeries.subst!(_4U_m_3X[i, j], tmp1890[3j - 2], tmp1892[3i - 2], ord)
                        TaylorSeries.mul!(tmp1895[3j - 1], 4, dq[3j - 1], ord)
                        TaylorSeries.mul!(tmp1897[3i - 1], 3, dq[3i - 1], ord)
                        TaylorSeries.subst!(_4V_m_3Y[i, j], tmp1895[3j - 1], tmp1897[3i - 1], ord)
                        TaylorSeries.mul!(tmp1900[3j], 4, dq[3j], ord)
                        TaylorSeries.mul!(tmp1902[3i], 3, dq[3i], ord)
                        TaylorSeries.subst!(_4W_m_3Z[i, j], tmp1900[3j], tmp1902[3i], ord)
                        TaylorSeries.mul!(pn2x[i, j], X[i, j], _4U_m_3X[i, j], ord)
                        TaylorSeries.mul!(pn2y[i, j], Y[i, j], _4V_m_3Y[i, j], ord)
                        TaylorSeries.mul!(pn2z[i, j], Z[i, j], _4W_m_3Z[i, j], ord)
                        TaylorSeries.mul!(UU[i, j], dq[3i - 2], dq[3j - 2], ord)
                        TaylorSeries.mul!(VV[i, j], dq[3i - 1], dq[3j - 1], ord)
                        TaylorSeries.mul!(WW[i, j], dq[3i], dq[3j], ord)
                        TaylorSeries.add!(tmp1910[i, j], UU[i, j], VV[i, j], ord)
                        TaylorSeries.add!(vi_dot_vj[i, j], tmp1910[i, j], WW[i, j], ord)
                        TaylorSeries.pow!(tmp1913[i, j], X[i, j], 2, ord)
                        TaylorSeries.pow!(tmp1915[i, j], Y[i, j], 2, ord)
                        TaylorSeries.add!(tmp1916[i, j], tmp1913[i, j], tmp1915[i, j], ord)
                        TaylorSeries.pow!(tmp1918[i, j], Z[i, j], 2, ord)
                        TaylorSeries.add!(r_p2[i, j], tmp1916[i, j], tmp1918[i, j], ord)
                        TaylorSeries.sqrt!(r_p1d2[i, j], r_p2[i, j], ord)
                        TaylorSeries.pow!(r_p3d2[i, j], r_p2[i, j], 1.5, ord)
                        TaylorSeries.pow!(r_p7d2[i, j], r_p2[i, j], 3.5, ord)
                        TaylorSeries.div!(newtonianCoeff[i, j], μ[i], r_p3d2[i, j], ord)
                        TaylorSeries.add!(tmp1926[i, j], pn2x[i, j], pn2y[i, j], ord)
                        TaylorSeries.add!(tmp1927[i, j], tmp1926[i, j], pn2z[i, j], ord)
                        TaylorSeries.mul!(pn2[i, j], newtonianCoeff[i, j], tmp1927[i, j], ord)
                        TaylorSeries.mul!(newton_acc_X[i, j], X[i, j], newtonianCoeff[i, j], ord)
                        TaylorSeries.mul!(newton_acc_Y[i, j], Y[i, j], newtonianCoeff[i, j], ord)
                        TaylorSeries.mul!(newton_acc_Z[i, j], Z[i, j], newtonianCoeff[i, j], ord)
                        TaylorSeries.div!(newtonian1b_Potential[i, j], μ[i], r_p1d2[i, j], ord)
                        TaylorSeries.mul!(pn3[i, j], 3.5, newtonian1b_Potential[i, j], ord)
                        TaylorSeries.mul!(U_t_pn2[i, j], pn2[i, j], U[i, j], ord)
                        TaylorSeries.mul!(V_t_pn2[i, j], pn2[i, j], V[i, j], ord)
                        TaylorSeries.mul!(W_t_pn2[i, j], pn2[i, j], W[i, j], ord)
                        TaylorSeries.mul!(tmp1938[i, j], X[i, j], newtonianCoeff[i, j], ord)
                        TaylorSeries.add!(temp_001[i, j], newtonX[j], tmp1938[i, j], ord)
                        TaylorSeries.identity!(newtonX[j], temp_001[i, j], ord)
                        TaylorSeries.mul!(tmp1940[i, j], Y[i, j], newtonianCoeff[i, j], ord)
                        TaylorSeries.add!(temp_002[i, j], newtonY[j], tmp1940[i, j], ord)
                        TaylorSeries.identity!(newtonY[j], temp_002[i, j], ord)
                        TaylorSeries.mul!(tmp1942[i, j], Z[i, j], newtonianCoeff[i, j], ord)
                        TaylorSeries.add!(temp_003[i, j], newtonZ[j], tmp1942[i, j], ord)
                        TaylorSeries.identity!(newtonZ[j], temp_003[i, j], ord)
                        TaylorSeries.add!(temp_004[i, j], newtonianNb_Potential[j], newtonian1b_Potential[i, j], ord)
                        TaylorSeries.identity!(newtonianNb_Potential[j], temp_004[i, j], ord)
                    end
                end
                TaylorSeries.pow!(tmp1946[3j - 2], dq[3j - 2], 2, ord)
                TaylorSeries.pow!(tmp1948[3j - 1], dq[3j - 1], 2, ord)
                TaylorSeries.add!(tmp1949[3j - 2], tmp1946[3j - 2], tmp1948[3j - 1], ord)
                TaylorSeries.pow!(tmp1951[3j], dq[3j], 2, ord)
                TaylorSeries.add!(v2[j], tmp1949[3j - 2], tmp1951[3j], ord)
            end
        TaylorSeries.add!(tmp1953, ITM_t[1, 1], ITM_t[2, 2], ord)
        TaylorSeries.div!(tmp1955, tmp1953, 2, ord)
        TaylorSeries.subst!(tmp1956, ITM_t[3, 3], tmp1955, ord)
        TaylorSeries.div!(J2M_t, tmp1956, μ[mo], ord)
        TaylorSeries.subst!(tmp1958, ITM_t[2, 2], ITM_t[1, 1], ord)
        TaylorSeries.div!(tmp1959, tmp1958, μ[mo], ord)
        TaylorSeries.div!(C22M_t, tmp1959, 4, ord)
        TaylorSeries.subst!(tmp1962, ITM_t[1, 3], ord)
        TaylorSeries.div!(C21M_t, tmp1962, μ[mo], ord)
        TaylorSeries.subst!(tmp1964, ITM_t[3, 2], ord)
        TaylorSeries.div!(S21M_t, tmp1964, μ[mo], ord)
        TaylorSeries.subst!(tmp1966, ITM_t[2, 1], ord)
        TaylorSeries.div!(tmp1967, tmp1966, μ[mo], ord)
        TaylorSeries.div!(S22M_t, tmp1967, 2, ord)
        TaylorSeries.identity!(J2_t[mo], J2M_t, ord)
        #= REPL[4]:265 =# Threads.@threads for j = 1:N_ext
                for i = 1:N_ext
                    if i == j
                        continue
                    else
                        if UJ_interaction[i, j]
                            TaylorSeries.mul!(X_bf_1[i, j], X[i, j], M_[1, 1, j], ord)
                            TaylorSeries.mul!(X_bf_2[i, j], Y[i, j], M_[1, 2, j], ord)
                            TaylorSeries.mul!(X_bf_3[i, j], Z[i, j], M_[1, 3, j], ord)
                            TaylorSeries.mul!(Y_bf_1[i, j], X[i, j], M_[2, 1, j], ord)
                            TaylorSeries.mul!(Y_bf_2[i, j], Y[i, j], M_[2, 2, j], ord)
                            TaylorSeries.mul!(Y_bf_3[i, j], Z[i, j], M_[2, 3, j], ord)
                            TaylorSeries.mul!(Z_bf_1[i, j], X[i, j], M_[3, 1, j], ord)
                            TaylorSeries.mul!(Z_bf_2[i, j], Y[i, j], M_[3, 2, j], ord)
                            TaylorSeries.mul!(Z_bf_3[i, j], Z[i, j], M_[3, 3, j], ord)
                            TaylorSeries.add!(tmp1979[i, j], X_bf_1[i, j], X_bf_2[i, j], ord)
                            TaylorSeries.add!(X_bf[i, j], tmp1979[i, j], X_bf_3[i, j], ord)
                            TaylorSeries.add!(tmp1981[i, j], Y_bf_1[i, j], Y_bf_2[i, j], ord)
                            TaylorSeries.add!(Y_bf[i, j], tmp1981[i, j], Y_bf_3[i, j], ord)
                            TaylorSeries.add!(tmp1983[i, j], Z_bf_1[i, j], Z_bf_2[i, j], ord)
                            TaylorSeries.add!(Z_bf[i, j], tmp1983[i, j], Z_bf_3[i, j], ord)
                            TaylorSeries.div!(sin_ϕ[i, j], Z_bf[i, j], r_p1d2[i, j], ord)
                            TaylorSeries.pow!(tmp1987[i, j], X_bf[i, j], 2, ord)
                            TaylorSeries.pow!(tmp1989[i, j], Y_bf[i, j], 2, ord)
                            TaylorSeries.add!(tmp1990[i, j], tmp1987[i, j], tmp1989[i, j], ord)
                            TaylorSeries.sqrt!(r_xy[i, j], tmp1990[i, j], ord)
                            TaylorSeries.div!(cos_ϕ[i, j], r_xy[i, j], r_p1d2[i, j], ord)
                            TaylorSeries.div!(sin_λ[i, j], Y_bf[i, j], r_xy[i, j], ord)
                            TaylorSeries.div!(cos_λ[i, j], X_bf[i, j], r_xy[i, j], ord)
                            TaylorSeries.identity!(P_n[i, j, 1], one_t, ord)
                            TaylorSeries.identity!(P_n[i, j, 2], sin_ϕ[i, j], ord)
                            TaylorSeries.identity!(dP_n[i, j, 1], zero_q_1, ord)
                            TaylorSeries.identity!(dP_n[i, j, 2], one_t, ord)
                            for n = 2:n1SEM[j]
                                TaylorSeries.mul!(tmp1995[i, j, n], P_n[i, j, n], sin_ϕ[i, j], ord)
                                TaylorSeries.mul!(tmp1996[i, j, n], tmp1995[i, j, n], fact1_jsem[n], ord)
                                TaylorSeries.mul!(tmp1997[i, j, n - 1], P_n[i, j, n - 1], fact2_jsem[n], ord)
                                TaylorSeries.subst!(P_n[i, j, n + 1], tmp1996[i, j, n], tmp1997[i, j, n - 1], ord)
                                TaylorSeries.mul!(tmp1999[i, j, n], dP_n[i, j, n], sin_ϕ[i, j], ord)
                                TaylorSeries.mul!(tmp2000[i, j, n], P_n[i, j, n], fact3_jsem[n], ord)
                                TaylorSeries.add!(dP_n[i, j, n + 1], tmp1999[i, j, n], tmp2000[i, j, n], ord)
                                TaylorSeries.pow!(temp_rn[i, j, n], r_p1d2[i, j], fact5_jsem[n], ord)
                            end
                            TaylorSeries.pow!(r_p4[i, j], r_p2[i, j], 2, ord)
                            TaylorSeries.mul!(tmp2005[i, j, 3], P_n[i, j, 3], fact4_jsem[2], ord)
                            TaylorSeries.mul!(tmp2006[i, j, 3], tmp2005[i, j, 3], J2_t[j], ord)
                            TaylorSeries.div!(F_J_ξ[i, j], tmp2006[i, j, 3], r_p4[i, j], ord)
                            TaylorSeries.subst!(tmp2008[i, j, 3], dP_n[i, j, 3], ord)
                            TaylorSeries.mul!(tmp2009[i, j, 3], tmp2008[i, j, 3], cos_ϕ[i, j], ord)
                            TaylorSeries.mul!(tmp2010[i, j, 3], tmp2009[i, j, 3], J2_t[j], ord)
                            TaylorSeries.div!(F_J_ζ[i, j], tmp2010[i, j, 3], r_p4[i, j], ord)
                            TaylorSeries.identity!(F_J_ξ_36[i, j], zero_q_1, ord)
                            TaylorSeries.identity!(F_J_ζ_36[i, j], zero_q_1, ord)
                            for n = 3:n1SEM[j]
                                TaylorSeries.mul!(tmp2012[i, j, n + 1], P_n[i, j, n + 1], fact4_jsem[n], ord)
                                TaylorSeries.mul!(tmp2013[i, j, n + 1], tmp2012[i, j, n + 1], JSEM[j, n], ord)
                                TaylorSeries.div!(tmp2014[i, j, n + 1], tmp2013[i, j, n + 1], temp_rn[i, j, n], ord)
                                TaylorSeries.add!(temp_fjξ[i, j, n], tmp2014[i, j, n + 1], F_J_ξ_36[i, j], ord)
                                TaylorSeries.subst!(tmp2016[i, j, n + 1], dP_n[i, j, n + 1], ord)
                                TaylorSeries.mul!(tmp2017[i, j, n + 1], tmp2016[i, j, n + 1], cos_ϕ[i, j], ord)
                                TaylorSeries.mul!(tmp2018[i, j, n + 1], tmp2017[i, j, n + 1], JSEM[j, n], ord)
                                TaylorSeries.div!(tmp2019[i, j, n + 1], tmp2018[i, j, n + 1], temp_rn[i, j, n], ord)
                                TaylorSeries.add!(temp_fjζ[i, j, n], tmp2019[i, j, n + 1], F_J_ζ_36[i, j], ord)
                                TaylorSeries.identity!(F_J_ξ_36[i, j], temp_fjξ[i, j, n], ord)
                                TaylorSeries.identity!(F_J_ζ_36[i, j], temp_fjζ[i, j, n], ord)
                            end
                            if j == mo
                                for m = 1:n1SEM[mo]
                                    if m == 1
                                        TaylorSeries.identity!(sin_mλ[i, j, 1], sin_λ[i, j], ord)
                                        TaylorSeries.identity!(cos_mλ[i, j, 1], cos_λ[i, j], ord)
                                        TaylorSeries.identity!(secϕ_P_nm[i, j, 1, 1], one_t, ord)
                                    else
                                        TaylorSeries.mul!(tmp2021[i, j, 1], sin_mλ[i, j, 1], cos_mλ[i, j, m - 1], ord)
                                        TaylorSeries.mul!(tmp2022[i, j, 1], cos_mλ[i, j, 1], sin_mλ[i, j, m - 1], ord)
                                        TaylorSeries.add!(sin_mλ[i, j, m], tmp2021[i, j, 1], tmp2022[i, j, 1], ord)
                                        TaylorSeries.mul!(tmp2024[i, j, 1], cos_mλ[i, j, 1], cos_mλ[i, j, m - 1], ord)
                                        TaylorSeries.mul!(tmp2025[i, j, 1], sin_mλ[i, j, 1], sin_mλ[i, j, m - 1], ord)
                                        TaylorSeries.subst!(cos_mλ[i, j, m], tmp2024[i, j, 1], tmp2025[i, j, 1], ord)
                                        TaylorSeries.mul!(tmp2027[i, j, m - 1, m - 1], secϕ_P_nm[i, j, m - 1, m - 1], cos_ϕ[i, j], ord)
                                        TaylorSeries.mul!(secϕ_P_nm[i, j, m, m], tmp2027[i, j, m - 1, m - 1], lnm5[m], ord)
                                        TaylorSeries.mul!(P_nm[i, j, m, m], secϕ_P_nm[i, j, m, m], cos_ϕ[i, j], ord)
                                        TaylorSeries.mul!(tmp2030[i, j, m, m], secϕ_P_nm[i, j, m, m], sin_ϕ[i, j], ord)
                                        TaylorSeries.mul!(cosϕ_dP_nm[i, j, m, m], tmp2030[i, j, m, m], lnm3[m], ord)
                                    end
                                    for n = m + 1:n1SEM[mo]
                                        if n == m + 1
                                            TaylorSeries.mul!(tmp2032[i, j, n - 1, m], secϕ_P_nm[i, j, n - 1, m], sin_ϕ[i, j], ord)
                                            TaylorSeries.mul!(secϕ_P_nm[i, j, n, m], tmp2032[i, j, n - 1, m], lnm1[n, m], ord)
                                        else
                                            TaylorSeries.mul!(tmp2034[i, j, n - 1, m], secϕ_P_nm[i, j, n - 1, m], sin_ϕ[i, j], ord)
                                            TaylorSeries.mul!(tmp2035[i, j, n - 1, m], tmp2034[i, j, n - 1, m], lnm1[n, m], ord)
                                            TaylorSeries.mul!(tmp2036[i, j, n - 2, m], secϕ_P_nm[i, j, n - 2, m], lnm2[n, m], ord)
                                            TaylorSeries.add!(secϕ_P_nm[i, j, n, m], tmp2035[i, j, n - 1, m], tmp2036[i, j, n - 2, m], ord)
                                        end
                                        TaylorSeries.mul!(P_nm[i, j, n, m], secϕ_P_nm[i, j, n, m], cos_ϕ[i, j], ord)
                                        TaylorSeries.mul!(tmp2039[i, j, n, m], secϕ_P_nm[i, j, n, m], sin_ϕ[i, j], ord)
                                        TaylorSeries.mul!(tmp2040[i, j, n, m], tmp2039[i, j, n, m], lnm3[n], ord)
                                        TaylorSeries.mul!(tmp2041[i, j, n - 1, m], secϕ_P_nm[i, j, n - 1, m], lnm4[n, m], ord)
                                        TaylorSeries.add!(cosϕ_dP_nm[i, j, n, m], tmp2040[i, j, n, m], tmp2041[i, j, n - 1, m], ord)
                                    end
                                end
                                TaylorSeries.mul!(tmp2043[i, j, 2, 1], P_nm[i, j, 2, 1], lnm6[2], ord)
                                TaylorSeries.mul!(tmp2044[i, j, 1], C21M_t, cos_mλ[i, j, 1], ord)
                                TaylorSeries.mul!(tmp2045[i, j, 1], S21M_t, sin_mλ[i, j, 1], ord)
                                TaylorSeries.add!(tmp2046[i, j, 1], tmp2044[i, j, 1], tmp2045[i, j, 1], ord)
                                TaylorSeries.mul!(tmp2047[i, j, 2, 1], tmp2043[i, j, 2, 1], tmp2046[i, j, 1], ord)
                                TaylorSeries.mul!(tmp2048[i, j, 2, 2], P_nm[i, j, 2, 2], lnm6[2], ord)
                                TaylorSeries.mul!(tmp2049[i, j, 2], C22M_t, cos_mλ[i, j, 2], ord)
                                TaylorSeries.mul!(tmp2050[i, j, 2], S22M_t, sin_mλ[i, j, 2], ord)
                                TaylorSeries.add!(tmp2051[i, j, 2], tmp2049[i, j, 2], tmp2050[i, j, 2], ord)
                                TaylorSeries.mul!(tmp2052[i, j, 2, 2], tmp2048[i, j, 2, 2], tmp2051[i, j, 2], ord)
                                TaylorSeries.add!(tmp2053[i, j, 2, 1], tmp2047[i, j, 2, 1], tmp2052[i, j, 2, 2], ord)
                                TaylorSeries.div!(F_CS_ξ[i, j], tmp2053[i, j, 2, 1], r_p4[i, j], ord)
                                TaylorSeries.mul!(tmp2055[i, j, 2, 1], secϕ_P_nm[i, j, 2, 1], lnm7[1], ord)
                                TaylorSeries.mul!(tmp2056[i, j, 1], S21M_t, cos_mλ[i, j, 1], ord)
                                TaylorSeries.mul!(tmp2057[i, j, 1], C21M_t, sin_mλ[i, j, 1], ord)
                                TaylorSeries.subst!(tmp2058[i, j, 1], tmp2056[i, j, 1], tmp2057[i, j, 1], ord)
                                TaylorSeries.mul!(tmp2059[i, j, 2, 1], tmp2055[i, j, 2, 1], tmp2058[i, j, 1], ord)
                                TaylorSeries.mul!(tmp2060[i, j, 2, 2], secϕ_P_nm[i, j, 2, 2], lnm7[2], ord)
                                TaylorSeries.mul!(tmp2061[i, j, 2], S22M_t, cos_mλ[i, j, 2], ord)
                                TaylorSeries.mul!(tmp2062[i, j, 2], C22M_t, sin_mλ[i, j, 2], ord)
                                TaylorSeries.subst!(tmp2063[i, j, 2], tmp2061[i, j, 2], tmp2062[i, j, 2], ord)
                                TaylorSeries.mul!(tmp2064[i, j, 2, 2], tmp2060[i, j, 2, 2], tmp2063[i, j, 2], ord)
                                TaylorSeries.add!(tmp2065[i, j, 2, 1], tmp2059[i, j, 2, 1], tmp2064[i, j, 2, 2], ord)
                                TaylorSeries.div!(F_CS_η[i, j], tmp2065[i, j, 2, 1], r_p4[i, j], ord)
                                TaylorSeries.mul!(tmp2067[i, j, 1], C21M_t, cos_mλ[i, j, 1], ord)
                                TaylorSeries.mul!(tmp2068[i, j, 1], S21M_t, sin_mλ[i, j, 1], ord)
                                TaylorSeries.add!(tmp2069[i, j, 1], tmp2067[i, j, 1], tmp2068[i, j, 1], ord)
                                TaylorSeries.mul!(tmp2070[i, j, 2, 1], cosϕ_dP_nm[i, j, 2, 1], tmp2069[i, j, 1], ord)
                                TaylorSeries.mul!(tmp2071[i, j, 2], C22M_t, cos_mλ[i, j, 2], ord)
                                TaylorSeries.mul!(tmp2072[i, j, 2], S22M_t, sin_mλ[i, j, 2], ord)
                                TaylorSeries.add!(tmp2073[i, j, 2], tmp2071[i, j, 2], tmp2072[i, j, 2], ord)
                                TaylorSeries.mul!(tmp2074[i, j, 2, 2], cosϕ_dP_nm[i, j, 2, 2], tmp2073[i, j, 2], ord)
                                TaylorSeries.add!(tmp2075[i, j, 2, 1], tmp2070[i, j, 2, 1], tmp2074[i, j, 2, 2], ord)
                                TaylorSeries.div!(F_CS_ζ[i, j], tmp2075[i, j, 2, 1], r_p4[i, j], ord)
                                TaylorSeries.identity!(F_CS_ξ_36[i, j], zero_q_1, ord)
                                TaylorSeries.identity!(F_CS_η_36[i, j], zero_q_1, ord)
                                TaylorSeries.identity!(F_CS_ζ_36[i, j], zero_q_1, ord)
                                for n = 3:n1SEM[mo]
                                    for m = 1:n
                                        TaylorSeries.mul!(tmp2077[i, j, n, m], P_nm[i, j, n, m], lnm6[n], ord)
                                        TaylorSeries.mul!(tmp2078[i, j, m], cos_mλ[i, j, m], CM[n, m], ord)
                                        TaylorSeries.mul!(tmp2079[i, j, m], sin_mλ[i, j, m], SM[n, m], ord)
                                        TaylorSeries.add!(tmp2080[i, j, m], tmp2078[i, j, m], tmp2079[i, j, m], ord)
                                        TaylorSeries.mul!(tmp2081[i, j, n, m], tmp2077[i, j, n, m], tmp2080[i, j, m], ord)
                                        TaylorSeries.div!(tmp2082[i, j, n, m], tmp2081[i, j, n, m], temp_rn[i, j, n], ord)
                                        TaylorSeries.add!(temp_CS_ξ[i, j, n, m], tmp2082[i, j, n, m], F_CS_ξ_36[i, j], ord)
                                        TaylorSeries.mul!(tmp2084[i, j, n, m], secϕ_P_nm[i, j, n, m], lnm7[m], ord)
                                        TaylorSeries.mul!(tmp2085[i, j, m], cos_mλ[i, j, m], SM[n, m], ord)
                                        TaylorSeries.mul!(tmp2086[i, j, m], sin_mλ[i, j, m], CM[n, m], ord)
                                        TaylorSeries.subst!(tmp2087[i, j, m], tmp2085[i, j, m], tmp2086[i, j, m], ord)
                                        TaylorSeries.mul!(tmp2088[i, j, n, m], tmp2084[i, j, n, m], tmp2087[i, j, m], ord)
                                        TaylorSeries.div!(tmp2089[i, j, n, m], tmp2088[i, j, n, m], temp_rn[i, j, n], ord)
                                        TaylorSeries.add!(temp_CS_η[i, j, n, m], tmp2089[i, j, n, m], F_CS_η_36[i, j], ord)
                                        TaylorSeries.mul!(tmp2091[i, j, m], cos_mλ[i, j, m], CM[n, m], ord)
                                        TaylorSeries.mul!(tmp2092[i, j, m], sin_mλ[i, j, m], SM[n, m], ord)
                                        TaylorSeries.add!(tmp2093[i, j, m], tmp2091[i, j, m], tmp2092[i, j, m], ord)
                                        TaylorSeries.mul!(tmp2094[i, j, n, m], cosϕ_dP_nm[i, j, n, m], tmp2093[i, j, m], ord)
                                        TaylorSeries.div!(tmp2095[i, j, n, m], tmp2094[i, j, n, m], temp_rn[i, j, n], ord)
                                        TaylorSeries.add!(temp_CS_ζ[i, j, n, m], tmp2095[i, j, n, m], F_CS_ζ_36[i, j], ord)
                                        TaylorSeries.identity!(F_CS_ξ_36[i, j], temp_CS_ξ[i, j, n, m], ord)
                                        TaylorSeries.identity!(F_CS_η_36[i, j], temp_CS_η[i, j, n, m], ord)
                                        TaylorSeries.identity!(F_CS_ζ_36[i, j], temp_CS_ζ[i, j, n, m], ord)
                                    end
                                end
                                TaylorSeries.add!(tmp2097[i, j], F_J_ξ[i, j], F_J_ξ_36[i, j], ord)
                                TaylorSeries.add!(tmp2098[i, j], F_CS_ξ[i, j], F_CS_ξ_36[i, j], ord)
                                TaylorSeries.add!(F_JCS_ξ[i, j], tmp2097[i, j], tmp2098[i, j], ord)
                                TaylorSeries.add!(F_JCS_η[i, j], F_CS_η[i, j], F_CS_η_36[i, j], ord)
                                TaylorSeries.add!(tmp2101[i, j], F_J_ζ[i, j], F_J_ζ_36[i, j], ord)
                                TaylorSeries.add!(tmp2102[i, j], F_CS_ζ[i, j], F_CS_ζ_36[i, j], ord)
                                TaylorSeries.add!(F_JCS_ζ[i, j], tmp2101[i, j], tmp2102[i, j], ord)
                            else
                                TaylorSeries.add!(F_JCS_ξ[i, j], F_J_ξ[i, j], F_J_ξ_36[i, j], ord)
                                TaylorSeries.identity!(F_JCS_η[i, j], zero_q_1, ord)
                                TaylorSeries.add!(F_JCS_ζ[i, j], F_J_ζ[i, j], F_J_ζ_36[i, j], ord)
                            end
                            TaylorSeries.mul!(Rb2p[i, j, 1, 1], cos_ϕ[i, j], cos_λ[i, j], ord)
                            TaylorSeries.subst!(Rb2p[i, j, 2, 1], sin_λ[i, j], ord)
                            TaylorSeries.subst!(tmp2108[i, j], sin_ϕ[i, j], ord)
                            TaylorSeries.mul!(Rb2p[i, j, 3, 1], tmp2108[i, j], cos_λ[i, j], ord)
                            TaylorSeries.mul!(Rb2p[i, j, 1, 2], cos_ϕ[i, j], sin_λ[i, j], ord)
                            TaylorSeries.identity!(Rb2p[i, j, 2, 2], cos_λ[i, j], ord)
                            TaylorSeries.subst!(tmp2111[i, j], sin_ϕ[i, j], ord)
                            TaylorSeries.mul!(Rb2p[i, j, 3, 2], tmp2111[i, j], sin_λ[i, j], ord)
                            TaylorSeries.identity!(Rb2p[i, j, 1, 3], sin_ϕ[i, j], ord)
                            TaylorSeries.identity!(Rb2p[i, j, 2, 3], zero_q_1, ord)
                            TaylorSeries.identity!(Rb2p[i, j, 3, 3], cos_ϕ[i, j], ord)
                            TaylorSeries.mul!(tmp2113[i, j, 1, 1], Rb2p[i, j, 1, 1], M_[1, 1, j], ord)
                            TaylorSeries.mul!(tmp2114[i, j, 1, 2], Rb2p[i, j, 1, 2], M_[2, 1, j], ord)
                            TaylorSeries.add!(tmp2115[i, j, 1, 1], tmp2113[i, j, 1, 1], tmp2114[i, j, 1, 2], ord)
                            TaylorSeries.mul!(tmp2116[i, j, 1, 3], Rb2p[i, j, 1, 3], M_[3, 1, j], ord)
                            TaylorSeries.add!(Gc2p[i, j, 1, 1], tmp2115[i, j, 1, 1], tmp2116[i, j, 1, 3], ord)
                            TaylorSeries.mul!(tmp2118[i, j, 2, 1], Rb2p[i, j, 2, 1], M_[1, 1, j], ord)
                            TaylorSeries.mul!(tmp2119[i, j, 2, 2], Rb2p[i, j, 2, 2], M_[2, 1, j], ord)
                            TaylorSeries.add!(tmp2120[i, j, 2, 1], tmp2118[i, j, 2, 1], tmp2119[i, j, 2, 2], ord)
                            TaylorSeries.mul!(tmp2121[i, j, 2, 3], Rb2p[i, j, 2, 3], M_[3, 1, j], ord)
                            TaylorSeries.add!(Gc2p[i, j, 2, 1], tmp2120[i, j, 2, 1], tmp2121[i, j, 2, 3], ord)
                            TaylorSeries.mul!(tmp2123[i, j, 3, 1], Rb2p[i, j, 3, 1], M_[1, 1, j], ord)
                            TaylorSeries.mul!(tmp2124[i, j, 3, 2], Rb2p[i, j, 3, 2], M_[2, 1, j], ord)
                            TaylorSeries.add!(tmp2125[i, j, 3, 1], tmp2123[i, j, 3, 1], tmp2124[i, j, 3, 2], ord)
                            TaylorSeries.mul!(tmp2126[i, j, 3, 3], Rb2p[i, j, 3, 3], M_[3, 1, j], ord)
                            TaylorSeries.add!(Gc2p[i, j, 3, 1], tmp2125[i, j, 3, 1], tmp2126[i, j, 3, 3], ord)
                            TaylorSeries.mul!(tmp2128[i, j, 1, 1], Rb2p[i, j, 1, 1], M_[1, 2, j], ord)
                            TaylorSeries.mul!(tmp2129[i, j, 1, 2], Rb2p[i, j, 1, 2], M_[2, 2, j], ord)
                            TaylorSeries.add!(tmp2130[i, j, 1, 1], tmp2128[i, j, 1, 1], tmp2129[i, j, 1, 2], ord)
                            TaylorSeries.mul!(tmp2131[i, j, 1, 3], Rb2p[i, j, 1, 3], M_[3, 2, j], ord)
                            TaylorSeries.add!(Gc2p[i, j, 1, 2], tmp2130[i, j, 1, 1], tmp2131[i, j, 1, 3], ord)
                            TaylorSeries.mul!(tmp2133[i, j, 2, 1], Rb2p[i, j, 2, 1], M_[1, 2, j], ord)
                            TaylorSeries.mul!(tmp2134[i, j, 2, 2], Rb2p[i, j, 2, 2], M_[2, 2, j], ord)
                            TaylorSeries.add!(tmp2135[i, j, 2, 1], tmp2133[i, j, 2, 1], tmp2134[i, j, 2, 2], ord)
                            TaylorSeries.mul!(tmp2136[i, j, 2, 3], Rb2p[i, j, 2, 3], M_[3, 2, j], ord)
                            TaylorSeries.add!(Gc2p[i, j, 2, 2], tmp2135[i, j, 2, 1], tmp2136[i, j, 2, 3], ord)
                            TaylorSeries.mul!(tmp2138[i, j, 3, 1], Rb2p[i, j, 3, 1], M_[1, 2, j], ord)
                            TaylorSeries.mul!(tmp2139[i, j, 3, 2], Rb2p[i, j, 3, 2], M_[2, 2, j], ord)
                            TaylorSeries.add!(tmp2140[i, j, 3, 1], tmp2138[i, j, 3, 1], tmp2139[i, j, 3, 2], ord)
                            TaylorSeries.mul!(tmp2141[i, j, 3, 3], Rb2p[i, j, 3, 3], M_[3, 2, j], ord)
                            TaylorSeries.add!(Gc2p[i, j, 3, 2], tmp2140[i, j, 3, 1], tmp2141[i, j, 3, 3], ord)
                            TaylorSeries.mul!(tmp2143[i, j, 1, 1], Rb2p[i, j, 1, 1], M_[1, 3, j], ord)
                            TaylorSeries.mul!(tmp2144[i, j, 1, 2], Rb2p[i, j, 1, 2], M_[2, 3, j], ord)
                            TaylorSeries.add!(tmp2145[i, j, 1, 1], tmp2143[i, j, 1, 1], tmp2144[i, j, 1, 2], ord)
                            TaylorSeries.mul!(tmp2146[i, j, 1, 3], Rb2p[i, j, 1, 3], M_[3, 3, j], ord)
                            TaylorSeries.add!(Gc2p[i, j, 1, 3], tmp2145[i, j, 1, 1], tmp2146[i, j, 1, 3], ord)
                            TaylorSeries.mul!(tmp2148[i, j, 2, 1], Rb2p[i, j, 2, 1], M_[1, 3, j], ord)
                            TaylorSeries.mul!(tmp2149[i, j, 2, 2], Rb2p[i, j, 2, 2], M_[2, 3, j], ord)
                            TaylorSeries.add!(tmp2150[i, j, 2, 1], tmp2148[i, j, 2, 1], tmp2149[i, j, 2, 2], ord)
                            TaylorSeries.mul!(tmp2151[i, j, 2, 3], Rb2p[i, j, 2, 3], M_[3, 3, j], ord)
                            TaylorSeries.add!(Gc2p[i, j, 2, 3], tmp2150[i, j, 2, 1], tmp2151[i, j, 2, 3], ord)
                            TaylorSeries.mul!(tmp2153[i, j, 3, 1], Rb2p[i, j, 3, 1], M_[1, 3, j], ord)
                            TaylorSeries.mul!(tmp2154[i, j, 3, 2], Rb2p[i, j, 3, 2], M_[2, 3, j], ord)
                            TaylorSeries.add!(tmp2155[i, j, 3, 1], tmp2153[i, j, 3, 1], tmp2154[i, j, 3, 2], ord)
                            TaylorSeries.mul!(tmp2156[i, j, 3, 3], Rb2p[i, j, 3, 3], M_[3, 3, j], ord)
                            TaylorSeries.add!(Gc2p[i, j, 3, 3], tmp2155[i, j, 3, 1], tmp2156[i, j, 3, 3], ord)
                            TaylorSeries.mul!(tmp2158[i, j, 1, 1], F_JCS_ξ[i, j], Gc2p[i, j, 1, 1], ord)
                            TaylorSeries.mul!(tmp2159[i, j, 2, 1], F_JCS_η[i, j], Gc2p[i, j, 2, 1], ord)
                            TaylorSeries.add!(tmp2160[i, j, 1, 1], tmp2158[i, j, 1, 1], tmp2159[i, j, 2, 1], ord)
                            TaylorSeries.mul!(tmp2161[i, j, 3, 1], F_JCS_ζ[i, j], Gc2p[i, j, 3, 1], ord)
                            TaylorSeries.add!(F_JCS_x[i, j], tmp2160[i, j, 1, 1], tmp2161[i, j, 3, 1], ord)
                            TaylorSeries.mul!(tmp2163[i, j, 1, 2], F_JCS_ξ[i, j], Gc2p[i, j, 1, 2], ord)
                            TaylorSeries.mul!(tmp2164[i, j, 2, 2], F_JCS_η[i, j], Gc2p[i, j, 2, 2], ord)
                            TaylorSeries.add!(tmp2165[i, j, 1, 2], tmp2163[i, j, 1, 2], tmp2164[i, j, 2, 2], ord)
                            TaylorSeries.mul!(tmp2166[i, j, 3, 2], F_JCS_ζ[i, j], Gc2p[i, j, 3, 2], ord)
                            TaylorSeries.add!(F_JCS_y[i, j], tmp2165[i, j, 1, 2], tmp2166[i, j, 3, 2], ord)
                            TaylorSeries.mul!(tmp2168[i, j, 1, 3], F_JCS_ξ[i, j], Gc2p[i, j, 1, 3], ord)
                            TaylorSeries.mul!(tmp2169[i, j, 2, 3], F_JCS_η[i, j], Gc2p[i, j, 2, 3], ord)
                            TaylorSeries.add!(tmp2170[i, j, 1, 3], tmp2168[i, j, 1, 3], tmp2169[i, j, 2, 3], ord)
                            TaylorSeries.mul!(tmp2171[i, j, 3, 3], F_JCS_ζ[i, j], Gc2p[i, j, 3, 3], ord)
                            TaylorSeries.add!(F_JCS_z[i, j], tmp2170[i, j, 1, 3], tmp2171[i, j, 3, 3], ord)
                        end
                    end
                end
            end
        for j = 1:N_ext
            for i = 1:N_ext
                if i == j
                    continue
                else
                    if UJ_interaction[i, j]
                        TaylorSeries.mul!(tmp2173[i, j], μ[i], F_JCS_x[i, j], ord)
                        TaylorSeries.subst!(temp_accX_j[i, j], accX[j], tmp2173[i, j], ord)
                        TaylorSeries.identity!(accX[j], temp_accX_j[i, j], ord)
                        TaylorSeries.mul!(tmp2175[i, j], μ[i], F_JCS_y[i, j], ord)
                        TaylorSeries.subst!(temp_accY_j[i, j], accY[j], tmp2175[i, j], ord)
                        TaylorSeries.identity!(accY[j], temp_accY_j[i, j], ord)
                        TaylorSeries.mul!(tmp2177[i, j], μ[i], F_JCS_z[i, j], ord)
                        TaylorSeries.subst!(temp_accZ_j[i, j], accZ[j], tmp2177[i, j], ord)
                        TaylorSeries.identity!(accZ[j], temp_accZ_j[i, j], ord)
                        TaylorSeries.mul!(tmp2179[i, j], μ[j], F_JCS_x[i, j], ord)
                        TaylorSeries.add!(temp_accX_i[i, j], accX[i], tmp2179[i, j], ord)
                        TaylorSeries.identity!(accX[i], temp_accX_i[i, j], ord)
                        TaylorSeries.mul!(tmp2181[i, j], μ[j], F_JCS_y[i, j], ord)
                        TaylorSeries.add!(temp_accY_i[i, j], accY[i], tmp2181[i, j], ord)
                        TaylorSeries.identity!(accY[i], temp_accY_i[i, j], ord)
                        TaylorSeries.mul!(tmp2183[i, j], μ[j], F_JCS_z[i, j], ord)
                        TaylorSeries.add!(temp_accZ_i[i, j], accZ[i], tmp2183[i, j], ord)
                        TaylorSeries.identity!(accZ[i], temp_accZ_i[i, j], ord)
                    end
                end
            end
        end
        #= REPL[4]:439 =# Threads.@threads for j = 1:N
                for i = 1:N
                    if i == j
                        continue
                    else
                        TaylorSeries.mul!(_4ϕj[i, j], 4, newtonianNb_Potential[j], ord)
                        TaylorSeries.add!(ϕi_plus_4ϕj[i, j], newtonianNb_Potential[i], _4ϕj[i, j], ord)
                        TaylorSeries.mul!(tmp2189[i], 2, v2[i], ord)
                        TaylorSeries.add!(tmp2190[j], v2[j], tmp2189[i], ord)
                        TaylorSeries.mul!(tmp2192[i, j], 4, vi_dot_vj[i, j], ord)
                        TaylorSeries.subst!(sj2_plus_2si2_minus_4vivj[i, j], tmp2190[j], tmp2192[i, j], ord)
                        TaylorSeries.subst!(ϕs_and_vs[i, j], sj2_plus_2si2_minus_4vivj[i, j], ϕi_plus_4ϕj[i, j], ord)
                        TaylorSeries.mul!(Xij_t_Ui[i, j], X[i, j], dq[3i - 2], ord)
                        TaylorSeries.mul!(Yij_t_Vi[i, j], Y[i, j], dq[3i - 1], ord)
                        TaylorSeries.mul!(Zij_t_Wi[i, j], Z[i, j], dq[3i], ord)
                        TaylorSeries.add!(tmp2198[i, j], Xij_t_Ui[i, j], Yij_t_Vi[i, j], ord)
                        TaylorSeries.add!(Rij_dot_Vi[i, j], tmp2198[i, j], Zij_t_Wi[i, j], ord)
                        TaylorSeries.pow!(tmp2201[i, j], Rij_dot_Vi[i, j], 2, ord)
                        TaylorSeries.div!(pn1t7[i, j], tmp2201[i, j], r_p2[i, j], ord)
                        TaylorSeries.mul!(tmp2204[i, j], 1.5, pn1t7[i, j], ord)
                        TaylorSeries.subst!(pn1t2_7[i, j], ϕs_and_vs[i, j], tmp2204[i, j], ord)
                        TaylorSeries.add!(pn1t1_7[i, j], c_p2, pn1t2_7[i, j], ord)
                        for k = 1:postnewton_iter
                            TaylorSeries.identity!(pn1[i, j, k], zero_q_1, ord)
                            TaylorSeries.identity!(X_t_pn1[i, j, k], zero_q_1, ord)
                            TaylorSeries.identity!(Y_t_pn1[i, j, k], zero_q_1, ord)
                            TaylorSeries.identity!(Z_t_pn1[i, j, k], zero_q_1, ord)
                            TaylorSeries.identity!(pNX_t_pn3[i, j, k], zero_q_1, ord)
                            TaylorSeries.identity!(pNY_t_pn3[i, j, k], zero_q_1, ord)
                            TaylorSeries.identity!(pNZ_t_pn3[i, j, k], zero_q_1, ord)
                            TaylorSeries.identity!(pNX_t_X[i, j, k], zero_q_1, ord)
                            TaylorSeries.identity!(pNY_t_Y[i, j, k], zero_q_1, ord)
                            TaylorSeries.identity!(pNZ_t_Z[i, j, k], zero_q_1, ord)
                        end
                    end
                end
                TaylorSeries.identity!(postNewtonX[j, 1], newtonX[j], ord)
                TaylorSeries.identity!(postNewtonY[j, 1], newtonY[j], ord)
                TaylorSeries.identity!(postNewtonZ[j, 1], newtonZ[j], ord)
                for k = 1:postnewton_iter
                    TaylorSeries.identity!(pntempX[j, k], zero_q_1, ord)
                    TaylorSeries.identity!(pntempY[j, k], zero_q_1, ord)
                    TaylorSeries.identity!(pntempZ[j, k], zero_q_1, ord)
                end
            end
        for k = 1:postnewton_iter
            #= REPL[4]:484 =# Threads.@threads for j = 1:N
                    for i = 1:N
                        if i == j
                            continue
                        else
                            TaylorSeries.mul!(pNX_t_X[i, j, k], postNewtonX[i, k], X[i, j], ord)
                            TaylorSeries.mul!(pNY_t_Y[i, j, k], postNewtonY[i, k], Y[i, j], ord)
                            TaylorSeries.mul!(pNZ_t_Z[i, j, k], postNewtonZ[i, k], Z[i, j], ord)
                            TaylorSeries.add!(tmp2211[i, j, k], pNX_t_X[i, j, k], pNY_t_Y[i, j, k], ord)
                            TaylorSeries.add!(tmp2212[i, j, k], tmp2211[i, j, k], pNZ_t_Z[i, j, k], ord)
                            TaylorSeries.mul!(tmp2213[i, j, k], 0.5, tmp2212[i, j, k], ord)
                            TaylorSeries.add!(pn1[i, j, k], pn1t1_7[i, j], tmp2213[i, j, k], ord)
                            TaylorSeries.mul!(X_t_pn1[i, j, k], newton_acc_X[i, j], pn1[i, j, k], ord)
                            TaylorSeries.mul!(Y_t_pn1[i, j, k], newton_acc_Y[i, j], pn1[i, j, k], ord)
                            TaylorSeries.mul!(Z_t_pn1[i, j, k], newton_acc_Z[i, j], pn1[i, j, k], ord)
                            TaylorSeries.mul!(pNX_t_pn3[i, j, k], postNewtonX[i, k], pn3[i, j], ord)
                            TaylorSeries.mul!(pNY_t_pn3[i, j, k], postNewtonY[i, k], pn3[i, j], ord)
                            TaylorSeries.mul!(pNZ_t_pn3[i, j, k], postNewtonZ[i, k], pn3[i, j], ord)
                            TaylorSeries.add!(tmp2221[i, j, k], U_t_pn2[i, j], pNX_t_pn3[i, j, k], ord)
                            TaylorSeries.add!(termpnx[i, j, k], X_t_pn1[i, j, k], tmp2221[i, j, k], ord)
                            TaylorSeries.add!(sumpnx[i, j, k], pntempX[j, k], termpnx[i, j, k], ord)
                            TaylorSeries.identity!(pntempX[j, k], sumpnx[i, j, k], ord)
                            TaylorSeries.add!(tmp2224[i, j, k], V_t_pn2[i, j], pNY_t_pn3[i, j, k], ord)
                            TaylorSeries.add!(termpny[i, j, k], Y_t_pn1[i, j, k], tmp2224[i, j, k], ord)
                            TaylorSeries.add!(sumpny[i, j, k], pntempY[j, k], termpny[i, j, k], ord)
                            TaylorSeries.identity!(pntempY[j, k], sumpny[i, j, k], ord)
                            TaylorSeries.add!(tmp2227[i, j, k], W_t_pn2[i, j], pNZ_t_pn3[i, j, k], ord)
                            TaylorSeries.add!(termpnz[i, j, k], Z_t_pn1[i, j, k], tmp2227[i, j, k], ord)
                            TaylorSeries.add!(sumpnz[i, j, k], pntempZ[j, k], termpnz[i, j, k], ord)
                            TaylorSeries.identity!(pntempZ[j, k], sumpnz[i, j, k], ord)
                        end
                    end
                    TaylorSeries.mul!(postNewtonX[j, k + 1], pntempX[j, k], c_m2, ord)
                    TaylorSeries.mul!(postNewtonY[j, k + 1], pntempY[j, k], c_m2, ord)
                    TaylorSeries.mul!(postNewtonZ[j, k + 1], pntempZ[j, k], c_m2, ord)
                end
        end
        #= REPL[4]:521 =# Threads.@threads for i = 1:N_ext
                TaylorSeries.add!(dq[3 * (N + i) - 2], postNewtonX[i, postnewton_iter + 1], accX[i], ord)
                TaylorSeries.add!(dq[3 * (N + i) - 1], postNewtonY[i, postnewton_iter + 1], accY[i], ord)
                TaylorSeries.add!(dq[3 * (N + i)], postNewtonZ[i, postnewton_iter + 1], accZ[i], ord)
            end
        #= REPL[4]:526 =# Threads.@threads for i = N_ext + 1:N
                TaylorSeries.identity!(dq[3 * (N + i) - 2], postNewtonX[i, postnewton_iter + 1], ord)
                TaylorSeries.identity!(dq[3 * (N + i) - 1], postNewtonY[i, postnewton_iter + 1], ord)
                TaylorSeries.identity!(dq[3 * (N + i)], postNewtonZ[i, postnewton_iter + 1], ord)
            end
        for __idx = eachindex(q)
            (q[__idx]).coeffs[ordnext + 1] = (dq[__idx]).coeffs[ordnext] / ordnext
        end
    end
    return nothing
end

function TaylorIntegration.jetcoeffs!(::Val{DE430!}, t::Taylor1{_T}, q::AbstractVector{Taylor1{_S}}, dq::AbstractVector{Taylor1{_S}}, params) where {_T <: Real, _S <: Number}
    order = t.order
    local (N, S, eulang_de430_, jd0) = params
    local N_ext = 11
    local N_back = 11
    local params_back = (N_back, S, eulang_de430_, jd0)
    local qq_ = Taylor1.(constant_term.(q[union(1:3N_back, 3N + 1:3N + 3N_back)]), t.order)
    local dqq_ = similar(qq_)
    local jtcffs = TaylorIntegration.jetcoeffs!(Val(NBP_pN_A_J23E_J23M_J2S_threads!), t, qq_, dqq_, params_back)
    local __t = Taylor1(t.order)
    local q_del_τ_M = qq_(__t - τ_M)
    local q_del_τ_0 = qq_(__t - τ_0p)
    local q_del_τ_1 = qq_(__t - τ_1p)
    local q_del_τ_2 = qq_(__t - τ_2p)
    local dsj2k = t + (jd0 - 2.451545e6)
    local eulang_t = eulang_de430_(dsj2k * daysec)
    local eulang_t_del = eulang_de430_((dsj2k - τ_M) * daysec)
    local postnewton_iter = 1
    local c_p2 = 29979.063823897606
    local c_m2 = 3.3356611996764786e-5
    local zero_q_1 = zero(q[1])
    local one_t = one(t)
    X = Array{Taylor1{S}}(undef, N, N)
    Y = Array{Taylor1{S}}(undef, N, N)
    Z = Array{Taylor1{S}}(undef, N, N)
    r_p2 = Array{Taylor1{S}}(undef, N, N)
    r_p3d2 = Array{Taylor1{S}}(undef, N, N)
    r_p7d2 = Array{Taylor1{S}}(undef, N, N)
    newtonX = Array{Taylor1{S}}(undef, N)
    newtonY = Array{Taylor1{S}}(undef, N)
    newtonZ = Array{Taylor1{S}}(undef, N)
    newtonianCoeff = Array{Taylor1{S}}(undef, N, N)
    U = Array{Taylor1{S}}(undef, N, N)
    V = Array{Taylor1{S}}(undef, N, N)
    W = Array{Taylor1{S}}(undef, N, N)
    _4U_m_3X = Array{Taylor1{S}}(undef, N, N)
    _4V_m_3Y = Array{Taylor1{S}}(undef, N, N)
    _4W_m_3Z = Array{Taylor1{S}}(undef, N, N)
    UU = Array{Taylor1{S}}(undef, N, N)
    VV = Array{Taylor1{S}}(undef, N, N)
    WW = Array{Taylor1{S}}(undef, N, N)
    r_p1d2 = Array{Taylor1{S}}(undef, N, N)
    newtonianNb_Potential = Array{Taylor1{S}}(undef, N)
    newtonian1b_Potential = Array{Taylor1{S}}(undef, N, N)
    newton_acc_X = Array{Taylor1{S}}(undef, N, N)
    newton_acc_Y = Array{Taylor1{S}}(undef, N, N)
    newton_acc_Z = Array{Taylor1{S}}(undef, N, N)
    v2 = Array{Taylor1{S}}(undef, N)
    vi_dot_vj = Array{Taylor1{S}}(undef, N, N)
    pn2 = Array{Taylor1{S}}(undef, N, N)
    pn3 = Array{Taylor1{S}}(undef, N, N)
    _4ϕj = Array{Taylor1{S}}(undef, N, N)
    ϕi_plus_4ϕj = Array{Taylor1{S}}(undef, N, N)
    sj2_plus_2si2_minus_4vivj = Array{Taylor1{S}}(undef, N, N)
    ϕs_and_vs = Array{Taylor1{S}}(undef, N, N)
    U_t_pn2 = Array{Taylor1{S}}(undef, N, N)
    V_t_pn2 = Array{Taylor1{S}}(undef, N, N)
    W_t_pn2 = Array{Taylor1{S}}(undef, N, N)
    pn1t1_7 = Array{Taylor1{S}}(undef, N, N)
    pntempX = Array{Taylor1{S}}(undef, N, postnewton_iter)
    pntempY = Array{Taylor1{S}}(undef, N, postnewton_iter)
    pntempZ = Array{Taylor1{S}}(undef, N, postnewton_iter)
    pn1 = Array{Taylor1{S}}(undef, N, N, postnewton_iter)
    X_t_pn1 = Array{Taylor1{S}}(undef, N, N, postnewton_iter)
    Y_t_pn1 = Array{Taylor1{S}}(undef, N, N, postnewton_iter)
    Z_t_pn1 = Array{Taylor1{S}}(undef, N, N, postnewton_iter)
    pNX_t_pn3 = Array{Taylor1{S}}(undef, N, N, postnewton_iter)
    pNY_t_pn3 = Array{Taylor1{S}}(undef, N, N, postnewton_iter)
    pNZ_t_pn3 = Array{Taylor1{S}}(undef, N, N, postnewton_iter)
    pNX_t_X = Array{Taylor1{S}}(undef, N, N, postnewton_iter)
    pNY_t_Y = Array{Taylor1{S}}(undef, N, N, postnewton_iter)
    pNZ_t_Z = Array{Taylor1{S}}(undef, N, N, postnewton_iter)
    postNewtonX = Array{Taylor1{S}}(undef, N, postnewton_iter + 1)
    postNewtonY = Array{Taylor1{S}}(undef, N, postnewton_iter + 1)
    postNewtonZ = Array{Taylor1{S}}(undef, N, postnewton_iter + 1)
    X_bf_1 = Array{Taylor1{S}}(undef, N_ext, N_ext)
    Y_bf_1 = Array{Taylor1{S}}(undef, N_ext, N_ext)
    Z_bf_1 = Array{Taylor1{S}}(undef, N_ext, N_ext)
    X_bf_2 = Array{Taylor1{S}}(undef, N_ext, N_ext)
    Y_bf_2 = Array{Taylor1{S}}(undef, N_ext, N_ext)
    Z_bf_2 = Array{Taylor1{S}}(undef, N_ext, N_ext)
    X_bf_3 = Array{Taylor1{S}}(undef, N_ext, N_ext)
    Y_bf_3 = Array{Taylor1{S}}(undef, N_ext, N_ext)
    Z_bf_3 = Array{Taylor1{S}}(undef, N_ext, N_ext)
    X_bf = Array{Taylor1{S}}(undef, N_ext, N_ext)
    Y_bf = Array{Taylor1{S}}(undef, N_ext, N_ext)
    Z_bf = Array{Taylor1{S}}(undef, N_ext, N_ext)
    F_JCS_x = Array{Taylor1{S}}(undef, N_ext, N_ext)
    F_JCS_y = Array{Taylor1{S}}(undef, N_ext, N_ext)
    F_JCS_z = Array{Taylor1{S}}(undef, N_ext, N_ext)
    temp_accX_j = Array{Taylor1{S}}(undef, N_ext, N_ext)
    temp_accY_j = Array{Taylor1{S}}(undef, N_ext, N_ext)
    temp_accZ_j = Array{Taylor1{S}}(undef, N_ext, N_ext)
    temp_accX_i = Array{Taylor1{S}}(undef, N_ext, N_ext)
    temp_accY_i = Array{Taylor1{S}}(undef, N_ext, N_ext)
    temp_accZ_i = Array{Taylor1{S}}(undef, N_ext, N_ext)
    sin_ϕ = Array{Taylor1{S}}(undef, N_ext, N_ext)
    cos_ϕ = Array{Taylor1{S}}(undef, N_ext, N_ext)
    sin_λ = Array{Taylor1{S}}(undef, N_ext, N_ext)
    cos_λ = Array{Taylor1{S}}(undef, N_ext, N_ext)
    r_xy = Array{Taylor1{S}}(undef, N_ext, N_ext)
    r_p4 = Array{Taylor1{S}}(undef, N_ext, N_ext)
    P_n = Array{Taylor1{S}}(undef, N_ext, N_ext, maximum(n1SEM) + 1)
    dP_n = Array{Taylor1{S}}(undef, N_ext, N_ext, maximum(n1SEM) + 1)
    temp_fjξ = Array{Taylor1{S}}(undef, N_ext, N_ext, maximum(n1SEM) + 1)
    temp_fjζ = Array{Taylor1{S}}(undef, N_ext, N_ext, maximum(n1SEM) + 1)
    temp_rn = Array{Taylor1{S}}(undef, N_ext, N_ext, maximum(n1SEM) + 1)
    F_CS_ξ_36 = Array{Taylor1{S}}(undef, N_ext, N_ext)
    F_CS_η_36 = Array{Taylor1{S}}(undef, N_ext, N_ext)
    F_CS_ζ_36 = Array{Taylor1{S}}(undef, N_ext, N_ext)
    F_J_ξ_36 = Array{Taylor1{S}}(undef, N_ext, N_ext)
    F_J_ζ_36 = Array{Taylor1{S}}(undef, N_ext, N_ext)
    sin_mλ = Array{Taylor1{S}}(undef, N_ext, N_ext, n1SEM[mo])
    cos_mλ = Array{Taylor1{S}}(undef, N_ext, N_ext, n1SEM[mo])
    secϕ_P_nm = Array{Taylor1{S}}(undef, N_ext, N_ext, n1SEM[mo] + 1, n1SEM[mo] + 1)
    P_nm = Array{Taylor1{S}}(undef, N_ext, N_ext, n1SEM[mo] + 1, n1SEM[mo] + 1)
    cosϕ_dP_nm = Array{Taylor1{S}}(undef, N_ext, N_ext, n1SEM[mo] + 1, n1SEM[mo] + 1)
    F_J_ξ = Array{Taylor1{S}}(undef, N_ext, N_ext)
    F_J_η = Array{Taylor1{S}}(undef, N_ext, N_ext)
    F_J_ζ = Array{Taylor1{S}}(undef, N_ext, N_ext)
    F_CS_ξ = Array{Taylor1{S}}(undef, N_ext, N_ext)
    F_CS_η = Array{Taylor1{S}}(undef, N_ext, N_ext)
    F_CS_ζ = Array{Taylor1{S}}(undef, N_ext, N_ext)
    F_JCS_ξ = Array{Taylor1{S}}(undef, N_ext, N_ext)
    F_JCS_η = Array{Taylor1{S}}(undef, N_ext, N_ext)
    F_JCS_ζ = Array{Taylor1{S}}(undef, N_ext, N_ext)
    Rb2p = Array{Taylor1{S}}(undef, N_ext, N_ext, 3, 3)
    Gc2p = Array{Taylor1{S}}(undef, N_ext, N_ext, 3, 3)
    accX = Array{Taylor1{S}}(undef, N_ext)
    accY = Array{Taylor1{S}}(undef, N_ext)
    accZ = Array{Taylor1{S}}(undef, N_ext)
    r_star_M_0 = Array{Taylor1{S}}(undef, 3)
    r_star_S_0 = Array{Taylor1{S}}(undef, 3)
    r_star_M_1 = Array{Taylor1{S}}(undef, 3)
    r_star_S_1 = Array{Taylor1{S}}(undef, 3)
    r_star_M_2 = Array{Taylor1{S}}(undef, 3)
    r_star_S_2 = Array{Taylor1{S}}(undef, 3)
    local αs = deg2rad(α_p_sun * one_t)
    local δs = deg2rad(δ_p_sun * one_t)
    local αm = eulang_t[1] - pi / 2
    local δm = pi / 2 - eulang_t[2]
    local Wm = eulang_t[3]
    local M_ = Array{Taylor1{S}}(undef, 3, 3, 5)
    local M_[:, :, ea] = c2t_jpl_de430(dsj2k)
    local M_[:, :, su] = pole_rotation(αs, δs)
    local M_[:, :, mo] = pole_rotation(αm, δm, Wm)
    local M_del_mo = pole_rotation(eulang_t_del[1] - pi / 2, pi / 2 - eulang_t_del[2], eulang_t_del[3])
    ITM_t = Array{Taylor1{S}}(undef, 3, 3)
    ITM2_t = Array{Taylor1{S}}(undef, 3, 3)
    local ITM2_t = ITM_und .* one_t + ITM2(eulang_t_del[1], eulang_t_del[2], eulang_t_del[3])
    local fact_num = -4.5257273867882326e-36
    local fact1_jsem = [(2n - 1) / n for n = 1:maximum(n1SEM)]
    local fact2_jsem = [(n - 1) / n for n = 1:maximum(n1SEM)]
    local fact3_jsem = [n for n = 1:maximum(n1SEM)]
    local fact4_jsem = [n + 1 for n = 1:maximum(n1SEM)]
    local fact5_jsem = [n + 2 for n = 1:maximum(n1SEM)]
    local lnm1 = [(2n - 1) / (n - m) for n = 1:6, m = 1:6]
    local lnm2 = [-(((n + m) - 1)) / (n - m) for n = 1:6, m = 1:6]
    local lnm3 = [-n for n = 1:6]
    local lnm4 = [n + m for n = 1:6, m = 1:6]
    local lnm5 = [2n - 1 for n = 1:6]
    local lnm6 = [-((n + 1)) for n = 1:6]
    local lnm7 = [m for m = 1:6]
    local RE_au = RE / au
    local J2E_t = (J2E + J2EDOT * (dsj2k / yr)) * RE_au ^ 2
    local J2S_t = JSEM[su, 2] * one_t
    local J2_t = Array{Taylor1{S}}(undef, 5)
    J2_t[su] = Taylor1(identity(constant_term(J2S_t)), order)
    J2_t[ea] = Taylor1(identity(constant_term(J2E_t)), order)
    local R30 = c2t_jpl_de430(dsj2k - τ_0p)
    local R31 = Rz(-ω_E * τ_1) * c2t_jpl_de430(dsj2k - τ_1p)
    local R32 = Rz(-ω_E * τ_2) * c2t_jpl_de430(dsj2k - τ_2p)
    local tid_num_coeff = 1.5 * (1.0 + μ[mo] / μ[ea])
    #= REPL[2]:211 =# Threads.@threads for j = 1:N
            newtonX[j] = Taylor1(identity(constant_term(zero_q_1)), order)
            newtonY[j] = Taylor1(identity(constant_term(zero_q_1)), order)
            newtonZ[j] = Taylor1(identity(constant_term(zero_q_1)), order)
            newtonianNb_Potential[j] = Taylor1(identity(constant_term(zero_q_1)), order)
            dq[3j - 2] = Taylor1(identity(constant_term(q[3 * (N + j) - 2])), order)
            dq[3j - 1] = Taylor1(identity(constant_term(q[3 * (N + j) - 1])), order)
            dq[3j] = Taylor1(identity(constant_term(q[3 * (N + j)])), order)
        end
    #= REPL[2]:221 =# Threads.@threads for j = 1:N_ext
            accX[j] = Taylor1(identity(constant_term(zero_q_1)), order)
            accY[j] = Taylor1(identity(constant_term(zero_q_1)), order)
            accZ[j] = Taylor1(identity(constant_term(zero_q_1)), order)
        end
    tmp1128 = Array{Taylor1{_S}}(undef, size(dq))
    tmp1128 .= Taylor1(zero(_S), order)
    tmp1130 = Array{Taylor1{_S}}(undef, size(dq))
    tmp1130 .= Taylor1(zero(_S), order)
    tmp1133 = Array{Taylor1{_S}}(undef, size(dq))
    tmp1133 .= Taylor1(zero(_S), order)
    tmp1135 = Array{Taylor1{_S}}(undef, size(dq))
    tmp1135 .= Taylor1(zero(_S), order)
    tmp1138 = Array{Taylor1{_S}}(undef, size(dq))
    tmp1138 .= Taylor1(zero(_S), order)
    tmp1140 = Array{Taylor1{_S}}(undef, size(dq))
    tmp1140 .= Taylor1(zero(_S), order)
    pn2x = Array{Taylor1{_S}}(undef, size(X))
    pn2x .= Taylor1(zero(_S), order)
    pn2y = Array{Taylor1{_S}}(undef, size(Y))
    pn2y .= Taylor1(zero(_S), order)
    pn2z = Array{Taylor1{_S}}(undef, size(Z))
    pn2z .= Taylor1(zero(_S), order)
    tmp1148 = Array{Taylor1{_S}}(undef, size(UU))
    tmp1148 .= Taylor1(zero(_S), order)
    tmp1151 = Array{Taylor1{_S}}(undef, size(X))
    tmp1151 .= Taylor1(zero(_S), order)
    tmp1153 = Array{Taylor1{_S}}(undef, size(Y))
    tmp1153 .= Taylor1(zero(_S), order)
    tmp1154 = Array{Taylor1{_S}}(undef, size(tmp1151))
    tmp1154 .= Taylor1(zero(_S), order)
    tmp1156 = Array{Taylor1{_S}}(undef, size(Z))
    tmp1156 .= Taylor1(zero(_S), order)
    tmp1164 = Array{Taylor1{_S}}(undef, size(pn2x))
    tmp1164 .= Taylor1(zero(_S), order)
    tmp1165 = Array{Taylor1{_S}}(undef, size(tmp1164))
    tmp1165 .= Taylor1(zero(_S), order)
    tmp1176 = Array{Taylor1{_S}}(undef, size(X))
    tmp1176 .= Taylor1(zero(_S), order)
    temp_001 = Array{Taylor1{_S}}(undef, size(tmp1176))
    temp_001 .= Taylor1(zero(_S), order)
    tmp1178 = Array{Taylor1{_S}}(undef, size(Y))
    tmp1178 .= Taylor1(zero(_S), order)
    temp_002 = Array{Taylor1{_S}}(undef, size(tmp1178))
    temp_002 .= Taylor1(zero(_S), order)
    tmp1180 = Array{Taylor1{_S}}(undef, size(Z))
    tmp1180 .= Taylor1(zero(_S), order)
    temp_003 = Array{Taylor1{_S}}(undef, size(tmp1180))
    temp_003 .= Taylor1(zero(_S), order)
    temp_004 = Array{Taylor1{_S}}(undef, size(newtonian1b_Potential))
    temp_004 .= Taylor1(zero(_S), order)
    tmp1184 = Array{Taylor1{_S}}(undef, size(dq))
    tmp1184 .= Taylor1(zero(_S), order)
    tmp1186 = Array{Taylor1{_S}}(undef, size(dq))
    tmp1186 .= Taylor1(zero(_S), order)
    tmp1187 = Array{Taylor1{_S}}(undef, size(tmp1184))
    tmp1187 .= Taylor1(zero(_S), order)
    tmp1189 = Array{Taylor1{_S}}(undef, size(dq))
    tmp1189 .= Taylor1(zero(_S), order)
    #= REPL[2]:228 =# Threads.@threads for j = 1:N
            for i = 1:N
                if i == j
                    continue
                else
                    X[i, j] = Taylor1(constant_term(q[3i - 2]) - constant_term(q[3j - 2]), order)
                    Y[i, j] = Taylor1(constant_term(q[3i - 1]) - constant_term(q[3j - 1]), order)
                    Z[i, j] = Taylor1(constant_term(q[3i]) - constant_term(q[3j]), order)
                    U[i, j] = Taylor1(constant_term(dq[3i - 2]) - constant_term(dq[3j - 2]), order)
                    V[i, j] = Taylor1(constant_term(dq[3i - 1]) - constant_term(dq[3j - 1]), order)
                    W[i, j] = Taylor1(constant_term(dq[3i]) - constant_term(dq[3j]), order)
                    tmp1128[3j - 2] = Taylor1(constant_term(4) * constant_term(dq[3j - 2]), order)
                    tmp1130[3i - 2] = Taylor1(constant_term(3) * constant_term(dq[3i - 2]), order)
                    _4U_m_3X[i, j] = Taylor1(constant_term(tmp1128[3j - 2]) - constant_term(tmp1130[3i - 2]), order)
                    tmp1133[3j - 1] = Taylor1(constant_term(4) * constant_term(dq[3j - 1]), order)
                    tmp1135[3i - 1] = Taylor1(constant_term(3) * constant_term(dq[3i - 1]), order)
                    _4V_m_3Y[i, j] = Taylor1(constant_term(tmp1133[3j - 1]) - constant_term(tmp1135[3i - 1]), order)
                    tmp1138[3j] = Taylor1(constant_term(4) * constant_term(dq[3j]), order)
                    tmp1140[3i] = Taylor1(constant_term(3) * constant_term(dq[3i]), order)
                    _4W_m_3Z[i, j] = Taylor1(constant_term(tmp1138[3j]) - constant_term(tmp1140[3i]), order)
                    pn2x[i, j] = Taylor1(constant_term(X[i, j]) * constant_term(_4U_m_3X[i, j]), order)
                    pn2y[i, j] = Taylor1(constant_term(Y[i, j]) * constant_term(_4V_m_3Y[i, j]), order)
                    pn2z[i, j] = Taylor1(constant_term(Z[i, j]) * constant_term(_4W_m_3Z[i, j]), order)
                    UU[i, j] = Taylor1(constant_term(dq[3i - 2]) * constant_term(dq[3j - 2]), order)
                    VV[i, j] = Taylor1(constant_term(dq[3i - 1]) * constant_term(dq[3j - 1]), order)
                    WW[i, j] = Taylor1(constant_term(dq[3i]) * constant_term(dq[3j]), order)
                    tmp1148[i, j] = Taylor1(constant_term(UU[i, j]) + constant_term(VV[i, j]), order)
                    vi_dot_vj[i, j] = Taylor1(constant_term(tmp1148[i, j]) + constant_term(WW[i, j]), order)
                    tmp1151[i, j] = Taylor1(constant_term(X[i, j]) ^ constant_term(2), order)
                    tmp1153[i, j] = Taylor1(constant_term(Y[i, j]) ^ constant_term(2), order)
                    tmp1154[i, j] = Taylor1(constant_term(tmp1151[i, j]) + constant_term(tmp1153[i, j]), order)
                    tmp1156[i, j] = Taylor1(constant_term(Z[i, j]) ^ constant_term(2), order)
                    r_p2[i, j] = Taylor1(constant_term(tmp1154[i, j]) + constant_term(tmp1156[i, j]), order)
                    r_p1d2[i, j] = Taylor1(sqrt(constant_term(r_p2[i, j])), order)
                    r_p3d2[i, j] = Taylor1(constant_term(r_p2[i, j]) ^ constant_term(1.5), order)
                    r_p7d2[i, j] = Taylor1(constant_term(r_p2[i, j]) ^ constant_term(3.5), order)
                    newtonianCoeff[i, j] = Taylor1(constant_term(μ[i]) / constant_term(r_p3d2[i, j]), order)
                    tmp1164[i, j] = Taylor1(constant_term(pn2x[i, j]) + constant_term(pn2y[i, j]), order)
                    tmp1165[i, j] = Taylor1(constant_term(tmp1164[i, j]) + constant_term(pn2z[i, j]), order)
                    pn2[i, j] = Taylor1(constant_term(newtonianCoeff[i, j]) * constant_term(tmp1165[i, j]), order)
                    newton_acc_X[i, j] = Taylor1(constant_term(X[i, j]) * constant_term(newtonianCoeff[i, j]), order)
                    newton_acc_Y[i, j] = Taylor1(constant_term(Y[i, j]) * constant_term(newtonianCoeff[i, j]), order)
                    newton_acc_Z[i, j] = Taylor1(constant_term(Z[i, j]) * constant_term(newtonianCoeff[i, j]), order)
                    newtonian1b_Potential[i, j] = Taylor1(constant_term(μ[i]) / constant_term(r_p1d2[i, j]), order)
                    pn3[i, j] = Taylor1(constant_term(3.5) * constant_term(newtonian1b_Potential[i, j]), order)
                    U_t_pn2[i, j] = Taylor1(constant_term(pn2[i, j]) * constant_term(U[i, j]), order)
                    V_t_pn2[i, j] = Taylor1(constant_term(pn2[i, j]) * constant_term(V[i, j]), order)
                    W_t_pn2[i, j] = Taylor1(constant_term(pn2[i, j]) * constant_term(W[i, j]), order)
                    tmp1176[i, j] = Taylor1(constant_term(X[i, j]) * constant_term(newtonianCoeff[i, j]), order)
                    temp_001[i, j] = Taylor1(constant_term(newtonX[j]) + constant_term(tmp1176[i, j]), order)
                    newtonX[j] = Taylor1(identity(constant_term(temp_001[i, j])), order)
                    tmp1178[i, j] = Taylor1(constant_term(Y[i, j]) * constant_term(newtonianCoeff[i, j]), order)
                    temp_002[i, j] = Taylor1(constant_term(newtonY[j]) + constant_term(tmp1178[i, j]), order)
                    newtonY[j] = Taylor1(identity(constant_term(temp_002[i, j])), order)
                    tmp1180[i, j] = Taylor1(constant_term(Z[i, j]) * constant_term(newtonianCoeff[i, j]), order)
                    temp_003[i, j] = Taylor1(constant_term(newtonZ[j]) + constant_term(tmp1180[i, j]), order)
                    newtonZ[j] = Taylor1(identity(constant_term(temp_003[i, j])), order)
                    temp_004[i, j] = Taylor1(constant_term(newtonianNb_Potential[j]) + constant_term(newtonian1b_Potential[i, j]), order)
                    newtonianNb_Potential[j] = Taylor1(identity(constant_term(temp_004[i, j])), order)
                end
            end
            tmp1184[3j - 2] = Taylor1(constant_term(dq[3j - 2]) ^ constant_term(2), order)
            tmp1186[3j - 1] = Taylor1(constant_term(dq[3j - 1]) ^ constant_term(2), order)
            tmp1187[3j - 2] = Taylor1(constant_term(tmp1184[3j - 2]) + constant_term(tmp1186[3j - 1]), order)
            tmp1189[3j] = Taylor1(constant_term(dq[3j]) ^ constant_term(2), order)
            v2[j] = Taylor1(constant_term(tmp1187[3j - 2]) + constant_term(tmp1189[3j]), order)
        end
    X_me_del_τ_M = Taylor1(constant_term(q_del_τ_M[3mo - 2]) - constant_term(q_del_τ_M[3ea - 2]), order)
    Y_me_del_τ_M = Taylor1(constant_term(q_del_τ_M[3mo - 1]) - constant_term(q_del_τ_M[3ea - 1]), order)
    Z_me_del_τ_M = Taylor1(constant_term(q_del_τ_M[3mo]) - constant_term(q_del_τ_M[3ea]), order)
    tmp1194 = Taylor1(constant_term(M_del_mo[1, 1]) * constant_term(X_me_del_τ_M), order)
    tmp1195 = Taylor1(constant_term(M_del_mo[1, 2]) * constant_term(Y_me_del_τ_M), order)
    tmp1196 = Taylor1(constant_term(tmp1194) + constant_term(tmp1195), order)
    tmp1197 = Taylor1(constant_term(M_del_mo[1, 3]) * constant_term(Z_me_del_τ_M), order)
    xmed = Taylor1(constant_term(tmp1196) + constant_term(tmp1197), order)
    tmp1199 = Taylor1(constant_term(M_del_mo[2, 1]) * constant_term(X_me_del_τ_M), order)
    tmp1200 = Taylor1(constant_term(M_del_mo[2, 2]) * constant_term(Y_me_del_τ_M), order)
    tmp1201 = Taylor1(constant_term(tmp1199) + constant_term(tmp1200), order)
    tmp1202 = Taylor1(constant_term(M_del_mo[2, 3]) * constant_term(Z_me_del_τ_M), order)
    ymed = Taylor1(constant_term(tmp1201) + constant_term(tmp1202), order)
    tmp1204 = Taylor1(constant_term(M_del_mo[3, 1]) * constant_term(X_me_del_τ_M), order)
    tmp1205 = Taylor1(constant_term(M_del_mo[3, 2]) * constant_term(Y_me_del_τ_M), order)
    tmp1206 = Taylor1(constant_term(tmp1204) + constant_term(tmp1205), order)
    tmp1207 = Taylor1(constant_term(M_del_mo[3, 3]) * constant_term(Z_me_del_τ_M), order)
    zmed = Taylor1(constant_term(tmp1206) + constant_term(tmp1207), order)
    tmp1210 = Taylor1(constant_term(xmed) ^ constant_term(2), order)
    tmp1212 = Taylor1(constant_term(ymed) ^ constant_term(2), order)
    tmp1213 = Taylor1(constant_term(tmp1210) + constant_term(tmp1212), order)
    tmp1215 = Taylor1(constant_term(zmed) ^ constant_term(2), order)
    rmed2 = Taylor1(constant_term(tmp1213) + constant_term(tmp1215), order)
    tmp1218 = Taylor1(constant_term(rmed2) ^ constant_term(2.5), order)
    factmed = Taylor1(constant_term(fact_num) / constant_term(tmp1218), order)
    tmp1221 = Taylor1(constant_term(xmed) ^ constant_term(2), order)
    tmp1223 = Taylor1(constant_term(rmed2) / constant_term(3), order)
    tmp1224 = Taylor1(constant_term(tmp1221) - constant_term(tmp1223), order)
    tmp1225 = Taylor1(constant_term(factmed) * constant_term(tmp1224), order)
    ITM_t[1, 1] = Taylor1(constant_term(ITM2_t[1, 1]) + constant_term(tmp1225), order)
    tmp1228 = Taylor1(constant_term(ymed) ^ constant_term(2), order)
    tmp1230 = Taylor1(constant_term(rmed2) / constant_term(3), order)
    tmp1231 = Taylor1(constant_term(tmp1228) - constant_term(tmp1230), order)
    tmp1232 = Taylor1(constant_term(factmed) * constant_term(tmp1231), order)
    ITM_t[2, 2] = Taylor1(constant_term(ITM2_t[2, 2]) + constant_term(tmp1232), order)
    tmp1235 = Taylor1(constant_term(zmed) ^ constant_term(2), order)
    tmp1237 = Taylor1(constant_term(rmed2) / constant_term(3), order)
    tmp1238 = Taylor1(constant_term(tmp1235) - constant_term(tmp1237), order)
    tmp1239 = Taylor1(constant_term(factmed) * constant_term(tmp1238), order)
    ITM_t[3, 3] = Taylor1(constant_term(ITM2_t[3, 3]) + constant_term(tmp1239), order)
    tmp1241 = Taylor1(constant_term(factmed) * constant_term(xmed), order)
    tmp1242 = Taylor1(constant_term(tmp1241) * constant_term(ymed), order)
    ITM_t[1, 2] = Taylor1(constant_term(ITM2_t[1, 2]) + constant_term(tmp1242), order)
    ITM_t[2, 1] = Taylor1(identity(constant_term(ITM_t[1, 2])), order)
    tmp1244 = Taylor1(constant_term(factmed) * constant_term(xmed), order)
    tmp1245 = Taylor1(constant_term(tmp1244) * constant_term(zmed), order)
    ITM_t[1, 3] = Taylor1(constant_term(ITM2_t[1, 3]) + constant_term(tmp1245), order)
    ITM_t[3, 1] = Taylor1(identity(constant_term(ITM_t[1, 3])), order)
    tmp1247 = Taylor1(constant_term(factmed) * constant_term(ymed), order)
    tmp1248 = Taylor1(constant_term(tmp1247) * constant_term(zmed), order)
    ITM_t[2, 3] = Taylor1(constant_term(ITM2_t[2, 3]) + constant_term(tmp1248), order)
    ITM_t[3, 2] = Taylor1(identity(constant_term(ITM_t[2, 3])), order)
    tmp1250 = Taylor1(constant_term(ITM_t[1, 1]) + constant_term(ITM_t[2, 2]), order)
    tmp1252 = Taylor1(constant_term(tmp1250) / constant_term(2), order)
    tmp1253 = Taylor1(constant_term(ITM_t[3, 3]) - constant_term(tmp1252), order)
    J2M_t = Taylor1(constant_term(tmp1253) / constant_term(μ[mo]), order)
    tmp1255 = Taylor1(constant_term(ITM_t[2, 2]) - constant_term(ITM_t[1, 1]), order)
    tmp1256 = Taylor1(constant_term(tmp1255) / constant_term(μ[mo]), order)
    C22M_t = Taylor1(constant_term(tmp1256) / constant_term(4), order)
    tmp1259 = Taylor1(-(constant_term(ITM_t[1, 3])), order)
    C21M_t = Taylor1(constant_term(tmp1259) / constant_term(μ[mo]), order)
    tmp1261 = Taylor1(-(constant_term(ITM_t[3, 2])), order)
    S21M_t = Taylor1(constant_term(tmp1261) / constant_term(μ[mo]), order)
    tmp1263 = Taylor1(-(constant_term(ITM_t[2, 1])), order)
    tmp1264 = Taylor1(constant_term(tmp1263) / constant_term(μ[mo]), order)
    S22M_t = Taylor1(constant_term(tmp1264) / constant_term(2), order)
    J2_t[mo] = Taylor1(identity(constant_term(J2M_t)), order)
    tmp1276 = Array{Taylor1{_S}}(undef, size(X_bf_1))
    tmp1276 .= Taylor1(zero(_S), order)
    tmp1278 = Array{Taylor1{_S}}(undef, size(Y_bf_1))
    tmp1278 .= Taylor1(zero(_S), order)
    tmp1280 = Array{Taylor1{_S}}(undef, size(Z_bf_1))
    tmp1280 .= Taylor1(zero(_S), order)
    tmp1284 = Array{Taylor1{_S}}(undef, size(X_bf))
    tmp1284 .= Taylor1(zero(_S), order)
    tmp1286 = Array{Taylor1{_S}}(undef, size(Y_bf))
    tmp1286 .= Taylor1(zero(_S), order)
    tmp1287 = Array{Taylor1{_S}}(undef, size(tmp1284))
    tmp1287 .= Taylor1(zero(_S), order)
    tmp1292 = Array{Taylor1{_S}}(undef, size(P_n))
    tmp1292 .= Taylor1(zero(_S), order)
    tmp1293 = Array{Taylor1{_S}}(undef, size(tmp1292))
    tmp1293 .= Taylor1(zero(_S), order)
    tmp1294 = Array{Taylor1{_S}}(undef, size(P_n))
    tmp1294 .= Taylor1(zero(_S), order)
    tmp1296 = Array{Taylor1{_S}}(undef, size(dP_n))
    tmp1296 .= Taylor1(zero(_S), order)
    tmp1297 = Array{Taylor1{_S}}(undef, size(P_n))
    tmp1297 .= Taylor1(zero(_S), order)
    tmp1302 = Array{Taylor1{_S}}(undef, size(P_n))
    tmp1302 .= Taylor1(zero(_S), order)
    tmp1303 = Array{Taylor1{_S}}(undef, size(tmp1302))
    tmp1303 .= Taylor1(zero(_S), order)
    tmp1305 = Array{Taylor1{_S}}(undef, size(dP_n))
    tmp1305 .= Taylor1(zero(_S), order)
    tmp1306 = Array{Taylor1{_S}}(undef, size(tmp1305))
    tmp1306 .= Taylor1(zero(_S), order)
    tmp1307 = Array{Taylor1{_S}}(undef, size(tmp1306))
    tmp1307 .= Taylor1(zero(_S), order)
    tmp1309 = Array{Taylor1{_S}}(undef, size(P_n))
    tmp1309 .= Taylor1(zero(_S), order)
    tmp1310 = Array{Taylor1{_S}}(undef, size(tmp1309))
    tmp1310 .= Taylor1(zero(_S), order)
    tmp1311 = Array{Taylor1{_S}}(undef, size(tmp1310))
    tmp1311 .= Taylor1(zero(_S), order)
    tmp1313 = Array{Taylor1{_S}}(undef, size(dP_n))
    tmp1313 .= Taylor1(zero(_S), order)
    tmp1314 = Array{Taylor1{_S}}(undef, size(tmp1313))
    tmp1314 .= Taylor1(zero(_S), order)
    tmp1315 = Array{Taylor1{_S}}(undef, size(tmp1314))
    tmp1315 .= Taylor1(zero(_S), order)
    tmp1316 = Array{Taylor1{_S}}(undef, size(tmp1315))
    tmp1316 .= Taylor1(zero(_S), order)
    tmp1318 = Array{Taylor1{_S}}(undef, size(sin_mλ))
    tmp1318 .= Taylor1(zero(_S), order)
    tmp1319 = Array{Taylor1{_S}}(undef, size(cos_mλ))
    tmp1319 .= Taylor1(zero(_S), order)
    tmp1321 = Array{Taylor1{_S}}(undef, size(cos_mλ))
    tmp1321 .= Taylor1(zero(_S), order)
    tmp1322 = Array{Taylor1{_S}}(undef, size(sin_mλ))
    tmp1322 .= Taylor1(zero(_S), order)
    tmp1324 = Array{Taylor1{_S}}(undef, size(secϕ_P_nm))
    tmp1324 .= Taylor1(zero(_S), order)
    tmp1327 = Array{Taylor1{_S}}(undef, size(secϕ_P_nm))
    tmp1327 .= Taylor1(zero(_S), order)
    tmp1329 = Array{Taylor1{_S}}(undef, size(secϕ_P_nm))
    tmp1329 .= Taylor1(zero(_S), order)
    tmp1331 = Array{Taylor1{_S}}(undef, size(secϕ_P_nm))
    tmp1331 .= Taylor1(zero(_S), order)
    tmp1332 = Array{Taylor1{_S}}(undef, size(tmp1331))
    tmp1332 .= Taylor1(zero(_S), order)
    tmp1333 = Array{Taylor1{_S}}(undef, size(secϕ_P_nm))
    tmp1333 .= Taylor1(zero(_S), order)
    tmp1336 = Array{Taylor1{_S}}(undef, size(secϕ_P_nm))
    tmp1336 .= Taylor1(zero(_S), order)
    tmp1337 = Array{Taylor1{_S}}(undef, size(tmp1336))
    tmp1337 .= Taylor1(zero(_S), order)
    tmp1338 = Array{Taylor1{_S}}(undef, size(secϕ_P_nm))
    tmp1338 .= Taylor1(zero(_S), order)
    tmp1340 = Array{Taylor1{_S}}(undef, size(P_nm))
    tmp1340 .= Taylor1(zero(_S), order)
    tmp1341 = Array{Taylor1{_S}}(undef, size(cos_mλ))
    tmp1341 .= Taylor1(zero(_S), order)
    tmp1342 = Array{Taylor1{_S}}(undef, size(sin_mλ))
    tmp1342 .= Taylor1(zero(_S), order)
    tmp1343 = Array{Taylor1{_S}}(undef, size(tmp1341))
    tmp1343 .= Taylor1(zero(_S), order)
    tmp1344 = Array{Taylor1{_S}}(undef, size(tmp1340))
    tmp1344 .= Taylor1(zero(_S), order)
    tmp1345 = Array{Taylor1{_S}}(undef, size(P_nm))
    tmp1345 .= Taylor1(zero(_S), order)
    tmp1346 = Array{Taylor1{_S}}(undef, size(cos_mλ))
    tmp1346 .= Taylor1(zero(_S), order)
    tmp1347 = Array{Taylor1{_S}}(undef, size(sin_mλ))
    tmp1347 .= Taylor1(zero(_S), order)
    tmp1348 = Array{Taylor1{_S}}(undef, size(tmp1346))
    tmp1348 .= Taylor1(zero(_S), order)
    tmp1349 = Array{Taylor1{_S}}(undef, size(tmp1345))
    tmp1349 .= Taylor1(zero(_S), order)
    tmp1350 = Array{Taylor1{_S}}(undef, size(tmp1344))
    tmp1350 .= Taylor1(zero(_S), order)
    tmp1352 = Array{Taylor1{_S}}(undef, size(secϕ_P_nm))
    tmp1352 .= Taylor1(zero(_S), order)
    tmp1353 = Array{Taylor1{_S}}(undef, size(cos_mλ))
    tmp1353 .= Taylor1(zero(_S), order)
    tmp1354 = Array{Taylor1{_S}}(undef, size(sin_mλ))
    tmp1354 .= Taylor1(zero(_S), order)
    tmp1355 = Array{Taylor1{_S}}(undef, size(tmp1353))
    tmp1355 .= Taylor1(zero(_S), order)
    tmp1356 = Array{Taylor1{_S}}(undef, size(tmp1352))
    tmp1356 .= Taylor1(zero(_S), order)
    tmp1357 = Array{Taylor1{_S}}(undef, size(secϕ_P_nm))
    tmp1357 .= Taylor1(zero(_S), order)
    tmp1358 = Array{Taylor1{_S}}(undef, size(cos_mλ))
    tmp1358 .= Taylor1(zero(_S), order)
    tmp1359 = Array{Taylor1{_S}}(undef, size(sin_mλ))
    tmp1359 .= Taylor1(zero(_S), order)
    tmp1360 = Array{Taylor1{_S}}(undef, size(tmp1358))
    tmp1360 .= Taylor1(zero(_S), order)
    tmp1361 = Array{Taylor1{_S}}(undef, size(tmp1357))
    tmp1361 .= Taylor1(zero(_S), order)
    tmp1362 = Array{Taylor1{_S}}(undef, size(tmp1356))
    tmp1362 .= Taylor1(zero(_S), order)
    tmp1364 = Array{Taylor1{_S}}(undef, size(cos_mλ))
    tmp1364 .= Taylor1(zero(_S), order)
    tmp1365 = Array{Taylor1{_S}}(undef, size(sin_mλ))
    tmp1365 .= Taylor1(zero(_S), order)
    tmp1366 = Array{Taylor1{_S}}(undef, size(tmp1364))
    tmp1366 .= Taylor1(zero(_S), order)
    tmp1367 = Array{Taylor1{_S}}(undef, size(cosϕ_dP_nm))
    tmp1367 .= Taylor1(zero(_S), order)
    tmp1368 = Array{Taylor1{_S}}(undef, size(cos_mλ))
    tmp1368 .= Taylor1(zero(_S), order)
    tmp1369 = Array{Taylor1{_S}}(undef, size(sin_mλ))
    tmp1369 .= Taylor1(zero(_S), order)
    tmp1370 = Array{Taylor1{_S}}(undef, size(tmp1368))
    tmp1370 .= Taylor1(zero(_S), order)
    tmp1371 = Array{Taylor1{_S}}(undef, size(cosϕ_dP_nm))
    tmp1371 .= Taylor1(zero(_S), order)
    tmp1372 = Array{Taylor1{_S}}(undef, size(tmp1367))
    tmp1372 .= Taylor1(zero(_S), order)
    tmp1374 = Array{Taylor1{_S}}(undef, size(P_nm))
    tmp1374 .= Taylor1(zero(_S), order)
    tmp1375 = Array{Taylor1{_S}}(undef, size(cos_mλ))
    tmp1375 .= Taylor1(zero(_S), order)
    tmp1376 = Array{Taylor1{_S}}(undef, size(sin_mλ))
    tmp1376 .= Taylor1(zero(_S), order)
    tmp1377 = Array{Taylor1{_S}}(undef, size(tmp1375))
    tmp1377 .= Taylor1(zero(_S), order)
    tmp1378 = Array{Taylor1{_S}}(undef, size(tmp1374))
    tmp1378 .= Taylor1(zero(_S), order)
    tmp1379 = Array{Taylor1{_S}}(undef, size(tmp1378))
    tmp1379 .= Taylor1(zero(_S), order)
    temp_CS_ξ = Array{Taylor1{_S}}(undef, size(tmp1379))
    temp_CS_ξ .= Taylor1(zero(_S), order)
    tmp1381 = Array{Taylor1{_S}}(undef, size(secϕ_P_nm))
    tmp1381 .= Taylor1(zero(_S), order)
    tmp1382 = Array{Taylor1{_S}}(undef, size(cos_mλ))
    tmp1382 .= Taylor1(zero(_S), order)
    tmp1383 = Array{Taylor1{_S}}(undef, size(sin_mλ))
    tmp1383 .= Taylor1(zero(_S), order)
    tmp1384 = Array{Taylor1{_S}}(undef, size(tmp1382))
    tmp1384 .= Taylor1(zero(_S), order)
    tmp1385 = Array{Taylor1{_S}}(undef, size(tmp1381))
    tmp1385 .= Taylor1(zero(_S), order)
    tmp1386 = Array{Taylor1{_S}}(undef, size(tmp1385))
    tmp1386 .= Taylor1(zero(_S), order)
    temp_CS_η = Array{Taylor1{_S}}(undef, size(tmp1386))
    temp_CS_η .= Taylor1(zero(_S), order)
    tmp1388 = Array{Taylor1{_S}}(undef, size(cos_mλ))
    tmp1388 .= Taylor1(zero(_S), order)
    tmp1389 = Array{Taylor1{_S}}(undef, size(sin_mλ))
    tmp1389 .= Taylor1(zero(_S), order)
    tmp1390 = Array{Taylor1{_S}}(undef, size(tmp1388))
    tmp1390 .= Taylor1(zero(_S), order)
    tmp1391 = Array{Taylor1{_S}}(undef, size(cosϕ_dP_nm))
    tmp1391 .= Taylor1(zero(_S), order)
    tmp1392 = Array{Taylor1{_S}}(undef, size(tmp1391))
    tmp1392 .= Taylor1(zero(_S), order)
    temp_CS_ζ = Array{Taylor1{_S}}(undef, size(tmp1392))
    temp_CS_ζ .= Taylor1(zero(_S), order)
    tmp1394 = Array{Taylor1{_S}}(undef, size(F_J_ξ))
    tmp1394 .= Taylor1(zero(_S), order)
    tmp1395 = Array{Taylor1{_S}}(undef, size(F_CS_ξ))
    tmp1395 .= Taylor1(zero(_S), order)
    tmp1398 = Array{Taylor1{_S}}(undef, size(F_J_ζ))
    tmp1398 .= Taylor1(zero(_S), order)
    tmp1399 = Array{Taylor1{_S}}(undef, size(F_CS_ζ))
    tmp1399 .= Taylor1(zero(_S), order)
    tmp1405 = Array{Taylor1{_S}}(undef, size(sin_ϕ))
    tmp1405 .= Taylor1(zero(_S), order)
    tmp1408 = Array{Taylor1{_S}}(undef, size(sin_ϕ))
    tmp1408 .= Taylor1(zero(_S), order)
    tmp1410 = Array{Taylor1{_S}}(undef, size(Rb2p))
    tmp1410 .= Taylor1(zero(_S), order)
    tmp1411 = Array{Taylor1{_S}}(undef, size(Rb2p))
    tmp1411 .= Taylor1(zero(_S), order)
    tmp1412 = Array{Taylor1{_S}}(undef, size(tmp1410))
    tmp1412 .= Taylor1(zero(_S), order)
    tmp1413 = Array{Taylor1{_S}}(undef, size(Rb2p))
    tmp1413 .= Taylor1(zero(_S), order)
    tmp1415 = Array{Taylor1{_S}}(undef, size(Rb2p))
    tmp1415 .= Taylor1(zero(_S), order)
    tmp1416 = Array{Taylor1{_S}}(undef, size(Rb2p))
    tmp1416 .= Taylor1(zero(_S), order)
    tmp1417 = Array{Taylor1{_S}}(undef, size(tmp1415))
    tmp1417 .= Taylor1(zero(_S), order)
    tmp1418 = Array{Taylor1{_S}}(undef, size(Rb2p))
    tmp1418 .= Taylor1(zero(_S), order)
    tmp1420 = Array{Taylor1{_S}}(undef, size(Rb2p))
    tmp1420 .= Taylor1(zero(_S), order)
    tmp1421 = Array{Taylor1{_S}}(undef, size(Rb2p))
    tmp1421 .= Taylor1(zero(_S), order)
    tmp1422 = Array{Taylor1{_S}}(undef, size(tmp1420))
    tmp1422 .= Taylor1(zero(_S), order)
    tmp1423 = Array{Taylor1{_S}}(undef, size(Rb2p))
    tmp1423 .= Taylor1(zero(_S), order)
    tmp1425 = Array{Taylor1{_S}}(undef, size(Rb2p))
    tmp1425 .= Taylor1(zero(_S), order)
    tmp1426 = Array{Taylor1{_S}}(undef, size(Rb2p))
    tmp1426 .= Taylor1(zero(_S), order)
    tmp1427 = Array{Taylor1{_S}}(undef, size(tmp1425))
    tmp1427 .= Taylor1(zero(_S), order)
    tmp1428 = Array{Taylor1{_S}}(undef, size(Rb2p))
    tmp1428 .= Taylor1(zero(_S), order)
    tmp1430 = Array{Taylor1{_S}}(undef, size(Rb2p))
    tmp1430 .= Taylor1(zero(_S), order)
    tmp1431 = Array{Taylor1{_S}}(undef, size(Rb2p))
    tmp1431 .= Taylor1(zero(_S), order)
    tmp1432 = Array{Taylor1{_S}}(undef, size(tmp1430))
    tmp1432 .= Taylor1(zero(_S), order)
    tmp1433 = Array{Taylor1{_S}}(undef, size(Rb2p))
    tmp1433 .= Taylor1(zero(_S), order)
    tmp1435 = Array{Taylor1{_S}}(undef, size(Rb2p))
    tmp1435 .= Taylor1(zero(_S), order)
    tmp1436 = Array{Taylor1{_S}}(undef, size(Rb2p))
    tmp1436 .= Taylor1(zero(_S), order)
    tmp1437 = Array{Taylor1{_S}}(undef, size(tmp1435))
    tmp1437 .= Taylor1(zero(_S), order)
    tmp1438 = Array{Taylor1{_S}}(undef, size(Rb2p))
    tmp1438 .= Taylor1(zero(_S), order)
    tmp1440 = Array{Taylor1{_S}}(undef, size(Rb2p))
    tmp1440 .= Taylor1(zero(_S), order)
    tmp1441 = Array{Taylor1{_S}}(undef, size(Rb2p))
    tmp1441 .= Taylor1(zero(_S), order)
    tmp1442 = Array{Taylor1{_S}}(undef, size(tmp1440))
    tmp1442 .= Taylor1(zero(_S), order)
    tmp1443 = Array{Taylor1{_S}}(undef, size(Rb2p))
    tmp1443 .= Taylor1(zero(_S), order)
    tmp1445 = Array{Taylor1{_S}}(undef, size(Rb2p))
    tmp1445 .= Taylor1(zero(_S), order)
    tmp1446 = Array{Taylor1{_S}}(undef, size(Rb2p))
    tmp1446 .= Taylor1(zero(_S), order)
    tmp1447 = Array{Taylor1{_S}}(undef, size(tmp1445))
    tmp1447 .= Taylor1(zero(_S), order)
    tmp1448 = Array{Taylor1{_S}}(undef, size(Rb2p))
    tmp1448 .= Taylor1(zero(_S), order)
    tmp1450 = Array{Taylor1{_S}}(undef, size(Rb2p))
    tmp1450 .= Taylor1(zero(_S), order)
    tmp1451 = Array{Taylor1{_S}}(undef, size(Rb2p))
    tmp1451 .= Taylor1(zero(_S), order)
    tmp1452 = Array{Taylor1{_S}}(undef, size(tmp1450))
    tmp1452 .= Taylor1(zero(_S), order)
    tmp1453 = Array{Taylor1{_S}}(undef, size(Rb2p))
    tmp1453 .= Taylor1(zero(_S), order)
    tmp1455 = Array{Taylor1{_S}}(undef, size(Gc2p))
    tmp1455 .= Taylor1(zero(_S), order)
    tmp1456 = Array{Taylor1{_S}}(undef, size(Gc2p))
    tmp1456 .= Taylor1(zero(_S), order)
    tmp1457 = Array{Taylor1{_S}}(undef, size(tmp1455))
    tmp1457 .= Taylor1(zero(_S), order)
    tmp1458 = Array{Taylor1{_S}}(undef, size(Gc2p))
    tmp1458 .= Taylor1(zero(_S), order)
    tmp1460 = Array{Taylor1{_S}}(undef, size(Gc2p))
    tmp1460 .= Taylor1(zero(_S), order)
    tmp1461 = Array{Taylor1{_S}}(undef, size(Gc2p))
    tmp1461 .= Taylor1(zero(_S), order)
    tmp1462 = Array{Taylor1{_S}}(undef, size(tmp1460))
    tmp1462 .= Taylor1(zero(_S), order)
    tmp1463 = Array{Taylor1{_S}}(undef, size(Gc2p))
    tmp1463 .= Taylor1(zero(_S), order)
    tmp1465 = Array{Taylor1{_S}}(undef, size(Gc2p))
    tmp1465 .= Taylor1(zero(_S), order)
    tmp1466 = Array{Taylor1{_S}}(undef, size(Gc2p))
    tmp1466 .= Taylor1(zero(_S), order)
    tmp1467 = Array{Taylor1{_S}}(undef, size(tmp1465))
    tmp1467 .= Taylor1(zero(_S), order)
    tmp1468 = Array{Taylor1{_S}}(undef, size(Gc2p))
    tmp1468 .= Taylor1(zero(_S), order)
    #= REPL[2]:315 =# Threads.@threads for j = 1:N_ext
            for i = 1:N_ext
                if i == j
                    continue
                else
                    if UJ_interaction[i, j]
                        X_bf_1[i, j] = Taylor1(constant_term(X[i, j]) * constant_term(M_[1, 1, j]), order)
                        X_bf_2[i, j] = Taylor1(constant_term(Y[i, j]) * constant_term(M_[1, 2, j]), order)
                        X_bf_3[i, j] = Taylor1(constant_term(Z[i, j]) * constant_term(M_[1, 3, j]), order)
                        Y_bf_1[i, j] = Taylor1(constant_term(X[i, j]) * constant_term(M_[2, 1, j]), order)
                        Y_bf_2[i, j] = Taylor1(constant_term(Y[i, j]) * constant_term(M_[2, 2, j]), order)
                        Y_bf_3[i, j] = Taylor1(constant_term(Z[i, j]) * constant_term(M_[2, 3, j]), order)
                        Z_bf_1[i, j] = Taylor1(constant_term(X[i, j]) * constant_term(M_[3, 1, j]), order)
                        Z_bf_2[i, j] = Taylor1(constant_term(Y[i, j]) * constant_term(M_[3, 2, j]), order)
                        Z_bf_3[i, j] = Taylor1(constant_term(Z[i, j]) * constant_term(M_[3, 3, j]), order)
                        tmp1276[i, j] = Taylor1(constant_term(X_bf_1[i, j]) + constant_term(X_bf_2[i, j]), order)
                        X_bf[i, j] = Taylor1(constant_term(tmp1276[i, j]) + constant_term(X_bf_3[i, j]), order)
                        tmp1278[i, j] = Taylor1(constant_term(Y_bf_1[i, j]) + constant_term(Y_bf_2[i, j]), order)
                        Y_bf[i, j] = Taylor1(constant_term(tmp1278[i, j]) + constant_term(Y_bf_3[i, j]), order)
                        tmp1280[i, j] = Taylor1(constant_term(Z_bf_1[i, j]) + constant_term(Z_bf_2[i, j]), order)
                        Z_bf[i, j] = Taylor1(constant_term(tmp1280[i, j]) + constant_term(Z_bf_3[i, j]), order)
                        sin_ϕ[i, j] = Taylor1(constant_term(Z_bf[i, j]) / constant_term(r_p1d2[i, j]), order)
                        tmp1284[i, j] = Taylor1(constant_term(X_bf[i, j]) ^ constant_term(2), order)
                        tmp1286[i, j] = Taylor1(constant_term(Y_bf[i, j]) ^ constant_term(2), order)
                        tmp1287[i, j] = Taylor1(constant_term(tmp1284[i, j]) + constant_term(tmp1286[i, j]), order)
                        r_xy[i, j] = Taylor1(sqrt(constant_term(tmp1287[i, j])), order)
                        cos_ϕ[i, j] = Taylor1(constant_term(r_xy[i, j]) / constant_term(r_p1d2[i, j]), order)
                        sin_λ[i, j] = Taylor1(constant_term(Y_bf[i, j]) / constant_term(r_xy[i, j]), order)
                        cos_λ[i, j] = Taylor1(constant_term(X_bf[i, j]) / constant_term(r_xy[i, j]), order)
                        P_n[i, j, 1] = Taylor1(identity(constant_term(one_t)), order)
                        P_n[i, j, 2] = Taylor1(identity(constant_term(sin_ϕ[i, j])), order)
                        dP_n[i, j, 1] = Taylor1(identity(constant_term(zero_q_1)), order)
                        dP_n[i, j, 2] = Taylor1(identity(constant_term(one_t)), order)
                        for n = 2:n1SEM[j]
                            tmp1292[i, j, n] = Taylor1(constant_term(P_n[i, j, n]) * constant_term(sin_ϕ[i, j]), order)
                            tmp1293[i, j, n] = Taylor1(constant_term(tmp1292[i, j, n]) * constant_term(fact1_jsem[n]), order)
                            tmp1294[i, j, n - 1] = Taylor1(constant_term(P_n[i, j, n - 1]) * constant_term(fact2_jsem[n]), order)
                            P_n[i, j, n + 1] = Taylor1(constant_term(tmp1293[i, j, n]) - constant_term(tmp1294[i, j, n - 1]), order)
                            tmp1296[i, j, n] = Taylor1(constant_term(dP_n[i, j, n]) * constant_term(sin_ϕ[i, j]), order)
                            tmp1297[i, j, n] = Taylor1(constant_term(P_n[i, j, n]) * constant_term(fact3_jsem[n]), order)
                            dP_n[i, j, n + 1] = Taylor1(constant_term(tmp1296[i, j, n]) + constant_term(tmp1297[i, j, n]), order)
                            temp_rn[i, j, n] = Taylor1(constant_term(r_p1d2[i, j]) ^ constant_term(fact5_jsem[n]), order)
                        end
                        r_p4[i, j] = Taylor1(constant_term(r_p2[i, j]) ^ constant_term(2), order)
                        tmp1302[i, j, 3] = Taylor1(constant_term(P_n[i, j, 3]) * constant_term(fact4_jsem[2]), order)
                        tmp1303[i, j, 3] = Taylor1(constant_term(tmp1302[i, j, 3]) * constant_term(J2_t[j]), order)
                        F_J_ξ[i, j] = Taylor1(constant_term(tmp1303[i, j, 3]) / constant_term(r_p4[i, j]), order)
                        tmp1305[i, j, 3] = Taylor1(-(constant_term(dP_n[i, j, 3])), order)
                        tmp1306[i, j, 3] = Taylor1(constant_term(tmp1305[i, j, 3]) * constant_term(cos_ϕ[i, j]), order)
                        tmp1307[i, j, 3] = Taylor1(constant_term(tmp1306[i, j, 3]) * constant_term(J2_t[j]), order)
                        F_J_ζ[i, j] = Taylor1(constant_term(tmp1307[i, j, 3]) / constant_term(r_p4[i, j]), order)
                        F_J_ξ_36[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                        F_J_ζ_36[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                        for n = 3:n1SEM[j]
                            tmp1309[i, j, n + 1] = Taylor1(constant_term(P_n[i, j, n + 1]) * constant_term(fact4_jsem[n]), order)
                            tmp1310[i, j, n + 1] = Taylor1(constant_term(tmp1309[i, j, n + 1]) * constant_term(JSEM[j, n]), order)
                            tmp1311[i, j, n + 1] = Taylor1(constant_term(tmp1310[i, j, n + 1]) / constant_term(temp_rn[i, j, n]), order)
                            temp_fjξ[i, j, n] = Taylor1(constant_term(tmp1311[i, j, n + 1]) + constant_term(F_J_ξ_36[i, j]), order)
                            tmp1313[i, j, n + 1] = Taylor1(-(constant_term(dP_n[i, j, n + 1])), order)
                            tmp1314[i, j, n + 1] = Taylor1(constant_term(tmp1313[i, j, n + 1]) * constant_term(cos_ϕ[i, j]), order)
                            tmp1315[i, j, n + 1] = Taylor1(constant_term(tmp1314[i, j, n + 1]) * constant_term(JSEM[j, n]), order)
                            tmp1316[i, j, n + 1] = Taylor1(constant_term(tmp1315[i, j, n + 1]) / constant_term(temp_rn[i, j, n]), order)
                            temp_fjζ[i, j, n] = Taylor1(constant_term(tmp1316[i, j, n + 1]) + constant_term(F_J_ζ_36[i, j]), order)
                            F_J_ξ_36[i, j] = Taylor1(identity(constant_term(temp_fjξ[i, j, n])), order)
                            F_J_ζ_36[i, j] = Taylor1(identity(constant_term(temp_fjζ[i, j, n])), order)
                        end
                        if j == mo
                            for m = 1:n1SEM[mo]
                                if m == 1
                                    sin_mλ[i, j, 1] = Taylor1(identity(constant_term(sin_λ[i, j])), order)
                                    cos_mλ[i, j, 1] = Taylor1(identity(constant_term(cos_λ[i, j])), order)
                                    secϕ_P_nm[i, j, 1, 1] = Taylor1(identity(constant_term(one_t)), order)
                                else
                                    tmp1318[i, j, 1] = Taylor1(constant_term(sin_mλ[i, j, 1]) * constant_term(cos_mλ[i, j, m - 1]), order)
                                    tmp1319[i, j, 1] = Taylor1(constant_term(cos_mλ[i, j, 1]) * constant_term(sin_mλ[i, j, m - 1]), order)
                                    sin_mλ[i, j, m] = Taylor1(constant_term(tmp1318[i, j, 1]) + constant_term(tmp1319[i, j, 1]), order)
                                    tmp1321[i, j, 1] = Taylor1(constant_term(cos_mλ[i, j, 1]) * constant_term(cos_mλ[i, j, m - 1]), order)
                                    tmp1322[i, j, 1] = Taylor1(constant_term(sin_mλ[i, j, 1]) * constant_term(sin_mλ[i, j, m - 1]), order)
                                    cos_mλ[i, j, m] = Taylor1(constant_term(tmp1321[i, j, 1]) - constant_term(tmp1322[i, j, 1]), order)
                                    tmp1324[i, j, m - 1, m - 1] = Taylor1(constant_term(secϕ_P_nm[i, j, m - 1, m - 1]) * constant_term(cos_ϕ[i, j]), order)
                                    secϕ_P_nm[i, j, m, m] = Taylor1(constant_term(tmp1324[i, j, m - 1, m - 1]) * constant_term(lnm5[m]), order)
                                    P_nm[i, j, m, m] = Taylor1(constant_term(secϕ_P_nm[i, j, m, m]) * constant_term(cos_ϕ[i, j]), order)
                                    tmp1327[i, j, m, m] = Taylor1(constant_term(secϕ_P_nm[i, j, m, m]) * constant_term(sin_ϕ[i, j]), order)
                                    cosϕ_dP_nm[i, j, m, m] = Taylor1(constant_term(tmp1327[i, j, m, m]) * constant_term(lnm3[m]), order)
                                end
                                for n = m + 1:n1SEM[mo]
                                    if n == m + 1
                                        tmp1329[i, j, n - 1, m] = Taylor1(constant_term(secϕ_P_nm[i, j, n - 1, m]) * constant_term(sin_ϕ[i, j]), order)
                                        secϕ_P_nm[i, j, n, m] = Taylor1(constant_term(tmp1329[i, j, n - 1, m]) * constant_term(lnm1[n, m]), order)
                                    else
                                        tmp1331[i, j, n - 1, m] = Taylor1(constant_term(secϕ_P_nm[i, j, n - 1, m]) * constant_term(sin_ϕ[i, j]), order)
                                        tmp1332[i, j, n - 1, m] = Taylor1(constant_term(tmp1331[i, j, n - 1, m]) * constant_term(lnm1[n, m]), order)
                                        tmp1333[i, j, n - 2, m] = Taylor1(constant_term(secϕ_P_nm[i, j, n - 2, m]) * constant_term(lnm2[n, m]), order)
                                        secϕ_P_nm[i, j, n, m] = Taylor1(constant_term(tmp1332[i, j, n - 1, m]) + constant_term(tmp1333[i, j, n - 2, m]), order)
                                    end
                                    P_nm[i, j, n, m] = Taylor1(constant_term(secϕ_P_nm[i, j, n, m]) * constant_term(cos_ϕ[i, j]), order)
                                    tmp1336[i, j, n, m] = Taylor1(constant_term(secϕ_P_nm[i, j, n, m]) * constant_term(sin_ϕ[i, j]), order)
                                    tmp1337[i, j, n, m] = Taylor1(constant_term(tmp1336[i, j, n, m]) * constant_term(lnm3[n]), order)
                                    tmp1338[i, j, n - 1, m] = Taylor1(constant_term(secϕ_P_nm[i, j, n - 1, m]) * constant_term(lnm4[n, m]), order)
                                    cosϕ_dP_nm[i, j, n, m] = Taylor1(constant_term(tmp1337[i, j, n, m]) + constant_term(tmp1338[i, j, n - 1, m]), order)
                                end
                            end
                            tmp1340[i, j, 2, 1] = Taylor1(constant_term(P_nm[i, j, 2, 1]) * constant_term(lnm6[2]), order)
                            tmp1341[i, j, 1] = Taylor1(constant_term(C21M_t) * constant_term(cos_mλ[i, j, 1]), order)
                            tmp1342[i, j, 1] = Taylor1(constant_term(S21M_t) * constant_term(sin_mλ[i, j, 1]), order)
                            tmp1343[i, j, 1] = Taylor1(constant_term(tmp1341[i, j, 1]) + constant_term(tmp1342[i, j, 1]), order)
                            tmp1344[i, j, 2, 1] = Taylor1(constant_term(tmp1340[i, j, 2, 1]) * constant_term(tmp1343[i, j, 1]), order)
                            tmp1345[i, j, 2, 2] = Taylor1(constant_term(P_nm[i, j, 2, 2]) * constant_term(lnm6[2]), order)
                            tmp1346[i, j, 2] = Taylor1(constant_term(C22M_t) * constant_term(cos_mλ[i, j, 2]), order)
                            tmp1347[i, j, 2] = Taylor1(constant_term(S22M_t) * constant_term(sin_mλ[i, j, 2]), order)
                            tmp1348[i, j, 2] = Taylor1(constant_term(tmp1346[i, j, 2]) + constant_term(tmp1347[i, j, 2]), order)
                            tmp1349[i, j, 2, 2] = Taylor1(constant_term(tmp1345[i, j, 2, 2]) * constant_term(tmp1348[i, j, 2]), order)
                            tmp1350[i, j, 2, 1] = Taylor1(constant_term(tmp1344[i, j, 2, 1]) + constant_term(tmp1349[i, j, 2, 2]), order)
                            F_CS_ξ[i, j] = Taylor1(constant_term(tmp1350[i, j, 2, 1]) / constant_term(r_p4[i, j]), order)
                            tmp1352[i, j, 2, 1] = Taylor1(constant_term(secϕ_P_nm[i, j, 2, 1]) * constant_term(lnm7[1]), order)
                            tmp1353[i, j, 1] = Taylor1(constant_term(S21M_t) * constant_term(cos_mλ[i, j, 1]), order)
                            tmp1354[i, j, 1] = Taylor1(constant_term(C21M_t) * constant_term(sin_mλ[i, j, 1]), order)
                            tmp1355[i, j, 1] = Taylor1(constant_term(tmp1353[i, j, 1]) - constant_term(tmp1354[i, j, 1]), order)
                            tmp1356[i, j, 2, 1] = Taylor1(constant_term(tmp1352[i, j, 2, 1]) * constant_term(tmp1355[i, j, 1]), order)
                            tmp1357[i, j, 2, 2] = Taylor1(constant_term(secϕ_P_nm[i, j, 2, 2]) * constant_term(lnm7[2]), order)
                            tmp1358[i, j, 2] = Taylor1(constant_term(S22M_t) * constant_term(cos_mλ[i, j, 2]), order)
                            tmp1359[i, j, 2] = Taylor1(constant_term(C22M_t) * constant_term(sin_mλ[i, j, 2]), order)
                            tmp1360[i, j, 2] = Taylor1(constant_term(tmp1358[i, j, 2]) - constant_term(tmp1359[i, j, 2]), order)
                            tmp1361[i, j, 2, 2] = Taylor1(constant_term(tmp1357[i, j, 2, 2]) * constant_term(tmp1360[i, j, 2]), order)
                            tmp1362[i, j, 2, 1] = Taylor1(constant_term(tmp1356[i, j, 2, 1]) + constant_term(tmp1361[i, j, 2, 2]), order)
                            F_CS_η[i, j] = Taylor1(constant_term(tmp1362[i, j, 2, 1]) / constant_term(r_p4[i, j]), order)
                            tmp1364[i, j, 1] = Taylor1(constant_term(C21M_t) * constant_term(cos_mλ[i, j, 1]), order)
                            tmp1365[i, j, 1] = Taylor1(constant_term(S21M_t) * constant_term(sin_mλ[i, j, 1]), order)
                            tmp1366[i, j, 1] = Taylor1(constant_term(tmp1364[i, j, 1]) + constant_term(tmp1365[i, j, 1]), order)
                            tmp1367[i, j, 2, 1] = Taylor1(constant_term(cosϕ_dP_nm[i, j, 2, 1]) * constant_term(tmp1366[i, j, 1]), order)
                            tmp1368[i, j, 2] = Taylor1(constant_term(C22M_t) * constant_term(cos_mλ[i, j, 2]), order)
                            tmp1369[i, j, 2] = Taylor1(constant_term(S22M_t) * constant_term(sin_mλ[i, j, 2]), order)
                            tmp1370[i, j, 2] = Taylor1(constant_term(tmp1368[i, j, 2]) + constant_term(tmp1369[i, j, 2]), order)
                            tmp1371[i, j, 2, 2] = Taylor1(constant_term(cosϕ_dP_nm[i, j, 2, 2]) * constant_term(tmp1370[i, j, 2]), order)
                            tmp1372[i, j, 2, 1] = Taylor1(constant_term(tmp1367[i, j, 2, 1]) + constant_term(tmp1371[i, j, 2, 2]), order)
                            F_CS_ζ[i, j] = Taylor1(constant_term(tmp1372[i, j, 2, 1]) / constant_term(r_p4[i, j]), order)
                            F_CS_ξ_36[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                            F_CS_η_36[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                            F_CS_ζ_36[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                            for n = 3:n1SEM[mo]
                                for m = 1:n
                                    tmp1374[i, j, n, m] = Taylor1(constant_term(P_nm[i, j, n, m]) * constant_term(lnm6[n]), order)
                                    tmp1375[i, j, m] = Taylor1(constant_term(cos_mλ[i, j, m]) * constant_term(CM[n, m]), order)
                                    tmp1376[i, j, m] = Taylor1(constant_term(sin_mλ[i, j, m]) * constant_term(SM[n, m]), order)
                                    tmp1377[i, j, m] = Taylor1(constant_term(tmp1375[i, j, m]) + constant_term(tmp1376[i, j, m]), order)
                                    tmp1378[i, j, n, m] = Taylor1(constant_term(tmp1374[i, j, n, m]) * constant_term(tmp1377[i, j, m]), order)
                                    tmp1379[i, j, n, m] = Taylor1(constant_term(tmp1378[i, j, n, m]) / constant_term(temp_rn[i, j, n]), order)
                                    temp_CS_ξ[i, j, n, m] = Taylor1(constant_term(tmp1379[i, j, n, m]) + constant_term(F_CS_ξ_36[i, j]), order)
                                    tmp1381[i, j, n, m] = Taylor1(constant_term(secϕ_P_nm[i, j, n, m]) * constant_term(lnm7[m]), order)
                                    tmp1382[i, j, m] = Taylor1(constant_term(cos_mλ[i, j, m]) * constant_term(SM[n, m]), order)
                                    tmp1383[i, j, m] = Taylor1(constant_term(sin_mλ[i, j, m]) * constant_term(CM[n, m]), order)
                                    tmp1384[i, j, m] = Taylor1(constant_term(tmp1382[i, j, m]) - constant_term(tmp1383[i, j, m]), order)
                                    tmp1385[i, j, n, m] = Taylor1(constant_term(tmp1381[i, j, n, m]) * constant_term(tmp1384[i, j, m]), order)
                                    tmp1386[i, j, n, m] = Taylor1(constant_term(tmp1385[i, j, n, m]) / constant_term(temp_rn[i, j, n]), order)
                                    temp_CS_η[i, j, n, m] = Taylor1(constant_term(tmp1386[i, j, n, m]) + constant_term(F_CS_η_36[i, j]), order)
                                    tmp1388[i, j, m] = Taylor1(constant_term(cos_mλ[i, j, m]) * constant_term(CM[n, m]), order)
                                    tmp1389[i, j, m] = Taylor1(constant_term(sin_mλ[i, j, m]) * constant_term(SM[n, m]), order)
                                    tmp1390[i, j, m] = Taylor1(constant_term(tmp1388[i, j, m]) + constant_term(tmp1389[i, j, m]), order)
                                    tmp1391[i, j, n, m] = Taylor1(constant_term(cosϕ_dP_nm[i, j, n, m]) * constant_term(tmp1390[i, j, m]), order)
                                    tmp1392[i, j, n, m] = Taylor1(constant_term(tmp1391[i, j, n, m]) / constant_term(temp_rn[i, j, n]), order)
                                    temp_CS_ζ[i, j, n, m] = Taylor1(constant_term(tmp1392[i, j, n, m]) + constant_term(F_CS_ζ_36[i, j]), order)
                                    F_CS_ξ_36[i, j] = Taylor1(identity(constant_term(temp_CS_ξ[i, j, n, m])), order)
                                    F_CS_η_36[i, j] = Taylor1(identity(constant_term(temp_CS_η[i, j, n, m])), order)
                                    F_CS_ζ_36[i, j] = Taylor1(identity(constant_term(temp_CS_ζ[i, j, n, m])), order)
                                end
                            end
                            tmp1394[i, j] = Taylor1(constant_term(F_J_ξ[i, j]) + constant_term(F_J_ξ_36[i, j]), order)
                            tmp1395[i, j] = Taylor1(constant_term(F_CS_ξ[i, j]) + constant_term(F_CS_ξ_36[i, j]), order)
                            F_JCS_ξ[i, j] = Taylor1(constant_term(tmp1394[i, j]) + constant_term(tmp1395[i, j]), order)
                            F_JCS_η[i, j] = Taylor1(constant_term(F_CS_η[i, j]) + constant_term(F_CS_η_36[i, j]), order)
                            tmp1398[i, j] = Taylor1(constant_term(F_J_ζ[i, j]) + constant_term(F_J_ζ_36[i, j]), order)
                            tmp1399[i, j] = Taylor1(constant_term(F_CS_ζ[i, j]) + constant_term(F_CS_ζ_36[i, j]), order)
                            F_JCS_ζ[i, j] = Taylor1(constant_term(tmp1398[i, j]) + constant_term(tmp1399[i, j]), order)
                        else
                            F_JCS_ξ[i, j] = Taylor1(constant_term(F_J_ξ[i, j]) + constant_term(F_J_ξ_36[i, j]), order)
                            F_JCS_η[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                            F_JCS_ζ[i, j] = Taylor1(constant_term(F_J_ζ[i, j]) + constant_term(F_J_ζ_36[i, j]), order)
                        end
                        Rb2p[i, j, 1, 1] = Taylor1(constant_term(cos_ϕ[i, j]) * constant_term(cos_λ[i, j]), order)
                        Rb2p[i, j, 2, 1] = Taylor1(-(constant_term(sin_λ[i, j])), order)
                        tmp1405[i, j] = Taylor1(-(constant_term(sin_ϕ[i, j])), order)
                        Rb2p[i, j, 3, 1] = Taylor1(constant_term(tmp1405[i, j]) * constant_term(cos_λ[i, j]), order)
                        Rb2p[i, j, 1, 2] = Taylor1(constant_term(cos_ϕ[i, j]) * constant_term(sin_λ[i, j]), order)
                        Rb2p[i, j, 2, 2] = Taylor1(identity(constant_term(cos_λ[i, j])), order)
                        tmp1408[i, j] = Taylor1(-(constant_term(sin_ϕ[i, j])), order)
                        Rb2p[i, j, 3, 2] = Taylor1(constant_term(tmp1408[i, j]) * constant_term(sin_λ[i, j]), order)
                        Rb2p[i, j, 1, 3] = Taylor1(identity(constant_term(sin_ϕ[i, j])), order)
                        Rb2p[i, j, 2, 3] = Taylor1(identity(constant_term(zero_q_1)), order)
                        Rb2p[i, j, 3, 3] = Taylor1(identity(constant_term(cos_ϕ[i, j])), order)
                        tmp1410[i, j, 1, 1] = Taylor1(constant_term(Rb2p[i, j, 1, 1]) * constant_term(M_[1, 1, j]), order)
                        tmp1411[i, j, 1, 2] = Taylor1(constant_term(Rb2p[i, j, 1, 2]) * constant_term(M_[2, 1, j]), order)
                        tmp1412[i, j, 1, 1] = Taylor1(constant_term(tmp1410[i, j, 1, 1]) + constant_term(tmp1411[i, j, 1, 2]), order)
                        tmp1413[i, j, 1, 3] = Taylor1(constant_term(Rb2p[i, j, 1, 3]) * constant_term(M_[3, 1, j]), order)
                        Gc2p[i, j, 1, 1] = Taylor1(constant_term(tmp1412[i, j, 1, 1]) + constant_term(tmp1413[i, j, 1, 3]), order)
                        tmp1415[i, j, 2, 1] = Taylor1(constant_term(Rb2p[i, j, 2, 1]) * constant_term(M_[1, 1, j]), order)
                        tmp1416[i, j, 2, 2] = Taylor1(constant_term(Rb2p[i, j, 2, 2]) * constant_term(M_[2, 1, j]), order)
                        tmp1417[i, j, 2, 1] = Taylor1(constant_term(tmp1415[i, j, 2, 1]) + constant_term(tmp1416[i, j, 2, 2]), order)
                        tmp1418[i, j, 2, 3] = Taylor1(constant_term(Rb2p[i, j, 2, 3]) * constant_term(M_[3, 1, j]), order)
                        Gc2p[i, j, 2, 1] = Taylor1(constant_term(tmp1417[i, j, 2, 1]) + constant_term(tmp1418[i, j, 2, 3]), order)
                        tmp1420[i, j, 3, 1] = Taylor1(constant_term(Rb2p[i, j, 3, 1]) * constant_term(M_[1, 1, j]), order)
                        tmp1421[i, j, 3, 2] = Taylor1(constant_term(Rb2p[i, j, 3, 2]) * constant_term(M_[2, 1, j]), order)
                        tmp1422[i, j, 3, 1] = Taylor1(constant_term(tmp1420[i, j, 3, 1]) + constant_term(tmp1421[i, j, 3, 2]), order)
                        tmp1423[i, j, 3, 3] = Taylor1(constant_term(Rb2p[i, j, 3, 3]) * constant_term(M_[3, 1, j]), order)
                        Gc2p[i, j, 3, 1] = Taylor1(constant_term(tmp1422[i, j, 3, 1]) + constant_term(tmp1423[i, j, 3, 3]), order)
                        tmp1425[i, j, 1, 1] = Taylor1(constant_term(Rb2p[i, j, 1, 1]) * constant_term(M_[1, 2, j]), order)
                        tmp1426[i, j, 1, 2] = Taylor1(constant_term(Rb2p[i, j, 1, 2]) * constant_term(M_[2, 2, j]), order)
                        tmp1427[i, j, 1, 1] = Taylor1(constant_term(tmp1425[i, j, 1, 1]) + constant_term(tmp1426[i, j, 1, 2]), order)
                        tmp1428[i, j, 1, 3] = Taylor1(constant_term(Rb2p[i, j, 1, 3]) * constant_term(M_[3, 2, j]), order)
                        Gc2p[i, j, 1, 2] = Taylor1(constant_term(tmp1427[i, j, 1, 1]) + constant_term(tmp1428[i, j, 1, 3]), order)
                        tmp1430[i, j, 2, 1] = Taylor1(constant_term(Rb2p[i, j, 2, 1]) * constant_term(M_[1, 2, j]), order)
                        tmp1431[i, j, 2, 2] = Taylor1(constant_term(Rb2p[i, j, 2, 2]) * constant_term(M_[2, 2, j]), order)
                        tmp1432[i, j, 2, 1] = Taylor1(constant_term(tmp1430[i, j, 2, 1]) + constant_term(tmp1431[i, j, 2, 2]), order)
                        tmp1433[i, j, 2, 3] = Taylor1(constant_term(Rb2p[i, j, 2, 3]) * constant_term(M_[3, 2, j]), order)
                        Gc2p[i, j, 2, 2] = Taylor1(constant_term(tmp1432[i, j, 2, 1]) + constant_term(tmp1433[i, j, 2, 3]), order)
                        tmp1435[i, j, 3, 1] = Taylor1(constant_term(Rb2p[i, j, 3, 1]) * constant_term(M_[1, 2, j]), order)
                        tmp1436[i, j, 3, 2] = Taylor1(constant_term(Rb2p[i, j, 3, 2]) * constant_term(M_[2, 2, j]), order)
                        tmp1437[i, j, 3, 1] = Taylor1(constant_term(tmp1435[i, j, 3, 1]) + constant_term(tmp1436[i, j, 3, 2]), order)
                        tmp1438[i, j, 3, 3] = Taylor1(constant_term(Rb2p[i, j, 3, 3]) * constant_term(M_[3, 2, j]), order)
                        Gc2p[i, j, 3, 2] = Taylor1(constant_term(tmp1437[i, j, 3, 1]) + constant_term(tmp1438[i, j, 3, 3]), order)
                        tmp1440[i, j, 1, 1] = Taylor1(constant_term(Rb2p[i, j, 1, 1]) * constant_term(M_[1, 3, j]), order)
                        tmp1441[i, j, 1, 2] = Taylor1(constant_term(Rb2p[i, j, 1, 2]) * constant_term(M_[2, 3, j]), order)
                        tmp1442[i, j, 1, 1] = Taylor1(constant_term(tmp1440[i, j, 1, 1]) + constant_term(tmp1441[i, j, 1, 2]), order)
                        tmp1443[i, j, 1, 3] = Taylor1(constant_term(Rb2p[i, j, 1, 3]) * constant_term(M_[3, 3, j]), order)
                        Gc2p[i, j, 1, 3] = Taylor1(constant_term(tmp1442[i, j, 1, 1]) + constant_term(tmp1443[i, j, 1, 3]), order)
                        tmp1445[i, j, 2, 1] = Taylor1(constant_term(Rb2p[i, j, 2, 1]) * constant_term(M_[1, 3, j]), order)
                        tmp1446[i, j, 2, 2] = Taylor1(constant_term(Rb2p[i, j, 2, 2]) * constant_term(M_[2, 3, j]), order)
                        tmp1447[i, j, 2, 1] = Taylor1(constant_term(tmp1445[i, j, 2, 1]) + constant_term(tmp1446[i, j, 2, 2]), order)
                        tmp1448[i, j, 2, 3] = Taylor1(constant_term(Rb2p[i, j, 2, 3]) * constant_term(M_[3, 3, j]), order)
                        Gc2p[i, j, 2, 3] = Taylor1(constant_term(tmp1447[i, j, 2, 1]) + constant_term(tmp1448[i, j, 2, 3]), order)
                        tmp1450[i, j, 3, 1] = Taylor1(constant_term(Rb2p[i, j, 3, 1]) * constant_term(M_[1, 3, j]), order)
                        tmp1451[i, j, 3, 2] = Taylor1(constant_term(Rb2p[i, j, 3, 2]) * constant_term(M_[2, 3, j]), order)
                        tmp1452[i, j, 3, 1] = Taylor1(constant_term(tmp1450[i, j, 3, 1]) + constant_term(tmp1451[i, j, 3, 2]), order)
                        tmp1453[i, j, 3, 3] = Taylor1(constant_term(Rb2p[i, j, 3, 3]) * constant_term(M_[3, 3, j]), order)
                        Gc2p[i, j, 3, 3] = Taylor1(constant_term(tmp1452[i, j, 3, 1]) + constant_term(tmp1453[i, j, 3, 3]), order)
                        tmp1455[i, j, 1, 1] = Taylor1(constant_term(F_JCS_ξ[i, j]) * constant_term(Gc2p[i, j, 1, 1]), order)
                        tmp1456[i, j, 2, 1] = Taylor1(constant_term(F_JCS_η[i, j]) * constant_term(Gc2p[i, j, 2, 1]), order)
                        tmp1457[i, j, 1, 1] = Taylor1(constant_term(tmp1455[i, j, 1, 1]) + constant_term(tmp1456[i, j, 2, 1]), order)
                        tmp1458[i, j, 3, 1] = Taylor1(constant_term(F_JCS_ζ[i, j]) * constant_term(Gc2p[i, j, 3, 1]), order)
                        F_JCS_x[i, j] = Taylor1(constant_term(tmp1457[i, j, 1, 1]) + constant_term(tmp1458[i, j, 3, 1]), order)
                        tmp1460[i, j, 1, 2] = Taylor1(constant_term(F_JCS_ξ[i, j]) * constant_term(Gc2p[i, j, 1, 2]), order)
                        tmp1461[i, j, 2, 2] = Taylor1(constant_term(F_JCS_η[i, j]) * constant_term(Gc2p[i, j, 2, 2]), order)
                        tmp1462[i, j, 1, 2] = Taylor1(constant_term(tmp1460[i, j, 1, 2]) + constant_term(tmp1461[i, j, 2, 2]), order)
                        tmp1463[i, j, 3, 2] = Taylor1(constant_term(F_JCS_ζ[i, j]) * constant_term(Gc2p[i, j, 3, 2]), order)
                        F_JCS_y[i, j] = Taylor1(constant_term(tmp1462[i, j, 1, 2]) + constant_term(tmp1463[i, j, 3, 2]), order)
                        tmp1465[i, j, 1, 3] = Taylor1(constant_term(F_JCS_ξ[i, j]) * constant_term(Gc2p[i, j, 1, 3]), order)
                        tmp1466[i, j, 2, 3] = Taylor1(constant_term(F_JCS_η[i, j]) * constant_term(Gc2p[i, j, 2, 3]), order)
                        tmp1467[i, j, 1, 3] = Taylor1(constant_term(tmp1465[i, j, 1, 3]) + constant_term(tmp1466[i, j, 2, 3]), order)
                        tmp1468[i, j, 3, 3] = Taylor1(constant_term(F_JCS_ζ[i, j]) * constant_term(Gc2p[i, j, 3, 3]), order)
                        F_JCS_z[i, j] = Taylor1(constant_term(tmp1467[i, j, 1, 3]) + constant_term(tmp1468[i, j, 3, 3]), order)
                    end
                end
            end
        end
    tmp1470 = Array{Taylor1{_S}}(undef, size(F_JCS_x))
    tmp1470 .= Taylor1(zero(_S), order)
    tmp1472 = Array{Taylor1{_S}}(undef, size(F_JCS_y))
    tmp1472 .= Taylor1(zero(_S), order)
    tmp1474 = Array{Taylor1{_S}}(undef, size(F_JCS_z))
    tmp1474 .= Taylor1(zero(_S), order)
    tmp1476 = Array{Taylor1{_S}}(undef, size(F_JCS_x))
    tmp1476 .= Taylor1(zero(_S), order)
    tmp1478 = Array{Taylor1{_S}}(undef, size(F_JCS_y))
    tmp1478 .= Taylor1(zero(_S), order)
    tmp1480 = Array{Taylor1{_S}}(undef, size(F_JCS_z))
    tmp1480 .= Taylor1(zero(_S), order)
    for j = 1:N_ext
        for i = 1:N_ext
            if i == j
                continue
            else
                if UJ_interaction[i, j]
                    tmp1470[i, j] = Taylor1(constant_term(μ[i]) * constant_term(F_JCS_x[i, j]), order)
                    temp_accX_j[i, j] = Taylor1(constant_term(accX[j]) - constant_term(tmp1470[i, j]), order)
                    accX[j] = Taylor1(identity(constant_term(temp_accX_j[i, j])), order)
                    tmp1472[i, j] = Taylor1(constant_term(μ[i]) * constant_term(F_JCS_y[i, j]), order)
                    temp_accY_j[i, j] = Taylor1(constant_term(accY[j]) - constant_term(tmp1472[i, j]), order)
                    accY[j] = Taylor1(identity(constant_term(temp_accY_j[i, j])), order)
                    tmp1474[i, j] = Taylor1(constant_term(μ[i]) * constant_term(F_JCS_z[i, j]), order)
                    temp_accZ_j[i, j] = Taylor1(constant_term(accZ[j]) - constant_term(tmp1474[i, j]), order)
                    accZ[j] = Taylor1(identity(constant_term(temp_accZ_j[i, j])), order)
                    tmp1476[i, j] = Taylor1(constant_term(μ[j]) * constant_term(F_JCS_x[i, j]), order)
                    temp_accX_i[i, j] = Taylor1(constant_term(accX[i]) + constant_term(tmp1476[i, j]), order)
                    accX[i] = Taylor1(identity(constant_term(temp_accX_i[i, j])), order)
                    tmp1478[i, j] = Taylor1(constant_term(μ[j]) * constant_term(F_JCS_y[i, j]), order)
                    temp_accY_i[i, j] = Taylor1(constant_term(accY[i]) + constant_term(tmp1478[i, j]), order)
                    accY[i] = Taylor1(identity(constant_term(temp_accY_i[i, j])), order)
                    tmp1480[i, j] = Taylor1(constant_term(μ[j]) * constant_term(F_JCS_z[i, j]), order)
                    temp_accZ_i[i, j] = Taylor1(constant_term(accZ[i]) + constant_term(tmp1480[i, j]), order)
                    accZ[i] = Taylor1(identity(constant_term(temp_accZ_i[i, j])), order)
                end
            end
        end
    end
    tmp1486 = Array{Taylor1{_S}}(undef, size(v2))
    tmp1486 .= Taylor1(zero(_S), order)
    tmp1487 = Array{Taylor1{_S}}(undef, size(v2))
    tmp1487 .= Taylor1(zero(_S), order)
    tmp1489 = Array{Taylor1{_S}}(undef, size(vi_dot_vj))
    tmp1489 .= Taylor1(zero(_S), order)
    Xij_t_Ui = Array{Taylor1{_S}}(undef, size(X))
    Xij_t_Ui .= Taylor1(zero(_S), order)
    Yij_t_Vi = Array{Taylor1{_S}}(undef, size(Y))
    Yij_t_Vi .= Taylor1(zero(_S), order)
    Zij_t_Wi = Array{Taylor1{_S}}(undef, size(Z))
    Zij_t_Wi .= Taylor1(zero(_S), order)
    tmp1495 = Array{Taylor1{_S}}(undef, size(Xij_t_Ui))
    tmp1495 .= Taylor1(zero(_S), order)
    Rij_dot_Vi = Array{Taylor1{_S}}(undef, size(tmp1495))
    Rij_dot_Vi .= Taylor1(zero(_S), order)
    tmp1498 = Array{Taylor1{_S}}(undef, size(Rij_dot_Vi))
    tmp1498 .= Taylor1(zero(_S), order)
    pn1t7 = Array{Taylor1{_S}}(undef, size(tmp1498))
    pn1t7 .= Taylor1(zero(_S), order)
    tmp1501 = Array{Taylor1{_S}}(undef, size(pn1t7))
    tmp1501 .= Taylor1(zero(_S), order)
    pn1t2_7 = Array{Taylor1{_S}}(undef, size(ϕs_and_vs))
    pn1t2_7 .= Taylor1(zero(_S), order)
    #= REPL[2]:489 =# Threads.@threads for j = 1:N
            for i = 1:N
                if i == j
                    continue
                else
                    _4ϕj[i, j] = Taylor1(constant_term(4) * constant_term(newtonianNb_Potential[j]), order)
                    ϕi_plus_4ϕj[i, j] = Taylor1(constant_term(newtonianNb_Potential[i]) + constant_term(_4ϕj[i, j]), order)
                    tmp1486[i] = Taylor1(constant_term(2) * constant_term(v2[i]), order)
                    tmp1487[j] = Taylor1(constant_term(v2[j]) + constant_term(tmp1486[i]), order)
                    tmp1489[i, j] = Taylor1(constant_term(4) * constant_term(vi_dot_vj[i, j]), order)
                    sj2_plus_2si2_minus_4vivj[i, j] = Taylor1(constant_term(tmp1487[j]) - constant_term(tmp1489[i, j]), order)
                    ϕs_and_vs[i, j] = Taylor1(constant_term(sj2_plus_2si2_minus_4vivj[i, j]) - constant_term(ϕi_plus_4ϕj[i, j]), order)
                    Xij_t_Ui[i, j] = Taylor1(constant_term(X[i, j]) * constant_term(dq[3i - 2]), order)
                    Yij_t_Vi[i, j] = Taylor1(constant_term(Y[i, j]) * constant_term(dq[3i - 1]), order)
                    Zij_t_Wi[i, j] = Taylor1(constant_term(Z[i, j]) * constant_term(dq[3i]), order)
                    tmp1495[i, j] = Taylor1(constant_term(Xij_t_Ui[i, j]) + constant_term(Yij_t_Vi[i, j]), order)
                    Rij_dot_Vi[i, j] = Taylor1(constant_term(tmp1495[i, j]) + constant_term(Zij_t_Wi[i, j]), order)
                    tmp1498[i, j] = Taylor1(constant_term(Rij_dot_Vi[i, j]) ^ constant_term(2), order)
                    pn1t7[i, j] = Taylor1(constant_term(tmp1498[i, j]) / constant_term(r_p2[i, j]), order)
                    tmp1501[i, j] = Taylor1(constant_term(1.5) * constant_term(pn1t7[i, j]), order)
                    pn1t2_7[i, j] = Taylor1(constant_term(ϕs_and_vs[i, j]) - constant_term(tmp1501[i, j]), order)
                    pn1t1_7[i, j] = Taylor1(constant_term(c_p2) + constant_term(pn1t2_7[i, j]), order)
                    for k = 1:postnewton_iter
                        pn1[i, j, k] = Taylor1(identity(constant_term(zero_q_1)), order)
                        X_t_pn1[i, j, k] = Taylor1(identity(constant_term(zero_q_1)), order)
                        Y_t_pn1[i, j, k] = Taylor1(identity(constant_term(zero_q_1)), order)
                        Z_t_pn1[i, j, k] = Taylor1(identity(constant_term(zero_q_1)), order)
                        pNX_t_pn3[i, j, k] = Taylor1(identity(constant_term(zero_q_1)), order)
                        pNY_t_pn3[i, j, k] = Taylor1(identity(constant_term(zero_q_1)), order)
                        pNZ_t_pn3[i, j, k] = Taylor1(identity(constant_term(zero_q_1)), order)
                        pNX_t_X[i, j, k] = Taylor1(identity(constant_term(zero_q_1)), order)
                        pNY_t_Y[i, j, k] = Taylor1(identity(constant_term(zero_q_1)), order)
                        pNZ_t_Z[i, j, k] = Taylor1(identity(constant_term(zero_q_1)), order)
                    end
                end
            end
            postNewtonX[j, 1] = Taylor1(identity(constant_term(newtonX[j])), order)
            postNewtonY[j, 1] = Taylor1(identity(constant_term(newtonY[j])), order)
            postNewtonZ[j, 1] = Taylor1(identity(constant_term(newtonZ[j])), order)
            for k = 1:postnewton_iter
                pntempX[j, k] = Taylor1(identity(constant_term(zero_q_1)), order)
                pntempY[j, k] = Taylor1(identity(constant_term(zero_q_1)), order)
                pntempZ[j, k] = Taylor1(identity(constant_term(zero_q_1)), order)
            end
        end
    tmp1508 = Array{Taylor1{_S}}(undef, size(pNX_t_X))
    tmp1508 .= Taylor1(zero(_S), order)
    tmp1509 = Array{Taylor1{_S}}(undef, size(tmp1508))
    tmp1509 .= Taylor1(zero(_S), order)
    tmp1510 = Array{Taylor1{_S}}(undef, size(tmp1509))
    tmp1510 .= Taylor1(zero(_S), order)
    tmp1518 = Array{Taylor1{_S}}(undef, size(pNX_t_pn3))
    tmp1518 .= Taylor1(zero(_S), order)
    termpnx = Array{Taylor1{_S}}(undef, size(X_t_pn1))
    termpnx .= Taylor1(zero(_S), order)
    sumpnx = Array{Taylor1{_S}}(undef, size(termpnx))
    sumpnx .= Taylor1(zero(_S), order)
    tmp1521 = Array{Taylor1{_S}}(undef, size(pNY_t_pn3))
    tmp1521 .= Taylor1(zero(_S), order)
    termpny = Array{Taylor1{_S}}(undef, size(Y_t_pn1))
    termpny .= Taylor1(zero(_S), order)
    sumpny = Array{Taylor1{_S}}(undef, size(termpny))
    sumpny .= Taylor1(zero(_S), order)
    tmp1524 = Array{Taylor1{_S}}(undef, size(pNZ_t_pn3))
    tmp1524 .= Taylor1(zero(_S), order)
    termpnz = Array{Taylor1{_S}}(undef, size(Z_t_pn1))
    termpnz .= Taylor1(zero(_S), order)
    sumpnz = Array{Taylor1{_S}}(undef, size(termpnz))
    sumpnz .= Taylor1(zero(_S), order)
    for k = 1:postnewton_iter
        #= REPL[2]:534 =# Threads.@threads for j = 1:N
                for i = 1:N
                    if i == j
                        continue
                    else
                        pNX_t_X[i, j, k] = Taylor1(constant_term(postNewtonX[i, k]) * constant_term(X[i, j]), order)
                        pNY_t_Y[i, j, k] = Taylor1(constant_term(postNewtonY[i, k]) * constant_term(Y[i, j]), order)
                        pNZ_t_Z[i, j, k] = Taylor1(constant_term(postNewtonZ[i, k]) * constant_term(Z[i, j]), order)
                        tmp1508[i, j, k] = Taylor1(constant_term(pNX_t_X[i, j, k]) + constant_term(pNY_t_Y[i, j, k]), order)
                        tmp1509[i, j, k] = Taylor1(constant_term(tmp1508[i, j, k]) + constant_term(pNZ_t_Z[i, j, k]), order)
                        tmp1510[i, j, k] = Taylor1(constant_term(0.5) * constant_term(tmp1509[i, j, k]), order)
                        pn1[i, j, k] = Taylor1(constant_term(pn1t1_7[i, j]) + constant_term(tmp1510[i, j, k]), order)
                        X_t_pn1[i, j, k] = Taylor1(constant_term(newton_acc_X[i, j]) * constant_term(pn1[i, j, k]), order)
                        Y_t_pn1[i, j, k] = Taylor1(constant_term(newton_acc_Y[i, j]) * constant_term(pn1[i, j, k]), order)
                        Z_t_pn1[i, j, k] = Taylor1(constant_term(newton_acc_Z[i, j]) * constant_term(pn1[i, j, k]), order)
                        pNX_t_pn3[i, j, k] = Taylor1(constant_term(postNewtonX[i, k]) * constant_term(pn3[i, j]), order)
                        pNY_t_pn3[i, j, k] = Taylor1(constant_term(postNewtonY[i, k]) * constant_term(pn3[i, j]), order)
                        pNZ_t_pn3[i, j, k] = Taylor1(constant_term(postNewtonZ[i, k]) * constant_term(pn3[i, j]), order)
                        tmp1518[i, j, k] = Taylor1(constant_term(U_t_pn2[i, j]) + constant_term(pNX_t_pn3[i, j, k]), order)
                        termpnx[i, j, k] = Taylor1(constant_term(X_t_pn1[i, j, k]) + constant_term(tmp1518[i, j, k]), order)
                        sumpnx[i, j, k] = Taylor1(constant_term(pntempX[j, k]) + constant_term(termpnx[i, j, k]), order)
                        pntempX[j, k] = Taylor1(identity(constant_term(sumpnx[i, j, k])), order)
                        tmp1521[i, j, k] = Taylor1(constant_term(V_t_pn2[i, j]) + constant_term(pNY_t_pn3[i, j, k]), order)
                        termpny[i, j, k] = Taylor1(constant_term(Y_t_pn1[i, j, k]) + constant_term(tmp1521[i, j, k]), order)
                        sumpny[i, j, k] = Taylor1(constant_term(pntempY[j, k]) + constant_term(termpny[i, j, k]), order)
                        pntempY[j, k] = Taylor1(identity(constant_term(sumpny[i, j, k])), order)
                        tmp1524[i, j, k] = Taylor1(constant_term(W_t_pn2[i, j]) + constant_term(pNZ_t_pn3[i, j, k]), order)
                        termpnz[i, j, k] = Taylor1(constant_term(Z_t_pn1[i, j, k]) + constant_term(tmp1524[i, j, k]), order)
                        sumpnz[i, j, k] = Taylor1(constant_term(pntempZ[j, k]) + constant_term(termpnz[i, j, k]), order)
                        pntempZ[j, k] = Taylor1(identity(constant_term(sumpnz[i, j, k])), order)
                    end
                end
                postNewtonX[j, k + 1] = Taylor1(constant_term(pntempX[j, k]) * constant_term(c_m2), order)
                postNewtonY[j, k + 1] = Taylor1(constant_term(pntempY[j, k]) * constant_term(c_m2), order)
                postNewtonZ[j, k + 1] = Taylor1(constant_term(pntempZ[j, k]) * constant_term(c_m2), order)
            end
    end
    X_E_τ_0 = Taylor1(identity(constant_term(q_del_τ_0[3ea - 2])), order)
    Y_E_τ_0 = Taylor1(identity(constant_term(q_del_τ_0[3ea - 1])), order)
    Z_E_τ_0 = Taylor1(identity(constant_term(q_del_τ_0[3ea])), order)
    X_E_τ_1 = Taylor1(identity(constant_term(q_del_τ_1[3ea - 2])), order)
    Y_E_τ_1 = Taylor1(identity(constant_term(q_del_τ_1[3ea - 1])), order)
    Z_E_τ_1 = Taylor1(identity(constant_term(q_del_τ_1[3ea])), order)
    X_E_τ_2 = Taylor1(identity(constant_term(q_del_τ_2[3ea - 2])), order)
    Y_E_τ_2 = Taylor1(identity(constant_term(q_del_τ_2[3ea - 1])), order)
    Z_E_τ_2 = Taylor1(identity(constant_term(q_del_τ_2[3ea])), order)
    X_ME_τ_0 = Taylor1(constant_term(q_del_τ_0[3mo - 2]) - constant_term(X_E_τ_0), order)
    Y_ME_τ_0 = Taylor1(constant_term(q_del_τ_0[3mo - 1]) - constant_term(Y_E_τ_0), order)
    Z_ME_τ_0 = Taylor1(constant_term(q_del_τ_0[3mo]) - constant_term(Z_E_τ_0), order)
    X_ME_τ_1 = Taylor1(constant_term(q_del_τ_1[3mo - 2]) - constant_term(X_E_τ_1), order)
    Y_ME_τ_1 = Taylor1(constant_term(q_del_τ_1[3mo - 1]) - constant_term(Y_E_τ_1), order)
    Z_ME_τ_1 = Taylor1(constant_term(q_del_τ_1[3mo]) - constant_term(Z_E_τ_1), order)
    X_ME_τ_2 = Taylor1(constant_term(q_del_τ_2[3mo - 2]) - constant_term(X_E_τ_2), order)
    Y_ME_τ_2 = Taylor1(constant_term(q_del_τ_2[3mo - 1]) - constant_term(Y_E_τ_2), order)
    Z_ME_τ_2 = Taylor1(constant_term(q_del_τ_2[3mo]) - constant_term(Z_E_τ_2), order)
    X_SE_τ_0 = Taylor1(constant_term(q_del_τ_0[3su - 2]) - constant_term(X_E_τ_0), order)
    Y_SE_τ_0 = Taylor1(constant_term(q_del_τ_0[3su - 1]) - constant_term(Y_E_τ_0), order)
    Z_SE_τ_0 = Taylor1(constant_term(q_del_τ_0[3su]) - constant_term(Z_E_τ_0), order)
    X_SE_τ_1 = Taylor1(constant_term(q_del_τ_1[3su - 2]) - constant_term(X_E_τ_1), order)
    Y_SE_τ_1 = Taylor1(constant_term(q_del_τ_1[3su - 1]) - constant_term(Y_E_τ_1), order)
    Z_SE_τ_1 = Taylor1(constant_term(q_del_τ_1[3su]) - constant_term(Z_E_τ_1), order)
    X_SE_τ_2 = Taylor1(constant_term(q_del_τ_2[3su - 2]) - constant_term(X_E_τ_2), order)
    Y_SE_τ_2 = Taylor1(constant_term(q_del_τ_2[3su - 1]) - constant_term(Y_E_τ_2), order)
    Z_SE_τ_2 = Taylor1(constant_term(q_del_τ_2[3su]) - constant_term(Z_E_τ_2), order)
    tmp1548 = Taylor1(constant_term(R30[1, 1]) * constant_term(X_ME_τ_0), order)
    tmp1549 = Taylor1(constant_term(R30[1, 2]) * constant_term(Y_ME_τ_0), order)
    tmp1550 = Taylor1(constant_term(tmp1548) + constant_term(tmp1549), order)
    tmp1551 = Taylor1(constant_term(R30[1, 3]) * constant_term(Z_ME_τ_0), order)
    r_star_M_0[1] = Taylor1(constant_term(tmp1550) + constant_term(tmp1551), order)
    tmp1553 = Taylor1(constant_term(R30[2, 1]) * constant_term(X_ME_τ_0), order)
    tmp1554 = Taylor1(constant_term(R30[2, 2]) * constant_term(Y_ME_τ_0), order)
    tmp1555 = Taylor1(constant_term(tmp1553) + constant_term(tmp1554), order)
    tmp1556 = Taylor1(constant_term(R30[2, 3]) * constant_term(Z_ME_τ_0), order)
    r_star_M_0[2] = Taylor1(constant_term(tmp1555) + constant_term(tmp1556), order)
    tmp1558 = Taylor1(constant_term(R30[3, 1]) * constant_term(X_ME_τ_0), order)
    tmp1559 = Taylor1(constant_term(R30[3, 2]) * constant_term(Y_ME_τ_0), order)
    tmp1560 = Taylor1(constant_term(tmp1558) + constant_term(tmp1559), order)
    tmp1561 = Taylor1(constant_term(R30[3, 3]) * constant_term(Z_ME_τ_0), order)
    r_star_M_0[3] = Taylor1(constant_term(tmp1560) + constant_term(tmp1561), order)
    tmp1563 = Taylor1(constant_term(R31[1, 1]) * constant_term(X_ME_τ_1), order)
    tmp1564 = Taylor1(constant_term(R31[1, 2]) * constant_term(Y_ME_τ_1), order)
    tmp1565 = Taylor1(constant_term(tmp1563) + constant_term(tmp1564), order)
    tmp1566 = Taylor1(constant_term(R31[1, 3]) * constant_term(Z_ME_τ_1), order)
    r_star_M_1[1] = Taylor1(constant_term(tmp1565) + constant_term(tmp1566), order)
    tmp1568 = Taylor1(constant_term(R31[2, 1]) * constant_term(X_ME_τ_1), order)
    tmp1569 = Taylor1(constant_term(R31[2, 2]) * constant_term(Y_ME_τ_1), order)
    tmp1570 = Taylor1(constant_term(tmp1568) + constant_term(tmp1569), order)
    tmp1571 = Taylor1(constant_term(R31[2, 3]) * constant_term(Z_ME_τ_1), order)
    r_star_M_1[2] = Taylor1(constant_term(tmp1570) + constant_term(tmp1571), order)
    tmp1573 = Taylor1(constant_term(R31[3, 1]) * constant_term(X_ME_τ_1), order)
    tmp1574 = Taylor1(constant_term(R31[3, 2]) * constant_term(Y_ME_τ_1), order)
    tmp1575 = Taylor1(constant_term(tmp1573) + constant_term(tmp1574), order)
    tmp1576 = Taylor1(constant_term(R31[3, 3]) * constant_term(Z_ME_τ_1), order)
    r_star_M_1[3] = Taylor1(constant_term(tmp1575) + constant_term(tmp1576), order)
    tmp1578 = Taylor1(constant_term(R32[1, 1]) * constant_term(X_ME_τ_2), order)
    tmp1579 = Taylor1(constant_term(R32[1, 2]) * constant_term(Y_ME_τ_2), order)
    tmp1580 = Taylor1(constant_term(tmp1578) + constant_term(tmp1579), order)
    tmp1581 = Taylor1(constant_term(R32[1, 3]) * constant_term(Z_ME_τ_2), order)
    r_star_M_2[1] = Taylor1(constant_term(tmp1580) + constant_term(tmp1581), order)
    tmp1583 = Taylor1(constant_term(R32[2, 1]) * constant_term(X_ME_τ_2), order)
    tmp1584 = Taylor1(constant_term(R32[2, 2]) * constant_term(Y_ME_τ_2), order)
    tmp1585 = Taylor1(constant_term(tmp1583) + constant_term(tmp1584), order)
    tmp1586 = Taylor1(constant_term(R32[2, 3]) * constant_term(Z_ME_τ_2), order)
    r_star_M_2[2] = Taylor1(constant_term(tmp1585) + constant_term(tmp1586), order)
    tmp1588 = Taylor1(constant_term(R32[3, 1]) * constant_term(X_ME_τ_2), order)
    tmp1589 = Taylor1(constant_term(R32[3, 2]) * constant_term(Y_ME_τ_2), order)
    tmp1590 = Taylor1(constant_term(tmp1588) + constant_term(tmp1589), order)
    tmp1591 = Taylor1(constant_term(R32[3, 3]) * constant_term(Z_ME_τ_2), order)
    r_star_M_2[3] = Taylor1(constant_term(tmp1590) + constant_term(tmp1591), order)
    tmp1593 = Taylor1(constant_term(R30[1, 1]) * constant_term(X_SE_τ_0), order)
    tmp1594 = Taylor1(constant_term(R30[1, 2]) * constant_term(Y_SE_τ_0), order)
    tmp1595 = Taylor1(constant_term(tmp1593) + constant_term(tmp1594), order)
    tmp1596 = Taylor1(constant_term(R30[1, 3]) * constant_term(Z_SE_τ_0), order)
    r_star_S_0[1] = Taylor1(constant_term(tmp1595) + constant_term(tmp1596), order)
    tmp1598 = Taylor1(constant_term(R30[2, 1]) * constant_term(X_SE_τ_0), order)
    tmp1599 = Taylor1(constant_term(R30[2, 2]) * constant_term(Y_SE_τ_0), order)
    tmp1600 = Taylor1(constant_term(tmp1598) + constant_term(tmp1599), order)
    tmp1601 = Taylor1(constant_term(R30[2, 3]) * constant_term(Z_SE_τ_0), order)
    r_star_S_0[2] = Taylor1(constant_term(tmp1600) + constant_term(tmp1601), order)
    tmp1603 = Taylor1(constant_term(R30[3, 1]) * constant_term(X_SE_τ_0), order)
    tmp1604 = Taylor1(constant_term(R30[3, 2]) * constant_term(Y_SE_τ_0), order)
    tmp1605 = Taylor1(constant_term(tmp1603) + constant_term(tmp1604), order)
    tmp1606 = Taylor1(constant_term(R30[3, 3]) * constant_term(Z_SE_τ_0), order)
    r_star_S_0[3] = Taylor1(constant_term(tmp1605) + constant_term(tmp1606), order)
    tmp1608 = Taylor1(constant_term(R31[1, 1]) * constant_term(X_SE_τ_1), order)
    tmp1609 = Taylor1(constant_term(R31[1, 2]) * constant_term(Y_SE_τ_1), order)
    tmp1610 = Taylor1(constant_term(tmp1608) + constant_term(tmp1609), order)
    tmp1611 = Taylor1(constant_term(R31[1, 3]) * constant_term(Z_SE_τ_1), order)
    r_star_S_1[1] = Taylor1(constant_term(tmp1610) + constant_term(tmp1611), order)
    tmp1613 = Taylor1(constant_term(R31[2, 1]) * constant_term(X_SE_τ_1), order)
    tmp1614 = Taylor1(constant_term(R31[2, 2]) * constant_term(Y_SE_τ_1), order)
    tmp1615 = Taylor1(constant_term(tmp1613) + constant_term(tmp1614), order)
    tmp1616 = Taylor1(constant_term(R31[2, 3]) * constant_term(Z_SE_τ_1), order)
    r_star_S_1[2] = Taylor1(constant_term(tmp1615) + constant_term(tmp1616), order)
    tmp1618 = Taylor1(constant_term(R31[3, 1]) * constant_term(X_SE_τ_1), order)
    tmp1619 = Taylor1(constant_term(R31[3, 2]) * constant_term(Y_SE_τ_1), order)
    tmp1620 = Taylor1(constant_term(tmp1618) + constant_term(tmp1619), order)
    tmp1621 = Taylor1(constant_term(R31[3, 3]) * constant_term(Z_SE_τ_1), order)
    r_star_S_1[3] = Taylor1(constant_term(tmp1620) + constant_term(tmp1621), order)
    tmp1623 = Taylor1(constant_term(R32[1, 1]) * constant_term(X_SE_τ_2), order)
    tmp1624 = Taylor1(constant_term(R32[1, 2]) * constant_term(Y_SE_τ_2), order)
    tmp1625 = Taylor1(constant_term(tmp1623) + constant_term(tmp1624), order)
    tmp1626 = Taylor1(constant_term(R32[1, 3]) * constant_term(Z_SE_τ_2), order)
    r_star_S_2[1] = Taylor1(constant_term(tmp1625) + constant_term(tmp1626), order)
    tmp1628 = Taylor1(constant_term(R32[2, 1]) * constant_term(X_SE_τ_2), order)
    tmp1629 = Taylor1(constant_term(R32[2, 2]) * constant_term(Y_SE_τ_2), order)
    tmp1630 = Taylor1(constant_term(tmp1628) + constant_term(tmp1629), order)
    tmp1631 = Taylor1(constant_term(R32[2, 3]) * constant_term(Z_SE_τ_2), order)
    r_star_S_2[2] = Taylor1(constant_term(tmp1630) + constant_term(tmp1631), order)
    tmp1633 = Taylor1(constant_term(R32[3, 1]) * constant_term(X_SE_τ_2), order)
    tmp1634 = Taylor1(constant_term(R32[3, 2]) * constant_term(Y_SE_τ_2), order)
    tmp1635 = Taylor1(constant_term(tmp1633) + constant_term(tmp1634), order)
    tmp1636 = Taylor1(constant_term(R32[3, 3]) * constant_term(Z_SE_τ_2), order)
    r_star_S_2[3] = Taylor1(constant_term(tmp1635) + constant_term(tmp1636), order)
    tmp1639 = Taylor1(constant_term(r_star_M_0[1]) ^ constant_term(2), order)
    tmp1641 = Taylor1(constant_term(r_star_M_0[2]) ^ constant_term(2), order)
    ρ0s2_M = Taylor1(constant_term(tmp1639) + constant_term(tmp1641), order)
    ρ0s_M = Taylor1(sqrt(constant_term(ρ0s2_M)), order)
    z0s2_M = Taylor1(constant_term(r_star_M_0[3]) ^ constant_term(2), order)
    r0s2_M = Taylor1(constant_term(ρ0s2_M) + constant_term(z0s2_M), order)
    r0s_M = Taylor1(sqrt(constant_term(r0s2_M)), order)
    r0s5_M = Taylor1(constant_term(r0s_M) ^ constant_term(5), order)
    tmp1651 = Taylor1(constant_term(r_star_S_0[1]) ^ constant_term(2), order)
    tmp1653 = Taylor1(constant_term(r_star_S_0[2]) ^ constant_term(2), order)
    ρ0s2_S = Taylor1(constant_term(tmp1651) + constant_term(tmp1653), order)
    ρ0s_S = Taylor1(sqrt(constant_term(ρ0s2_S)), order)
    z0s2_S = Taylor1(constant_term(r_star_S_0[3]) ^ constant_term(2), order)
    r0s2_S = Taylor1(constant_term(ρ0s2_S) + constant_term(z0s2_S), order)
    r0s_S = Taylor1(sqrt(constant_term(r0s2_S)), order)
    r0s5_S = Taylor1(constant_term(r0s_S) ^ constant_term(5), order)
    tmp1663 = Taylor1(constant_term(Z_bf[mo, ea]) * constant_term(r_star_M_0[3]), order)
    tmp1665 = Taylor1(constant_term(tmp1663) ^ constant_term(2), order)
    tmp1667 = Taylor1(constant_term(r_xy[mo, ea]) * constant_term(ρ0s_M), order)
    tmp1669 = Taylor1(constant_term(tmp1667) ^ constant_term(2), order)
    tmp1670 = Taylor1(constant_term(0.5) * constant_term(tmp1669), order)
    tmp1671 = Taylor1(constant_term(tmp1665) + constant_term(tmp1670), order)
    tmp1672 = Taylor1(constant_term(tmp1671) / constant_term(r_p2[mo, ea]), order)
    tmp1673 = Taylor1(constant_term(5) * constant_term(tmp1672), order)
    coeff0_M = Taylor1(constant_term(r0s2_M) - constant_term(tmp1673), order)
    tmp1676 = Taylor1(constant_term(Z_bf[mo, ea]) * constant_term(r_star_S_0[3]), order)
    tmp1678 = Taylor1(constant_term(tmp1676) ^ constant_term(2), order)
    tmp1680 = Taylor1(constant_term(r_xy[mo, ea]) * constant_term(ρ0s_S), order)
    tmp1682 = Taylor1(constant_term(tmp1680) ^ constant_term(2), order)
    tmp1683 = Taylor1(constant_term(0.5) * constant_term(tmp1682), order)
    tmp1684 = Taylor1(constant_term(tmp1678) + constant_term(tmp1683), order)
    tmp1685 = Taylor1(constant_term(tmp1684) / constant_term(r_p2[mo, ea]), order)
    tmp1686 = Taylor1(constant_term(5) * constant_term(tmp1685), order)
    coeff0_S = Taylor1(constant_term(r0s2_S) - constant_term(tmp1686), order)
    k_20E_div_r0s5_M = Taylor1(constant_term(k_20E) / constant_term(r0s5_M), order)
    k_20E_div_r0s5_S = Taylor1(constant_term(k_20E) / constant_term(r0s5_S), order)
    tmp1690 = Taylor1(constant_term(ρ0s2_M) + constant_term(coeff0_M), order)
    tmp1691 = Taylor1(constant_term(k_20E_div_r0s5_M) * constant_term(tmp1690), order)
    aux0_M_x = Taylor1(constant_term(tmp1691) * constant_term(X_bf[mo, ea]), order)
    tmp1693 = Taylor1(constant_term(ρ0s2_M) + constant_term(coeff0_M), order)
    tmp1694 = Taylor1(constant_term(k_20E_div_r0s5_M) * constant_term(tmp1693), order)
    aux0_M_y = Taylor1(constant_term(tmp1694) * constant_term(Y_bf[mo, ea]), order)
    tmp1697 = Taylor1(constant_term(2) * constant_term(z0s2_M), order)
    tmp1698 = Taylor1(constant_term(tmp1697) + constant_term(coeff0_M), order)
    tmp1699 = Taylor1(constant_term(k_20E_div_r0s5_M) * constant_term(tmp1698), order)
    aux0_M_z = Taylor1(constant_term(tmp1699) * constant_term(Z_bf[mo, ea]), order)
    tmp1701 = Taylor1(constant_term(ρ0s2_S) + constant_term(coeff0_S), order)
    tmp1702 = Taylor1(constant_term(k_20E_div_r0s5_S) * constant_term(tmp1701), order)
    aux0_S_x = Taylor1(constant_term(tmp1702) * constant_term(X_bf[mo, ea]), order)
    tmp1704 = Taylor1(constant_term(ρ0s2_S) + constant_term(coeff0_S), order)
    tmp1705 = Taylor1(constant_term(k_20E_div_r0s5_S) * constant_term(tmp1704), order)
    aux0_S_y = Taylor1(constant_term(tmp1705) * constant_term(Y_bf[mo, ea]), order)
    tmp1708 = Taylor1(constant_term(2) * constant_term(z0s2_S), order)
    tmp1709 = Taylor1(constant_term(tmp1708) + constant_term(coeff0_S), order)
    tmp1710 = Taylor1(constant_term(k_20E_div_r0s5_S) * constant_term(tmp1709), order)
    aux0_S_z = Taylor1(constant_term(tmp1710) * constant_term(Z_bf[mo, ea]), order)
    tmp1713 = Taylor1(constant_term(r_star_M_1[1]) ^ constant_term(2), order)
    tmp1715 = Taylor1(constant_term(r_star_M_1[2]) ^ constant_term(2), order)
    ρ1s2_M = Taylor1(constant_term(tmp1713) + constant_term(tmp1715), order)
    ρ1s_M = Taylor1(sqrt(constant_term(ρ1s2_M)), order)
    z1s2_M = Taylor1(constant_term(r_star_M_1[3]) ^ constant_term(2), order)
    r1s2_M = Taylor1(constant_term(ρ1s2_M) + constant_term(z1s2_M), order)
    r1s_M = Taylor1(sqrt(constant_term(r1s2_M)), order)
    r1s5_M = Taylor1(constant_term(r1s_M) ^ constant_term(5), order)
    tmp1725 = Taylor1(constant_term(r_star_S_1[1]) ^ constant_term(2), order)
    tmp1727 = Taylor1(constant_term(r_star_S_1[2]) ^ constant_term(2), order)
    ρ1s2_S = Taylor1(constant_term(tmp1725) + constant_term(tmp1727), order)
    ρ1s_S = Taylor1(sqrt(constant_term(ρ1s2_S)), order)
    z1s2_S = Taylor1(constant_term(r_star_S_1[3]) ^ constant_term(2), order)
    r1s2_S = Taylor1(constant_term(ρ1s2_S) + constant_term(z1s2_S), order)
    r1s_S = Taylor1(sqrt(constant_term(r1s2_S)), order)
    r1s5_S = Taylor1(constant_term(r1s_S) ^ constant_term(5), order)
    tmp1736 = Taylor1(constant_term(X_bf[mo, ea]) * constant_term(r_star_M_1[1]), order)
    tmp1737 = Taylor1(constant_term(Y_bf[mo, ea]) * constant_term(r_star_M_1[2]), order)
    coeff1_1_M = Taylor1(constant_term(tmp1736) + constant_term(tmp1737), order)
    tmp1739 = Taylor1(constant_term(X_bf[mo, ea]) * constant_term(r_star_S_1[1]), order)
    tmp1740 = Taylor1(constant_term(Y_bf[mo, ea]) * constant_term(r_star_S_1[2]), order)
    coeff1_1_S = Taylor1(constant_term(tmp1739) + constant_term(tmp1740), order)
    coeff2_1_M = Taylor1(constant_term(Z_bf[mo, ea]) * constant_term(r_star_M_1[3]), order)
    coeff2_1_S = Taylor1(constant_term(Z_bf[mo, ea]) * constant_term(r_star_S_1[3]), order)
    tmp1745 = Taylor1(constant_term(10) * constant_term(coeff1_1_M), order)
    tmp1746 = Taylor1(constant_term(tmp1745) * constant_term(coeff2_1_M), order)
    coeff3_1_M = Taylor1(constant_term(tmp1746) / constant_term(r_p2[mo, ea]), order)
    tmp1749 = Taylor1(constant_term(10) * constant_term(coeff1_1_S), order)
    tmp1750 = Taylor1(constant_term(tmp1749) * constant_term(coeff2_1_S), order)
    coeff3_1_S = Taylor1(constant_term(tmp1750) / constant_term(r_p2[mo, ea]), order)
    k_21E_div_r1s5_M = Taylor1(constant_term(k_21E) / constant_term(r1s5_M), order)
    k_21E_div_r1s5_S = Taylor1(constant_term(k_21E) / constant_term(r1s5_S), order)
    tmp1755 = Taylor1(constant_term(2) * constant_term(coeff2_1_M), order)
    tmp1756 = Taylor1(constant_term(tmp1755) * constant_term(r_star_M_1[1]), order)
    tmp1757 = Taylor1(constant_term(coeff3_1_M) * constant_term(X_bf[mo, ea]), order)
    tmp1758 = Taylor1(constant_term(tmp1756) - constant_term(tmp1757), order)
    aux1_M_x = Taylor1(constant_term(k_21E_div_r1s5_M) * constant_term(tmp1758), order)
    tmp1761 = Taylor1(constant_term(2) * constant_term(coeff2_1_M), order)
    tmp1762 = Taylor1(constant_term(tmp1761) * constant_term(r_star_M_1[2]), order)
    tmp1763 = Taylor1(constant_term(coeff3_1_M) * constant_term(Y_bf[mo, ea]), order)
    tmp1764 = Taylor1(constant_term(tmp1762) - constant_term(tmp1763), order)
    aux1_M_y = Taylor1(constant_term(k_21E_div_r1s5_M) * constant_term(tmp1764), order)
    tmp1767 = Taylor1(constant_term(2) * constant_term(coeff1_1_M), order)
    tmp1768 = Taylor1(constant_term(tmp1767) * constant_term(r_star_M_1[3]), order)
    tmp1769 = Taylor1(constant_term(coeff3_1_M) * constant_term(Z_bf[mo, ea]), order)
    tmp1770 = Taylor1(constant_term(tmp1768) - constant_term(tmp1769), order)
    aux1_M_z = Taylor1(constant_term(k_21E_div_r1s5_M) * constant_term(tmp1770), order)
    tmp1773 = Taylor1(constant_term(2) * constant_term(coeff2_1_S), order)
    tmp1774 = Taylor1(constant_term(tmp1773) * constant_term(r_star_S_1[1]), order)
    tmp1775 = Taylor1(constant_term(coeff3_1_S) * constant_term(X_bf[mo, ea]), order)
    tmp1776 = Taylor1(constant_term(tmp1774) - constant_term(tmp1775), order)
    aux1_S_x = Taylor1(constant_term(k_21E_div_r1s5_S) * constant_term(tmp1776), order)
    tmp1779 = Taylor1(constant_term(2) * constant_term(coeff2_1_S), order)
    tmp1780 = Taylor1(constant_term(tmp1779) * constant_term(r_star_S_1[2]), order)
    tmp1781 = Taylor1(constant_term(coeff3_1_S) * constant_term(Y_bf[mo, ea]), order)
    tmp1782 = Taylor1(constant_term(tmp1780) - constant_term(tmp1781), order)
    aux1_S_y = Taylor1(constant_term(k_21E_div_r1s5_S) * constant_term(tmp1782), order)
    tmp1785 = Taylor1(constant_term(2) * constant_term(coeff1_1_S), order)
    tmp1786 = Taylor1(constant_term(tmp1785) * constant_term(r_star_S_1[3]), order)
    tmp1787 = Taylor1(constant_term(coeff3_1_S) * constant_term(Z_bf[mo, ea]), order)
    tmp1788 = Taylor1(constant_term(tmp1786) - constant_term(tmp1787), order)
    aux1_S_z = Taylor1(constant_term(k_21E_div_r1s5_S) * constant_term(tmp1788), order)
    tmp1791 = Taylor1(constant_term(r_star_M_2[1]) ^ constant_term(2), order)
    tmp1793 = Taylor1(constant_term(r_star_M_2[2]) ^ constant_term(2), order)
    ρ2s2_M = Taylor1(constant_term(tmp1791) + constant_term(tmp1793), order)
    ρ2s_M = Taylor1(sqrt(constant_term(ρ2s2_M)), order)
    z2s2_M = Taylor1(constant_term(r_star_M_2[3]) ^ constant_term(2), order)
    r2s2_M = Taylor1(constant_term(ρ2s2_M) + constant_term(z2s2_M), order)
    r2s_M = Taylor1(sqrt(constant_term(r2s2_M)), order)
    r2s5_M = Taylor1(constant_term(r2s_M) ^ constant_term(5), order)
    tmp1803 = Taylor1(constant_term(r_star_S_2[1]) ^ constant_term(2), order)
    tmp1805 = Taylor1(constant_term(r_star_S_2[2]) ^ constant_term(2), order)
    ρ2s2_S = Taylor1(constant_term(tmp1803) + constant_term(tmp1805), order)
    ρ2s_S = Taylor1(sqrt(constant_term(ρ2s2_S)), order)
    z2s2_S = Taylor1(constant_term(r_star_S_2[3]) ^ constant_term(2), order)
    r2s2_S = Taylor1(constant_term(ρ2s2_S) + constant_term(z2s2_S), order)
    r2s_S = Taylor1(sqrt(constant_term(r2s2_S)), order)
    r2s5_S = Taylor1(constant_term(r2s_S) ^ constant_term(5), order)
    tmp1814 = Taylor1(constant_term(X_bf[mo, ea]) * constant_term(r_star_M_2[1]), order)
    tmp1815 = Taylor1(constant_term(Y_bf[mo, ea]) * constant_term(r_star_M_2[2]), order)
    coeff1_2_M = Taylor1(constant_term(tmp1814) + constant_term(tmp1815), order)
    tmp1817 = Taylor1(constant_term(X_bf[mo, ea]) * constant_term(r_star_S_2[1]), order)
    tmp1818 = Taylor1(constant_term(Y_bf[mo, ea]) * constant_term(r_star_S_2[2]), order)
    coeff1_2_S = Taylor1(constant_term(tmp1817) + constant_term(tmp1818), order)
    tmp1822 = Taylor1(constant_term(coeff1_2_M) ^ constant_term(2), order)
    tmp1825 = Taylor1(constant_term(r_xy[mo, ea]) ^ constant_term(2), order)
    tmp1826 = Taylor1(constant_term(0.5) * constant_term(tmp1825), order)
    tmp1827 = Taylor1(constant_term(tmp1826) * constant_term(ρ2s2_M), order)
    tmp1828 = Taylor1(constant_term(tmp1822) - constant_term(tmp1827), order)
    tmp1829 = Taylor1(constant_term(5) * constant_term(tmp1828), order)
    coeff3_2_M = Taylor1(constant_term(tmp1829) / constant_term(r_p2[mo, ea]), order)
    tmp1833 = Taylor1(constant_term(coeff1_2_S) ^ constant_term(2), order)
    tmp1836 = Taylor1(constant_term(r_xy[mo, ea]) ^ constant_term(2), order)
    tmp1837 = Taylor1(constant_term(0.5) * constant_term(tmp1836), order)
    tmp1838 = Taylor1(constant_term(tmp1837) * constant_term(ρ2s2_S), order)
    tmp1839 = Taylor1(constant_term(tmp1833) - constant_term(tmp1838), order)
    tmp1840 = Taylor1(constant_term(5) * constant_term(tmp1839), order)
    coeff3_2_S = Taylor1(constant_term(tmp1840) / constant_term(r_p2[mo, ea]), order)
    k_22E_div_r2s5_M = Taylor1(constant_term(k_22E) / constant_term(r2s5_M), order)
    k_22E_div_r2s5_S = Taylor1(constant_term(k_22E) / constant_term(r2s5_S), order)
    tmp1845 = Taylor1(constant_term(2) * constant_term(coeff1_2_M), order)
    tmp1846 = Taylor1(constant_term(tmp1845) * constant_term(r_star_M_2[1]), order)
    tmp1847 = Taylor1(constant_term(ρ2s2_M) + constant_term(coeff3_2_M), order)
    tmp1848 = Taylor1(constant_term(tmp1847) * constant_term(X_bf[mo, ea]), order)
    tmp1849 = Taylor1(constant_term(tmp1846) - constant_term(tmp1848), order)
    aux2_M_x = Taylor1(constant_term(k_22E_div_r2s5_M) * constant_term(tmp1849), order)
    tmp1852 = Taylor1(constant_term(2) * constant_term(coeff1_2_M), order)
    tmp1853 = Taylor1(constant_term(tmp1852) * constant_term(r_star_M_2[2]), order)
    tmp1854 = Taylor1(constant_term(ρ2s2_M) + constant_term(coeff3_2_M), order)
    tmp1855 = Taylor1(constant_term(tmp1854) * constant_term(Y_bf[mo, ea]), order)
    tmp1856 = Taylor1(constant_term(tmp1853) - constant_term(tmp1855), order)
    aux2_M_y = Taylor1(constant_term(k_22E_div_r2s5_M) * constant_term(tmp1856), order)
    tmp1858 = Taylor1(-(constant_term(coeff3_2_M)), order)
    tmp1859 = Taylor1(constant_term(k_22E_div_r2s5_M) * constant_term(tmp1858), order)
    aux2_M_z = Taylor1(constant_term(tmp1859) * constant_term(Z_bf[mo, ea]), order)
    tmp1862 = Taylor1(constant_term(2) * constant_term(coeff1_2_S), order)
    tmp1863 = Taylor1(constant_term(tmp1862) * constant_term(r_star_S_2[1]), order)
    tmp1864 = Taylor1(constant_term(ρ2s2_S) + constant_term(coeff3_2_S), order)
    tmp1865 = Taylor1(constant_term(tmp1864) * constant_term(X_bf[mo, ea]), order)
    tmp1866 = Taylor1(constant_term(tmp1863) - constant_term(tmp1865), order)
    aux2_S_x = Taylor1(constant_term(k_22E_div_r2s5_S) * constant_term(tmp1866), order)
    tmp1869 = Taylor1(constant_term(2) * constant_term(coeff1_2_S), order)
    tmp1870 = Taylor1(constant_term(tmp1869) * constant_term(r_star_S_2[2]), order)
    tmp1871 = Taylor1(constant_term(ρ2s2_S) + constant_term(coeff3_2_S), order)
    tmp1872 = Taylor1(constant_term(tmp1871) * constant_term(Y_bf[mo, ea]), order)
    tmp1873 = Taylor1(constant_term(tmp1870) - constant_term(tmp1872), order)
    aux2_S_y = Taylor1(constant_term(k_22E_div_r2s5_S) * constant_term(tmp1873), order)
    tmp1875 = Taylor1(-(constant_term(coeff3_2_S)), order)
    tmp1876 = Taylor1(constant_term(k_22E_div_r2s5_S) * constant_term(tmp1875), order)
    aux2_S_z = Taylor1(constant_term(tmp1876) * constant_term(Z_bf[mo, ea]), order)
    tmp1878 = Taylor1(constant_term(RE_au) / constant_term(r_p1d2[mo, ea]), order)
    RE_div_r_p5 = Taylor1(constant_term(tmp1878) ^ constant_term(5), order)
    aux_tidacc = Taylor1(constant_term(tid_num_coeff) * constant_term(RE_div_r_p5), order)
    tide_acc_coeff_M = Taylor1(constant_term(μ[mo]) * constant_term(aux_tidacc), order)
    tide_acc_coeff_S = Taylor1(constant_term(μ[su]) * constant_term(aux_tidacc), order)
    tmp1884 = Taylor1(constant_term(aux0_M_x) + constant_term(aux1_M_x), order)
    tmp1885 = Taylor1(constant_term(tmp1884) + constant_term(aux2_M_x), order)
    tmp1886 = Taylor1(constant_term(tide_acc_coeff_M) * constant_term(tmp1885), order)
    tmp1887 = Taylor1(constant_term(aux0_S_x) + constant_term(aux1_S_x), order)
    tmp1888 = Taylor1(constant_term(tmp1887) + constant_term(aux2_S_x), order)
    tmp1889 = Taylor1(constant_term(tide_acc_coeff_S) * constant_term(tmp1888), order)
    tidal_bf_x = Taylor1(constant_term(tmp1886) + constant_term(tmp1889), order)
    tmp1891 = Taylor1(constant_term(aux0_M_y) + constant_term(aux1_M_y), order)
    tmp1892 = Taylor1(constant_term(tmp1891) + constant_term(aux2_M_y), order)
    tmp1893 = Taylor1(constant_term(tide_acc_coeff_M) * constant_term(tmp1892), order)
    tmp1894 = Taylor1(constant_term(aux0_S_y) + constant_term(aux1_S_y), order)
    tmp1895 = Taylor1(constant_term(tmp1894) + constant_term(aux2_S_y), order)
    tmp1896 = Taylor1(constant_term(tide_acc_coeff_S) * constant_term(tmp1895), order)
    tidal_bf_y = Taylor1(constant_term(tmp1893) + constant_term(tmp1896), order)
    tmp1898 = Taylor1(constant_term(aux0_M_z) + constant_term(aux1_M_z), order)
    tmp1899 = Taylor1(constant_term(tmp1898) + constant_term(aux2_M_z), order)
    tmp1900 = Taylor1(constant_term(tide_acc_coeff_M) * constant_term(tmp1899), order)
    tmp1901 = Taylor1(constant_term(aux0_S_z) + constant_term(aux1_S_z), order)
    tmp1902 = Taylor1(constant_term(tmp1901) + constant_term(aux2_S_z), order)
    tmp1903 = Taylor1(constant_term(tide_acc_coeff_S) * constant_term(tmp1902), order)
    tidal_bf_z = Taylor1(constant_term(tmp1900) + constant_term(tmp1903), order)
    tmp1905 = Taylor1(constant_term(M_[1, 1, ea]) * constant_term(tidal_bf_x), order)
    tmp1906 = Taylor1(constant_term(M_[2, 1, ea]) * constant_term(tidal_bf_y), order)
    tmp1907 = Taylor1(constant_term(tmp1905) + constant_term(tmp1906), order)
    tmp1908 = Taylor1(constant_term(M_[3, 1, ea]) * constant_term(tidal_bf_z), order)
    tidal_x = Taylor1(constant_term(tmp1907) + constant_term(tmp1908), order)
    tmp1910 = Taylor1(constant_term(M_[1, 2, ea]) * constant_term(tidal_bf_x), order)
    tmp1911 = Taylor1(constant_term(M_[2, 2, ea]) * constant_term(tidal_bf_y), order)
    tmp1912 = Taylor1(constant_term(tmp1910) + constant_term(tmp1911), order)
    tmp1913 = Taylor1(constant_term(M_[3, 2, ea]) * constant_term(tidal_bf_z), order)
    tidal_y = Taylor1(constant_term(tmp1912) + constant_term(tmp1913), order)
    tmp1915 = Taylor1(constant_term(M_[1, 3, ea]) * constant_term(tidal_bf_x), order)
    tmp1916 = Taylor1(constant_term(M_[2, 3, ea]) * constant_term(tidal_bf_y), order)
    tmp1917 = Taylor1(constant_term(tmp1915) + constant_term(tmp1916), order)
    tmp1918 = Taylor1(constant_term(M_[3, 3, ea]) * constant_term(tidal_bf_z), order)
    tidal_z = Taylor1(constant_term(tmp1917) + constant_term(tmp1918), order)
    accX_mo_tides = Taylor1(constant_term(accX[mo]) + constant_term(tidal_x), order)
    accY_mo_tides = Taylor1(constant_term(accY[mo]) + constant_term(tidal_y), order)
    accZ_mo_tides = Taylor1(constant_term(accZ[mo]) + constant_term(tidal_z), order)
    accX[mo] = Taylor1(identity(constant_term(accX_mo_tides)), order)
    accY[mo] = Taylor1(identity(constant_term(accY_mo_tides)), order)
    accZ[mo] = Taylor1(identity(constant_term(accZ_mo_tides)), order)
    #= REPL[2]:744 =# Threads.@threads for i = 1:N_ext
            dq[3 * (N + i) - 2] = Taylor1(constant_term(postNewtonX[i, postnewton_iter + 1]) + constant_term(accX[i]), order)
            dq[3 * (N + i) - 1] = Taylor1(constant_term(postNewtonY[i, postnewton_iter + 1]) + constant_term(accY[i]), order)
            dq[3 * (N + i)] = Taylor1(constant_term(postNewtonZ[i, postnewton_iter + 1]) + constant_term(accZ[i]), order)
        end
    #= REPL[2]:749 =# Threads.@threads for i = N_ext + 1:N
            dq[3 * (N + i) - 2] = Taylor1(identity(constant_term(postNewtonX[i, postnewton_iter + 1])), order)
            dq[3 * (N + i) - 1] = Taylor1(identity(constant_term(postNewtonY[i, postnewton_iter + 1])), order)
            dq[3 * (N + i)] = Taylor1(identity(constant_term(postNewtonZ[i, postnewton_iter + 1])), order)
        end
    for __idx = eachindex(q)
        (q[__idx]).coeffs[2] = (dq[__idx]).coeffs[1]
    end
    for ord = 1:order - 1
        ordnext = ord + 1
        TaylorSeries.identity!(J2_t[su], J2S_t, ord)
        TaylorSeries.identity!(J2_t[ea], J2E_t, ord)
        #= REPL[2]:211 =# Threads.@threads for j = 1:N
                TaylorSeries.identity!(newtonX[j], zero_q_1, ord)
                TaylorSeries.identity!(newtonY[j], zero_q_1, ord)
                TaylorSeries.identity!(newtonZ[j], zero_q_1, ord)
                TaylorSeries.identity!(newtonianNb_Potential[j], zero_q_1, ord)
                TaylorSeries.identity!(dq[3j - 2], q[3 * (N + j) - 2], ord)
                TaylorSeries.identity!(dq[3j - 1], q[3 * (N + j) - 1], ord)
                TaylorSeries.identity!(dq[3j], q[3 * (N + j)], ord)
            end
        #= REPL[2]:221 =# Threads.@threads for j = 1:N_ext
                TaylorSeries.identity!(accX[j], zero_q_1, ord)
                TaylorSeries.identity!(accY[j], zero_q_1, ord)
                TaylorSeries.identity!(accZ[j], zero_q_1, ord)
            end
        #= REPL[2]:228 =# Threads.@threads for j = 1:N
                for i = 1:N
                    if i == j
                        continue
                    else
                        TaylorSeries.subst!(X[i, j], q[3i - 2], q[3j - 2], ord)
                        TaylorSeries.subst!(Y[i, j], q[3i - 1], q[3j - 1], ord)
                        TaylorSeries.subst!(Z[i, j], q[3i], q[3j], ord)
                        TaylorSeries.subst!(U[i, j], dq[3i - 2], dq[3j - 2], ord)
                        TaylorSeries.subst!(V[i, j], dq[3i - 1], dq[3j - 1], ord)
                        TaylorSeries.subst!(W[i, j], dq[3i], dq[3j], ord)
                        TaylorSeries.mul!(tmp1128[3j - 2], 4, dq[3j - 2], ord)
                        TaylorSeries.mul!(tmp1130[3i - 2], 3, dq[3i - 2], ord)
                        TaylorSeries.subst!(_4U_m_3X[i, j], tmp1128[3j - 2], tmp1130[3i - 2], ord)
                        TaylorSeries.mul!(tmp1133[3j - 1], 4, dq[3j - 1], ord)
                        TaylorSeries.mul!(tmp1135[3i - 1], 3, dq[3i - 1], ord)
                        TaylorSeries.subst!(_4V_m_3Y[i, j], tmp1133[3j - 1], tmp1135[3i - 1], ord)
                        TaylorSeries.mul!(tmp1138[3j], 4, dq[3j], ord)
                        TaylorSeries.mul!(tmp1140[3i], 3, dq[3i], ord)
                        TaylorSeries.subst!(_4W_m_3Z[i, j], tmp1138[3j], tmp1140[3i], ord)
                        TaylorSeries.mul!(pn2x[i, j], X[i, j], _4U_m_3X[i, j], ord)
                        TaylorSeries.mul!(pn2y[i, j], Y[i, j], _4V_m_3Y[i, j], ord)
                        TaylorSeries.mul!(pn2z[i, j], Z[i, j], _4W_m_3Z[i, j], ord)
                        TaylorSeries.mul!(UU[i, j], dq[3i - 2], dq[3j - 2], ord)
                        TaylorSeries.mul!(VV[i, j], dq[3i - 1], dq[3j - 1], ord)
                        TaylorSeries.mul!(WW[i, j], dq[3i], dq[3j], ord)
                        TaylorSeries.add!(tmp1148[i, j], UU[i, j], VV[i, j], ord)
                        TaylorSeries.add!(vi_dot_vj[i, j], tmp1148[i, j], WW[i, j], ord)
                        TaylorSeries.pow!(tmp1151[i, j], X[i, j], 2, ord)
                        TaylorSeries.pow!(tmp1153[i, j], Y[i, j], 2, ord)
                        TaylorSeries.add!(tmp1154[i, j], tmp1151[i, j], tmp1153[i, j], ord)
                        TaylorSeries.pow!(tmp1156[i, j], Z[i, j], 2, ord)
                        TaylorSeries.add!(r_p2[i, j], tmp1154[i, j], tmp1156[i, j], ord)
                        TaylorSeries.sqrt!(r_p1d2[i, j], r_p2[i, j], ord)
                        TaylorSeries.pow!(r_p3d2[i, j], r_p2[i, j], 1.5, ord)
                        TaylorSeries.pow!(r_p7d2[i, j], r_p2[i, j], 3.5, ord)
                        TaylorSeries.div!(newtonianCoeff[i, j], μ[i], r_p3d2[i, j], ord)
                        TaylorSeries.add!(tmp1164[i, j], pn2x[i, j], pn2y[i, j], ord)
                        TaylorSeries.add!(tmp1165[i, j], tmp1164[i, j], pn2z[i, j], ord)
                        TaylorSeries.mul!(pn2[i, j], newtonianCoeff[i, j], tmp1165[i, j], ord)
                        TaylorSeries.mul!(newton_acc_X[i, j], X[i, j], newtonianCoeff[i, j], ord)
                        TaylorSeries.mul!(newton_acc_Y[i, j], Y[i, j], newtonianCoeff[i, j], ord)
                        TaylorSeries.mul!(newton_acc_Z[i, j], Z[i, j], newtonianCoeff[i, j], ord)
                        TaylorSeries.div!(newtonian1b_Potential[i, j], μ[i], r_p1d2[i, j], ord)
                        TaylorSeries.mul!(pn3[i, j], 3.5, newtonian1b_Potential[i, j], ord)
                        TaylorSeries.mul!(U_t_pn2[i, j], pn2[i, j], U[i, j], ord)
                        TaylorSeries.mul!(V_t_pn2[i, j], pn2[i, j], V[i, j], ord)
                        TaylorSeries.mul!(W_t_pn2[i, j], pn2[i, j], W[i, j], ord)
                        TaylorSeries.mul!(tmp1176[i, j], X[i, j], newtonianCoeff[i, j], ord)
                        TaylorSeries.add!(temp_001[i, j], newtonX[j], tmp1176[i, j], ord)
                        TaylorSeries.identity!(newtonX[j], temp_001[i, j], ord)
                        TaylorSeries.mul!(tmp1178[i, j], Y[i, j], newtonianCoeff[i, j], ord)
                        TaylorSeries.add!(temp_002[i, j], newtonY[j], tmp1178[i, j], ord)
                        TaylorSeries.identity!(newtonY[j], temp_002[i, j], ord)
                        TaylorSeries.mul!(tmp1180[i, j], Z[i, j], newtonianCoeff[i, j], ord)
                        TaylorSeries.add!(temp_003[i, j], newtonZ[j], tmp1180[i, j], ord)
                        TaylorSeries.identity!(newtonZ[j], temp_003[i, j], ord)
                        TaylorSeries.add!(temp_004[i, j], newtonianNb_Potential[j], newtonian1b_Potential[i, j], ord)
                        TaylorSeries.identity!(newtonianNb_Potential[j], temp_004[i, j], ord)
                    end
                end
                TaylorSeries.pow!(tmp1184[3j - 2], dq[3j - 2], 2, ord)
                TaylorSeries.pow!(tmp1186[3j - 1], dq[3j - 1], 2, ord)
                TaylorSeries.add!(tmp1187[3j - 2], tmp1184[3j - 2], tmp1186[3j - 1], ord)
                TaylorSeries.pow!(tmp1189[3j], dq[3j], 2, ord)
                TaylorSeries.add!(v2[j], tmp1187[3j - 2], tmp1189[3j], ord)
            end
        TaylorSeries.subst!(X_me_del_τ_M, q_del_τ_M[3mo - 2], q_del_τ_M[3ea - 2], ord)
        TaylorSeries.subst!(Y_me_del_τ_M, q_del_τ_M[3mo - 1], q_del_τ_M[3ea - 1], ord)
        TaylorSeries.subst!(Z_me_del_τ_M, q_del_τ_M[3mo], q_del_τ_M[3ea], ord)
        TaylorSeries.mul!(tmp1194, M_del_mo[1, 1], X_me_del_τ_M, ord)
        TaylorSeries.mul!(tmp1195, M_del_mo[1, 2], Y_me_del_τ_M, ord)
        TaylorSeries.add!(tmp1196, tmp1194, tmp1195, ord)
        TaylorSeries.mul!(tmp1197, M_del_mo[1, 3], Z_me_del_τ_M, ord)
        TaylorSeries.add!(xmed, tmp1196, tmp1197, ord)
        TaylorSeries.mul!(tmp1199, M_del_mo[2, 1], X_me_del_τ_M, ord)
        TaylorSeries.mul!(tmp1200, M_del_mo[2, 2], Y_me_del_τ_M, ord)
        TaylorSeries.add!(tmp1201, tmp1199, tmp1200, ord)
        TaylorSeries.mul!(tmp1202, M_del_mo[2, 3], Z_me_del_τ_M, ord)
        TaylorSeries.add!(ymed, tmp1201, tmp1202, ord)
        TaylorSeries.mul!(tmp1204, M_del_mo[3, 1], X_me_del_τ_M, ord)
        TaylorSeries.mul!(tmp1205, M_del_mo[3, 2], Y_me_del_τ_M, ord)
        TaylorSeries.add!(tmp1206, tmp1204, tmp1205, ord)
        TaylorSeries.mul!(tmp1207, M_del_mo[3, 3], Z_me_del_τ_M, ord)
        TaylorSeries.add!(zmed, tmp1206, tmp1207, ord)
        TaylorSeries.pow!(tmp1210, xmed, 2, ord)
        TaylorSeries.pow!(tmp1212, ymed, 2, ord)
        TaylorSeries.add!(tmp1213, tmp1210, tmp1212, ord)
        TaylorSeries.pow!(tmp1215, zmed, 2, ord)
        TaylorSeries.add!(rmed2, tmp1213, tmp1215, ord)
        TaylorSeries.pow!(tmp1218, rmed2, 2.5, ord)
        TaylorSeries.div!(factmed, fact_num, tmp1218, ord)
        TaylorSeries.pow!(tmp1221, xmed, 2, ord)
        TaylorSeries.div!(tmp1223, rmed2, 3, ord)
        TaylorSeries.subst!(tmp1224, tmp1221, tmp1223, ord)
        TaylorSeries.mul!(tmp1225, factmed, tmp1224, ord)
        TaylorSeries.add!(ITM_t[1, 1], ITM2_t[1, 1], tmp1225, ord)
        TaylorSeries.pow!(tmp1228, ymed, 2, ord)
        TaylorSeries.div!(tmp1230, rmed2, 3, ord)
        TaylorSeries.subst!(tmp1231, tmp1228, tmp1230, ord)
        TaylorSeries.mul!(tmp1232, factmed, tmp1231, ord)
        TaylorSeries.add!(ITM_t[2, 2], ITM2_t[2, 2], tmp1232, ord)
        TaylorSeries.pow!(tmp1235, zmed, 2, ord)
        TaylorSeries.div!(tmp1237, rmed2, 3, ord)
        TaylorSeries.subst!(tmp1238, tmp1235, tmp1237, ord)
        TaylorSeries.mul!(tmp1239, factmed, tmp1238, ord)
        TaylorSeries.add!(ITM_t[3, 3], ITM2_t[3, 3], tmp1239, ord)
        TaylorSeries.mul!(tmp1241, factmed, xmed, ord)
        TaylorSeries.mul!(tmp1242, tmp1241, ymed, ord)
        TaylorSeries.add!(ITM_t[1, 2], ITM2_t[1, 2], tmp1242, ord)
        TaylorSeries.identity!(ITM_t[2, 1], ITM_t[1, 2], ord)
        TaylorSeries.mul!(tmp1244, factmed, xmed, ord)
        TaylorSeries.mul!(tmp1245, tmp1244, zmed, ord)
        TaylorSeries.add!(ITM_t[1, 3], ITM2_t[1, 3], tmp1245, ord)
        TaylorSeries.identity!(ITM_t[3, 1], ITM_t[1, 3], ord)
        TaylorSeries.mul!(tmp1247, factmed, ymed, ord)
        TaylorSeries.mul!(tmp1248, tmp1247, zmed, ord)
        TaylorSeries.add!(ITM_t[2, 3], ITM2_t[2, 3], tmp1248, ord)
        TaylorSeries.identity!(ITM_t[3, 2], ITM_t[2, 3], ord)
        TaylorSeries.add!(tmp1250, ITM_t[1, 1], ITM_t[2, 2], ord)
        TaylorSeries.div!(tmp1252, tmp1250, 2, ord)
        TaylorSeries.subst!(tmp1253, ITM_t[3, 3], tmp1252, ord)
        TaylorSeries.div!(J2M_t, tmp1253, μ[mo], ord)
        TaylorSeries.subst!(tmp1255, ITM_t[2, 2], ITM_t[1, 1], ord)
        TaylorSeries.div!(tmp1256, tmp1255, μ[mo], ord)
        TaylorSeries.div!(C22M_t, tmp1256, 4, ord)
        TaylorSeries.subst!(tmp1259, ITM_t[1, 3], ord)
        TaylorSeries.div!(C21M_t, tmp1259, μ[mo], ord)
        TaylorSeries.subst!(tmp1261, ITM_t[3, 2], ord)
        TaylorSeries.div!(S21M_t, tmp1261, μ[mo], ord)
        TaylorSeries.subst!(tmp1263, ITM_t[2, 1], ord)
        TaylorSeries.div!(tmp1264, tmp1263, μ[mo], ord)
        TaylorSeries.div!(S22M_t, tmp1264, 2, ord)
        TaylorSeries.identity!(J2_t[mo], J2M_t, ord)
        #= REPL[2]:315 =# Threads.@threads for j = 1:N_ext
                for i = 1:N_ext
                    if i == j
                        continue
                    else
                        if UJ_interaction[i, j]
                            TaylorSeries.mul!(X_bf_1[i, j], X[i, j], M_[1, 1, j], ord)
                            TaylorSeries.mul!(X_bf_2[i, j], Y[i, j], M_[1, 2, j], ord)
                            TaylorSeries.mul!(X_bf_3[i, j], Z[i, j], M_[1, 3, j], ord)
                            TaylorSeries.mul!(Y_bf_1[i, j], X[i, j], M_[2, 1, j], ord)
                            TaylorSeries.mul!(Y_bf_2[i, j], Y[i, j], M_[2, 2, j], ord)
                            TaylorSeries.mul!(Y_bf_3[i, j], Z[i, j], M_[2, 3, j], ord)
                            TaylorSeries.mul!(Z_bf_1[i, j], X[i, j], M_[3, 1, j], ord)
                            TaylorSeries.mul!(Z_bf_2[i, j], Y[i, j], M_[3, 2, j], ord)
                            TaylorSeries.mul!(Z_bf_3[i, j], Z[i, j], M_[3, 3, j], ord)
                            TaylorSeries.add!(tmp1276[i, j], X_bf_1[i, j], X_bf_2[i, j], ord)
                            TaylorSeries.add!(X_bf[i, j], tmp1276[i, j], X_bf_3[i, j], ord)
                            TaylorSeries.add!(tmp1278[i, j], Y_bf_1[i, j], Y_bf_2[i, j], ord)
                            TaylorSeries.add!(Y_bf[i, j], tmp1278[i, j], Y_bf_3[i, j], ord)
                            TaylorSeries.add!(tmp1280[i, j], Z_bf_1[i, j], Z_bf_2[i, j], ord)
                            TaylorSeries.add!(Z_bf[i, j], tmp1280[i, j], Z_bf_3[i, j], ord)
                            TaylorSeries.div!(sin_ϕ[i, j], Z_bf[i, j], r_p1d2[i, j], ord)
                            TaylorSeries.pow!(tmp1284[i, j], X_bf[i, j], 2, ord)
                            TaylorSeries.pow!(tmp1286[i, j], Y_bf[i, j], 2, ord)
                            TaylorSeries.add!(tmp1287[i, j], tmp1284[i, j], tmp1286[i, j], ord)
                            TaylorSeries.sqrt!(r_xy[i, j], tmp1287[i, j], ord)
                            TaylorSeries.div!(cos_ϕ[i, j], r_xy[i, j], r_p1d2[i, j], ord)
                            TaylorSeries.div!(sin_λ[i, j], Y_bf[i, j], r_xy[i, j], ord)
                            TaylorSeries.div!(cos_λ[i, j], X_bf[i, j], r_xy[i, j], ord)
                            TaylorSeries.identity!(P_n[i, j, 1], one_t, ord)
                            TaylorSeries.identity!(P_n[i, j, 2], sin_ϕ[i, j], ord)
                            TaylorSeries.identity!(dP_n[i, j, 1], zero_q_1, ord)
                            TaylorSeries.identity!(dP_n[i, j, 2], one_t, ord)
                            for n = 2:n1SEM[j]
                                TaylorSeries.mul!(tmp1292[i, j, n], P_n[i, j, n], sin_ϕ[i, j], ord)
                                TaylorSeries.mul!(tmp1293[i, j, n], tmp1292[i, j, n], fact1_jsem[n], ord)
                                TaylorSeries.mul!(tmp1294[i, j, n - 1], P_n[i, j, n - 1], fact2_jsem[n], ord)
                                TaylorSeries.subst!(P_n[i, j, n + 1], tmp1293[i, j, n], tmp1294[i, j, n - 1], ord)
                                TaylorSeries.mul!(tmp1296[i, j, n], dP_n[i, j, n], sin_ϕ[i, j], ord)
                                TaylorSeries.mul!(tmp1297[i, j, n], P_n[i, j, n], fact3_jsem[n], ord)
                                TaylorSeries.add!(dP_n[i, j, n + 1], tmp1296[i, j, n], tmp1297[i, j, n], ord)
                                TaylorSeries.pow!(temp_rn[i, j, n], r_p1d2[i, j], fact5_jsem[n], ord)
                            end
                            TaylorSeries.pow!(r_p4[i, j], r_p2[i, j], 2, ord)
                            TaylorSeries.mul!(tmp1302[i, j, 3], P_n[i, j, 3], fact4_jsem[2], ord)
                            TaylorSeries.mul!(tmp1303[i, j, 3], tmp1302[i, j, 3], J2_t[j], ord)
                            TaylorSeries.div!(F_J_ξ[i, j], tmp1303[i, j, 3], r_p4[i, j], ord)
                            TaylorSeries.subst!(tmp1305[i, j, 3], dP_n[i, j, 3], ord)
                            TaylorSeries.mul!(tmp1306[i, j, 3], tmp1305[i, j, 3], cos_ϕ[i, j], ord)
                            TaylorSeries.mul!(tmp1307[i, j, 3], tmp1306[i, j, 3], J2_t[j], ord)
                            TaylorSeries.div!(F_J_ζ[i, j], tmp1307[i, j, 3], r_p4[i, j], ord)
                            TaylorSeries.identity!(F_J_ξ_36[i, j], zero_q_1, ord)
                            TaylorSeries.identity!(F_J_ζ_36[i, j], zero_q_1, ord)
                            for n = 3:n1SEM[j]
                                TaylorSeries.mul!(tmp1309[i, j, n + 1], P_n[i, j, n + 1], fact4_jsem[n], ord)
                                TaylorSeries.mul!(tmp1310[i, j, n + 1], tmp1309[i, j, n + 1], JSEM[j, n], ord)
                                TaylorSeries.div!(tmp1311[i, j, n + 1], tmp1310[i, j, n + 1], temp_rn[i, j, n], ord)
                                TaylorSeries.add!(temp_fjξ[i, j, n], tmp1311[i, j, n + 1], F_J_ξ_36[i, j], ord)
                                TaylorSeries.subst!(tmp1313[i, j, n + 1], dP_n[i, j, n + 1], ord)
                                TaylorSeries.mul!(tmp1314[i, j, n + 1], tmp1313[i, j, n + 1], cos_ϕ[i, j], ord)
                                TaylorSeries.mul!(tmp1315[i, j, n + 1], tmp1314[i, j, n + 1], JSEM[j, n], ord)
                                TaylorSeries.div!(tmp1316[i, j, n + 1], tmp1315[i, j, n + 1], temp_rn[i, j, n], ord)
                                TaylorSeries.add!(temp_fjζ[i, j, n], tmp1316[i, j, n + 1], F_J_ζ_36[i, j], ord)
                                TaylorSeries.identity!(F_J_ξ_36[i, j], temp_fjξ[i, j, n], ord)
                                TaylorSeries.identity!(F_J_ζ_36[i, j], temp_fjζ[i, j, n], ord)
                            end
                            if j == mo
                                for m = 1:n1SEM[mo]
                                    if m == 1
                                        TaylorSeries.identity!(sin_mλ[i, j, 1], sin_λ[i, j], ord)
                                        TaylorSeries.identity!(cos_mλ[i, j, 1], cos_λ[i, j], ord)
                                        TaylorSeries.identity!(secϕ_P_nm[i, j, 1, 1], one_t, ord)
                                    else
                                        TaylorSeries.mul!(tmp1318[i, j, 1], sin_mλ[i, j, 1], cos_mλ[i, j, m - 1], ord)
                                        TaylorSeries.mul!(tmp1319[i, j, 1], cos_mλ[i, j, 1], sin_mλ[i, j, m - 1], ord)
                                        TaylorSeries.add!(sin_mλ[i, j, m], tmp1318[i, j, 1], tmp1319[i, j, 1], ord)
                                        TaylorSeries.mul!(tmp1321[i, j, 1], cos_mλ[i, j, 1], cos_mλ[i, j, m - 1], ord)
                                        TaylorSeries.mul!(tmp1322[i, j, 1], sin_mλ[i, j, 1], sin_mλ[i, j, m - 1], ord)
                                        TaylorSeries.subst!(cos_mλ[i, j, m], tmp1321[i, j, 1], tmp1322[i, j, 1], ord)
                                        TaylorSeries.mul!(tmp1324[i, j, m - 1, m - 1], secϕ_P_nm[i, j, m - 1, m - 1], cos_ϕ[i, j], ord)
                                        TaylorSeries.mul!(secϕ_P_nm[i, j, m, m], tmp1324[i, j, m - 1, m - 1], lnm5[m], ord)
                                        TaylorSeries.mul!(P_nm[i, j, m, m], secϕ_P_nm[i, j, m, m], cos_ϕ[i, j], ord)
                                        TaylorSeries.mul!(tmp1327[i, j, m, m], secϕ_P_nm[i, j, m, m], sin_ϕ[i, j], ord)
                                        TaylorSeries.mul!(cosϕ_dP_nm[i, j, m, m], tmp1327[i, j, m, m], lnm3[m], ord)
                                    end
                                    for n = m + 1:n1SEM[mo]
                                        if n == m + 1
                                            TaylorSeries.mul!(tmp1329[i, j, n - 1, m], secϕ_P_nm[i, j, n - 1, m], sin_ϕ[i, j], ord)
                                            TaylorSeries.mul!(secϕ_P_nm[i, j, n, m], tmp1329[i, j, n - 1, m], lnm1[n, m], ord)
                                        else
                                            TaylorSeries.mul!(tmp1331[i, j, n - 1, m], secϕ_P_nm[i, j, n - 1, m], sin_ϕ[i, j], ord)
                                            TaylorSeries.mul!(tmp1332[i, j, n - 1, m], tmp1331[i, j, n - 1, m], lnm1[n, m], ord)
                                            TaylorSeries.mul!(tmp1333[i, j, n - 2, m], secϕ_P_nm[i, j, n - 2, m], lnm2[n, m], ord)
                                            TaylorSeries.add!(secϕ_P_nm[i, j, n, m], tmp1332[i, j, n - 1, m], tmp1333[i, j, n - 2, m], ord)
                                        end
                                        TaylorSeries.mul!(P_nm[i, j, n, m], secϕ_P_nm[i, j, n, m], cos_ϕ[i, j], ord)
                                        TaylorSeries.mul!(tmp1336[i, j, n, m], secϕ_P_nm[i, j, n, m], sin_ϕ[i, j], ord)
                                        TaylorSeries.mul!(tmp1337[i, j, n, m], tmp1336[i, j, n, m], lnm3[n], ord)
                                        TaylorSeries.mul!(tmp1338[i, j, n - 1, m], secϕ_P_nm[i, j, n - 1, m], lnm4[n, m], ord)
                                        TaylorSeries.add!(cosϕ_dP_nm[i, j, n, m], tmp1337[i, j, n, m], tmp1338[i, j, n - 1, m], ord)
                                    end
                                end
                                TaylorSeries.mul!(tmp1340[i, j, 2, 1], P_nm[i, j, 2, 1], lnm6[2], ord)
                                TaylorSeries.mul!(tmp1341[i, j, 1], C21M_t, cos_mλ[i, j, 1], ord)
                                TaylorSeries.mul!(tmp1342[i, j, 1], S21M_t, sin_mλ[i, j, 1], ord)
                                TaylorSeries.add!(tmp1343[i, j, 1], tmp1341[i, j, 1], tmp1342[i, j, 1], ord)
                                TaylorSeries.mul!(tmp1344[i, j, 2, 1], tmp1340[i, j, 2, 1], tmp1343[i, j, 1], ord)
                                TaylorSeries.mul!(tmp1345[i, j, 2, 2], P_nm[i, j, 2, 2], lnm6[2], ord)
                                TaylorSeries.mul!(tmp1346[i, j, 2], C22M_t, cos_mλ[i, j, 2], ord)
                                TaylorSeries.mul!(tmp1347[i, j, 2], S22M_t, sin_mλ[i, j, 2], ord)
                                TaylorSeries.add!(tmp1348[i, j, 2], tmp1346[i, j, 2], tmp1347[i, j, 2], ord)
                                TaylorSeries.mul!(tmp1349[i, j, 2, 2], tmp1345[i, j, 2, 2], tmp1348[i, j, 2], ord)
                                TaylorSeries.add!(tmp1350[i, j, 2, 1], tmp1344[i, j, 2, 1], tmp1349[i, j, 2, 2], ord)
                                TaylorSeries.div!(F_CS_ξ[i, j], tmp1350[i, j, 2, 1], r_p4[i, j], ord)
                                TaylorSeries.mul!(tmp1352[i, j, 2, 1], secϕ_P_nm[i, j, 2, 1], lnm7[1], ord)
                                TaylorSeries.mul!(tmp1353[i, j, 1], S21M_t, cos_mλ[i, j, 1], ord)
                                TaylorSeries.mul!(tmp1354[i, j, 1], C21M_t, sin_mλ[i, j, 1], ord)
                                TaylorSeries.subst!(tmp1355[i, j, 1], tmp1353[i, j, 1], tmp1354[i, j, 1], ord)
                                TaylorSeries.mul!(tmp1356[i, j, 2, 1], tmp1352[i, j, 2, 1], tmp1355[i, j, 1], ord)
                                TaylorSeries.mul!(tmp1357[i, j, 2, 2], secϕ_P_nm[i, j, 2, 2], lnm7[2], ord)
                                TaylorSeries.mul!(tmp1358[i, j, 2], S22M_t, cos_mλ[i, j, 2], ord)
                                TaylorSeries.mul!(tmp1359[i, j, 2], C22M_t, sin_mλ[i, j, 2], ord)
                                TaylorSeries.subst!(tmp1360[i, j, 2], tmp1358[i, j, 2], tmp1359[i, j, 2], ord)
                                TaylorSeries.mul!(tmp1361[i, j, 2, 2], tmp1357[i, j, 2, 2], tmp1360[i, j, 2], ord)
                                TaylorSeries.add!(tmp1362[i, j, 2, 1], tmp1356[i, j, 2, 1], tmp1361[i, j, 2, 2], ord)
                                TaylorSeries.div!(F_CS_η[i, j], tmp1362[i, j, 2, 1], r_p4[i, j], ord)
                                TaylorSeries.mul!(tmp1364[i, j, 1], C21M_t, cos_mλ[i, j, 1], ord)
                                TaylorSeries.mul!(tmp1365[i, j, 1], S21M_t, sin_mλ[i, j, 1], ord)
                                TaylorSeries.add!(tmp1366[i, j, 1], tmp1364[i, j, 1], tmp1365[i, j, 1], ord)
                                TaylorSeries.mul!(tmp1367[i, j, 2, 1], cosϕ_dP_nm[i, j, 2, 1], tmp1366[i, j, 1], ord)
                                TaylorSeries.mul!(tmp1368[i, j, 2], C22M_t, cos_mλ[i, j, 2], ord)
                                TaylorSeries.mul!(tmp1369[i, j, 2], S22M_t, sin_mλ[i, j, 2], ord)
                                TaylorSeries.add!(tmp1370[i, j, 2], tmp1368[i, j, 2], tmp1369[i, j, 2], ord)
                                TaylorSeries.mul!(tmp1371[i, j, 2, 2], cosϕ_dP_nm[i, j, 2, 2], tmp1370[i, j, 2], ord)
                                TaylorSeries.add!(tmp1372[i, j, 2, 1], tmp1367[i, j, 2, 1], tmp1371[i, j, 2, 2], ord)
                                TaylorSeries.div!(F_CS_ζ[i, j], tmp1372[i, j, 2, 1], r_p4[i, j], ord)
                                TaylorSeries.identity!(F_CS_ξ_36[i, j], zero_q_1, ord)
                                TaylorSeries.identity!(F_CS_η_36[i, j], zero_q_1, ord)
                                TaylorSeries.identity!(F_CS_ζ_36[i, j], zero_q_1, ord)
                                for n = 3:n1SEM[mo]
                                    for m = 1:n
                                        TaylorSeries.mul!(tmp1374[i, j, n, m], P_nm[i, j, n, m], lnm6[n], ord)
                                        TaylorSeries.mul!(tmp1375[i, j, m], cos_mλ[i, j, m], CM[n, m], ord)
                                        TaylorSeries.mul!(tmp1376[i, j, m], sin_mλ[i, j, m], SM[n, m], ord)
                                        TaylorSeries.add!(tmp1377[i, j, m], tmp1375[i, j, m], tmp1376[i, j, m], ord)
                                        TaylorSeries.mul!(tmp1378[i, j, n, m], tmp1374[i, j, n, m], tmp1377[i, j, m], ord)
                                        TaylorSeries.div!(tmp1379[i, j, n, m], tmp1378[i, j, n, m], temp_rn[i, j, n], ord)
                                        TaylorSeries.add!(temp_CS_ξ[i, j, n, m], tmp1379[i, j, n, m], F_CS_ξ_36[i, j], ord)
                                        TaylorSeries.mul!(tmp1381[i, j, n, m], secϕ_P_nm[i, j, n, m], lnm7[m], ord)
                                        TaylorSeries.mul!(tmp1382[i, j, m], cos_mλ[i, j, m], SM[n, m], ord)
                                        TaylorSeries.mul!(tmp1383[i, j, m], sin_mλ[i, j, m], CM[n, m], ord)
                                        TaylorSeries.subst!(tmp1384[i, j, m], tmp1382[i, j, m], tmp1383[i, j, m], ord)
                                        TaylorSeries.mul!(tmp1385[i, j, n, m], tmp1381[i, j, n, m], tmp1384[i, j, m], ord)
                                        TaylorSeries.div!(tmp1386[i, j, n, m], tmp1385[i, j, n, m], temp_rn[i, j, n], ord)
                                        TaylorSeries.add!(temp_CS_η[i, j, n, m], tmp1386[i, j, n, m], F_CS_η_36[i, j], ord)
                                        TaylorSeries.mul!(tmp1388[i, j, m], cos_mλ[i, j, m], CM[n, m], ord)
                                        TaylorSeries.mul!(tmp1389[i, j, m], sin_mλ[i, j, m], SM[n, m], ord)
                                        TaylorSeries.add!(tmp1390[i, j, m], tmp1388[i, j, m], tmp1389[i, j, m], ord)
                                        TaylorSeries.mul!(tmp1391[i, j, n, m], cosϕ_dP_nm[i, j, n, m], tmp1390[i, j, m], ord)
                                        TaylorSeries.div!(tmp1392[i, j, n, m], tmp1391[i, j, n, m], temp_rn[i, j, n], ord)
                                        TaylorSeries.add!(temp_CS_ζ[i, j, n, m], tmp1392[i, j, n, m], F_CS_ζ_36[i, j], ord)
                                        TaylorSeries.identity!(F_CS_ξ_36[i, j], temp_CS_ξ[i, j, n, m], ord)
                                        TaylorSeries.identity!(F_CS_η_36[i, j], temp_CS_η[i, j, n, m], ord)
                                        TaylorSeries.identity!(F_CS_ζ_36[i, j], temp_CS_ζ[i, j, n, m], ord)
                                    end
                                end
                                TaylorSeries.add!(tmp1394[i, j], F_J_ξ[i, j], F_J_ξ_36[i, j], ord)
                                TaylorSeries.add!(tmp1395[i, j], F_CS_ξ[i, j], F_CS_ξ_36[i, j], ord)
                                TaylorSeries.add!(F_JCS_ξ[i, j], tmp1394[i, j], tmp1395[i, j], ord)
                                TaylorSeries.add!(F_JCS_η[i, j], F_CS_η[i, j], F_CS_η_36[i, j], ord)
                                TaylorSeries.add!(tmp1398[i, j], F_J_ζ[i, j], F_J_ζ_36[i, j], ord)
                                TaylorSeries.add!(tmp1399[i, j], F_CS_ζ[i, j], F_CS_ζ_36[i, j], ord)
                                TaylorSeries.add!(F_JCS_ζ[i, j], tmp1398[i, j], tmp1399[i, j], ord)
                            else
                                TaylorSeries.add!(F_JCS_ξ[i, j], F_J_ξ[i, j], F_J_ξ_36[i, j], ord)
                                TaylorSeries.identity!(F_JCS_η[i, j], zero_q_1, ord)
                                TaylorSeries.add!(F_JCS_ζ[i, j], F_J_ζ[i, j], F_J_ζ_36[i, j], ord)
                            end
                            TaylorSeries.mul!(Rb2p[i, j, 1, 1], cos_ϕ[i, j], cos_λ[i, j], ord)
                            TaylorSeries.subst!(Rb2p[i, j, 2, 1], sin_λ[i, j], ord)
                            TaylorSeries.subst!(tmp1405[i, j], sin_ϕ[i, j], ord)
                            TaylorSeries.mul!(Rb2p[i, j, 3, 1], tmp1405[i, j], cos_λ[i, j], ord)
                            TaylorSeries.mul!(Rb2p[i, j, 1, 2], cos_ϕ[i, j], sin_λ[i, j], ord)
                            TaylorSeries.identity!(Rb2p[i, j, 2, 2], cos_λ[i, j], ord)
                            TaylorSeries.subst!(tmp1408[i, j], sin_ϕ[i, j], ord)
                            TaylorSeries.mul!(Rb2p[i, j, 3, 2], tmp1408[i, j], sin_λ[i, j], ord)
                            TaylorSeries.identity!(Rb2p[i, j, 1, 3], sin_ϕ[i, j], ord)
                            TaylorSeries.identity!(Rb2p[i, j, 2, 3], zero_q_1, ord)
                            TaylorSeries.identity!(Rb2p[i, j, 3, 3], cos_ϕ[i, j], ord)
                            TaylorSeries.mul!(tmp1410[i, j, 1, 1], Rb2p[i, j, 1, 1], M_[1, 1, j], ord)
                            TaylorSeries.mul!(tmp1411[i, j, 1, 2], Rb2p[i, j, 1, 2], M_[2, 1, j], ord)
                            TaylorSeries.add!(tmp1412[i, j, 1, 1], tmp1410[i, j, 1, 1], tmp1411[i, j, 1, 2], ord)
                            TaylorSeries.mul!(tmp1413[i, j, 1, 3], Rb2p[i, j, 1, 3], M_[3, 1, j], ord)
                            TaylorSeries.add!(Gc2p[i, j, 1, 1], tmp1412[i, j, 1, 1], tmp1413[i, j, 1, 3], ord)
                            TaylorSeries.mul!(tmp1415[i, j, 2, 1], Rb2p[i, j, 2, 1], M_[1, 1, j], ord)
                            TaylorSeries.mul!(tmp1416[i, j, 2, 2], Rb2p[i, j, 2, 2], M_[2, 1, j], ord)
                            TaylorSeries.add!(tmp1417[i, j, 2, 1], tmp1415[i, j, 2, 1], tmp1416[i, j, 2, 2], ord)
                            TaylorSeries.mul!(tmp1418[i, j, 2, 3], Rb2p[i, j, 2, 3], M_[3, 1, j], ord)
                            TaylorSeries.add!(Gc2p[i, j, 2, 1], tmp1417[i, j, 2, 1], tmp1418[i, j, 2, 3], ord)
                            TaylorSeries.mul!(tmp1420[i, j, 3, 1], Rb2p[i, j, 3, 1], M_[1, 1, j], ord)
                            TaylorSeries.mul!(tmp1421[i, j, 3, 2], Rb2p[i, j, 3, 2], M_[2, 1, j], ord)
                            TaylorSeries.add!(tmp1422[i, j, 3, 1], tmp1420[i, j, 3, 1], tmp1421[i, j, 3, 2], ord)
                            TaylorSeries.mul!(tmp1423[i, j, 3, 3], Rb2p[i, j, 3, 3], M_[3, 1, j], ord)
                            TaylorSeries.add!(Gc2p[i, j, 3, 1], tmp1422[i, j, 3, 1], tmp1423[i, j, 3, 3], ord)
                            TaylorSeries.mul!(tmp1425[i, j, 1, 1], Rb2p[i, j, 1, 1], M_[1, 2, j], ord)
                            TaylorSeries.mul!(tmp1426[i, j, 1, 2], Rb2p[i, j, 1, 2], M_[2, 2, j], ord)
                            TaylorSeries.add!(tmp1427[i, j, 1, 1], tmp1425[i, j, 1, 1], tmp1426[i, j, 1, 2], ord)
                            TaylorSeries.mul!(tmp1428[i, j, 1, 3], Rb2p[i, j, 1, 3], M_[3, 2, j], ord)
                            TaylorSeries.add!(Gc2p[i, j, 1, 2], tmp1427[i, j, 1, 1], tmp1428[i, j, 1, 3], ord)
                            TaylorSeries.mul!(tmp1430[i, j, 2, 1], Rb2p[i, j, 2, 1], M_[1, 2, j], ord)
                            TaylorSeries.mul!(tmp1431[i, j, 2, 2], Rb2p[i, j, 2, 2], M_[2, 2, j], ord)
                            TaylorSeries.add!(tmp1432[i, j, 2, 1], tmp1430[i, j, 2, 1], tmp1431[i, j, 2, 2], ord)
                            TaylorSeries.mul!(tmp1433[i, j, 2, 3], Rb2p[i, j, 2, 3], M_[3, 2, j], ord)
                            TaylorSeries.add!(Gc2p[i, j, 2, 2], tmp1432[i, j, 2, 1], tmp1433[i, j, 2, 3], ord)
                            TaylorSeries.mul!(tmp1435[i, j, 3, 1], Rb2p[i, j, 3, 1], M_[1, 2, j], ord)
                            TaylorSeries.mul!(tmp1436[i, j, 3, 2], Rb2p[i, j, 3, 2], M_[2, 2, j], ord)
                            TaylorSeries.add!(tmp1437[i, j, 3, 1], tmp1435[i, j, 3, 1], tmp1436[i, j, 3, 2], ord)
                            TaylorSeries.mul!(tmp1438[i, j, 3, 3], Rb2p[i, j, 3, 3], M_[3, 2, j], ord)
                            TaylorSeries.add!(Gc2p[i, j, 3, 2], tmp1437[i, j, 3, 1], tmp1438[i, j, 3, 3], ord)
                            TaylorSeries.mul!(tmp1440[i, j, 1, 1], Rb2p[i, j, 1, 1], M_[1, 3, j], ord)
                            TaylorSeries.mul!(tmp1441[i, j, 1, 2], Rb2p[i, j, 1, 2], M_[2, 3, j], ord)
                            TaylorSeries.add!(tmp1442[i, j, 1, 1], tmp1440[i, j, 1, 1], tmp1441[i, j, 1, 2], ord)
                            TaylorSeries.mul!(tmp1443[i, j, 1, 3], Rb2p[i, j, 1, 3], M_[3, 3, j], ord)
                            TaylorSeries.add!(Gc2p[i, j, 1, 3], tmp1442[i, j, 1, 1], tmp1443[i, j, 1, 3], ord)
                            TaylorSeries.mul!(tmp1445[i, j, 2, 1], Rb2p[i, j, 2, 1], M_[1, 3, j], ord)
                            TaylorSeries.mul!(tmp1446[i, j, 2, 2], Rb2p[i, j, 2, 2], M_[2, 3, j], ord)
                            TaylorSeries.add!(tmp1447[i, j, 2, 1], tmp1445[i, j, 2, 1], tmp1446[i, j, 2, 2], ord)
                            TaylorSeries.mul!(tmp1448[i, j, 2, 3], Rb2p[i, j, 2, 3], M_[3, 3, j], ord)
                            TaylorSeries.add!(Gc2p[i, j, 2, 3], tmp1447[i, j, 2, 1], tmp1448[i, j, 2, 3], ord)
                            TaylorSeries.mul!(tmp1450[i, j, 3, 1], Rb2p[i, j, 3, 1], M_[1, 3, j], ord)
                            TaylorSeries.mul!(tmp1451[i, j, 3, 2], Rb2p[i, j, 3, 2], M_[2, 3, j], ord)
                            TaylorSeries.add!(tmp1452[i, j, 3, 1], tmp1450[i, j, 3, 1], tmp1451[i, j, 3, 2], ord)
                            TaylorSeries.mul!(tmp1453[i, j, 3, 3], Rb2p[i, j, 3, 3], M_[3, 3, j], ord)
                            TaylorSeries.add!(Gc2p[i, j, 3, 3], tmp1452[i, j, 3, 1], tmp1453[i, j, 3, 3], ord)
                            TaylorSeries.mul!(tmp1455[i, j, 1, 1], F_JCS_ξ[i, j], Gc2p[i, j, 1, 1], ord)
                            TaylorSeries.mul!(tmp1456[i, j, 2, 1], F_JCS_η[i, j], Gc2p[i, j, 2, 1], ord)
                            TaylorSeries.add!(tmp1457[i, j, 1, 1], tmp1455[i, j, 1, 1], tmp1456[i, j, 2, 1], ord)
                            TaylorSeries.mul!(tmp1458[i, j, 3, 1], F_JCS_ζ[i, j], Gc2p[i, j, 3, 1], ord)
                            TaylorSeries.add!(F_JCS_x[i, j], tmp1457[i, j, 1, 1], tmp1458[i, j, 3, 1], ord)
                            TaylorSeries.mul!(tmp1460[i, j, 1, 2], F_JCS_ξ[i, j], Gc2p[i, j, 1, 2], ord)
                            TaylorSeries.mul!(tmp1461[i, j, 2, 2], F_JCS_η[i, j], Gc2p[i, j, 2, 2], ord)
                            TaylorSeries.add!(tmp1462[i, j, 1, 2], tmp1460[i, j, 1, 2], tmp1461[i, j, 2, 2], ord)
                            TaylorSeries.mul!(tmp1463[i, j, 3, 2], F_JCS_ζ[i, j], Gc2p[i, j, 3, 2], ord)
                            TaylorSeries.add!(F_JCS_y[i, j], tmp1462[i, j, 1, 2], tmp1463[i, j, 3, 2], ord)
                            TaylorSeries.mul!(tmp1465[i, j, 1, 3], F_JCS_ξ[i, j], Gc2p[i, j, 1, 3], ord)
                            TaylorSeries.mul!(tmp1466[i, j, 2, 3], F_JCS_η[i, j], Gc2p[i, j, 2, 3], ord)
                            TaylorSeries.add!(tmp1467[i, j, 1, 3], tmp1465[i, j, 1, 3], tmp1466[i, j, 2, 3], ord)
                            TaylorSeries.mul!(tmp1468[i, j, 3, 3], F_JCS_ζ[i, j], Gc2p[i, j, 3, 3], ord)
                            TaylorSeries.add!(F_JCS_z[i, j], tmp1467[i, j, 1, 3], tmp1468[i, j, 3, 3], ord)
                        end
                    end
                end
            end
        for j = 1:N_ext
            for i = 1:N_ext
                if i == j
                    continue
                else
                    if UJ_interaction[i, j]
                        TaylorSeries.mul!(tmp1470[i, j], μ[i], F_JCS_x[i, j], ord)
                        TaylorSeries.subst!(temp_accX_j[i, j], accX[j], tmp1470[i, j], ord)
                        TaylorSeries.identity!(accX[j], temp_accX_j[i, j], ord)
                        TaylorSeries.mul!(tmp1472[i, j], μ[i], F_JCS_y[i, j], ord)
                        TaylorSeries.subst!(temp_accY_j[i, j], accY[j], tmp1472[i, j], ord)
                        TaylorSeries.identity!(accY[j], temp_accY_j[i, j], ord)
                        TaylorSeries.mul!(tmp1474[i, j], μ[i], F_JCS_z[i, j], ord)
                        TaylorSeries.subst!(temp_accZ_j[i, j], accZ[j], tmp1474[i, j], ord)
                        TaylorSeries.identity!(accZ[j], temp_accZ_j[i, j], ord)
                        TaylorSeries.mul!(tmp1476[i, j], μ[j], F_JCS_x[i, j], ord)
                        TaylorSeries.add!(temp_accX_i[i, j], accX[i], tmp1476[i, j], ord)
                        TaylorSeries.identity!(accX[i], temp_accX_i[i, j], ord)
                        TaylorSeries.mul!(tmp1478[i, j], μ[j], F_JCS_y[i, j], ord)
                        TaylorSeries.add!(temp_accY_i[i, j], accY[i], tmp1478[i, j], ord)
                        TaylorSeries.identity!(accY[i], temp_accY_i[i, j], ord)
                        TaylorSeries.mul!(tmp1480[i, j], μ[j], F_JCS_z[i, j], ord)
                        TaylorSeries.add!(temp_accZ_i[i, j], accZ[i], tmp1480[i, j], ord)
                        TaylorSeries.identity!(accZ[i], temp_accZ_i[i, j], ord)
                    end
                end
            end
        end
        #= REPL[2]:489 =# Threads.@threads for j = 1:N
                for i = 1:N
                    if i == j
                        continue
                    else
                        TaylorSeries.mul!(_4ϕj[i, j], 4, newtonianNb_Potential[j], ord)
                        TaylorSeries.add!(ϕi_plus_4ϕj[i, j], newtonianNb_Potential[i], _4ϕj[i, j], ord)
                        TaylorSeries.mul!(tmp1486[i], 2, v2[i], ord)
                        TaylorSeries.add!(tmp1487[j], v2[j], tmp1486[i], ord)
                        TaylorSeries.mul!(tmp1489[i, j], 4, vi_dot_vj[i, j], ord)
                        TaylorSeries.subst!(sj2_plus_2si2_minus_4vivj[i, j], tmp1487[j], tmp1489[i, j], ord)
                        TaylorSeries.subst!(ϕs_and_vs[i, j], sj2_plus_2si2_minus_4vivj[i, j], ϕi_plus_4ϕj[i, j], ord)
                        TaylorSeries.mul!(Xij_t_Ui[i, j], X[i, j], dq[3i - 2], ord)
                        TaylorSeries.mul!(Yij_t_Vi[i, j], Y[i, j], dq[3i - 1], ord)
                        TaylorSeries.mul!(Zij_t_Wi[i, j], Z[i, j], dq[3i], ord)
                        TaylorSeries.add!(tmp1495[i, j], Xij_t_Ui[i, j], Yij_t_Vi[i, j], ord)
                        TaylorSeries.add!(Rij_dot_Vi[i, j], tmp1495[i, j], Zij_t_Wi[i, j], ord)
                        TaylorSeries.pow!(tmp1498[i, j], Rij_dot_Vi[i, j], 2, ord)
                        TaylorSeries.div!(pn1t7[i, j], tmp1498[i, j], r_p2[i, j], ord)
                        TaylorSeries.mul!(tmp1501[i, j], 1.5, pn1t7[i, j], ord)
                        TaylorSeries.subst!(pn1t2_7[i, j], ϕs_and_vs[i, j], tmp1501[i, j], ord)
                        TaylorSeries.add!(pn1t1_7[i, j], c_p2, pn1t2_7[i, j], ord)
                        for k = 1:postnewton_iter
                            TaylorSeries.identity!(pn1[i, j, k], zero_q_1, ord)
                            TaylorSeries.identity!(X_t_pn1[i, j, k], zero_q_1, ord)
                            TaylorSeries.identity!(Y_t_pn1[i, j, k], zero_q_1, ord)
                            TaylorSeries.identity!(Z_t_pn1[i, j, k], zero_q_1, ord)
                            TaylorSeries.identity!(pNX_t_pn3[i, j, k], zero_q_1, ord)
                            TaylorSeries.identity!(pNY_t_pn3[i, j, k], zero_q_1, ord)
                            TaylorSeries.identity!(pNZ_t_pn3[i, j, k], zero_q_1, ord)
                            TaylorSeries.identity!(pNX_t_X[i, j, k], zero_q_1, ord)
                            TaylorSeries.identity!(pNY_t_Y[i, j, k], zero_q_1, ord)
                            TaylorSeries.identity!(pNZ_t_Z[i, j, k], zero_q_1, ord)
                        end
                    end
                end
                TaylorSeries.identity!(postNewtonX[j, 1], newtonX[j], ord)
                TaylorSeries.identity!(postNewtonY[j, 1], newtonY[j], ord)
                TaylorSeries.identity!(postNewtonZ[j, 1], newtonZ[j], ord)
                for k = 1:postnewton_iter
                    TaylorSeries.identity!(pntempX[j, k], zero_q_1, ord)
                    TaylorSeries.identity!(pntempY[j, k], zero_q_1, ord)
                    TaylorSeries.identity!(pntempZ[j, k], zero_q_1, ord)
                end
            end
        for k = 1:postnewton_iter
            #= REPL[2]:534 =# Threads.@threads for j = 1:N
                    for i = 1:N
                        if i == j
                            continue
                        else
                            TaylorSeries.mul!(pNX_t_X[i, j, k], postNewtonX[i, k], X[i, j], ord)
                            TaylorSeries.mul!(pNY_t_Y[i, j, k], postNewtonY[i, k], Y[i, j], ord)
                            TaylorSeries.mul!(pNZ_t_Z[i, j, k], postNewtonZ[i, k], Z[i, j], ord)
                            TaylorSeries.add!(tmp1508[i, j, k], pNX_t_X[i, j, k], pNY_t_Y[i, j, k], ord)
                            TaylorSeries.add!(tmp1509[i, j, k], tmp1508[i, j, k], pNZ_t_Z[i, j, k], ord)
                            TaylorSeries.mul!(tmp1510[i, j, k], 0.5, tmp1509[i, j, k], ord)
                            TaylorSeries.add!(pn1[i, j, k], pn1t1_7[i, j], tmp1510[i, j, k], ord)
                            TaylorSeries.mul!(X_t_pn1[i, j, k], newton_acc_X[i, j], pn1[i, j, k], ord)
                            TaylorSeries.mul!(Y_t_pn1[i, j, k], newton_acc_Y[i, j], pn1[i, j, k], ord)
                            TaylorSeries.mul!(Z_t_pn1[i, j, k], newton_acc_Z[i, j], pn1[i, j, k], ord)
                            TaylorSeries.mul!(pNX_t_pn3[i, j, k], postNewtonX[i, k], pn3[i, j], ord)
                            TaylorSeries.mul!(pNY_t_pn3[i, j, k], postNewtonY[i, k], pn3[i, j], ord)
                            TaylorSeries.mul!(pNZ_t_pn3[i, j, k], postNewtonZ[i, k], pn3[i, j], ord)
                            TaylorSeries.add!(tmp1518[i, j, k], U_t_pn2[i, j], pNX_t_pn3[i, j, k], ord)
                            TaylorSeries.add!(termpnx[i, j, k], X_t_pn1[i, j, k], tmp1518[i, j, k], ord)
                            TaylorSeries.add!(sumpnx[i, j, k], pntempX[j, k], termpnx[i, j, k], ord)
                            TaylorSeries.identity!(pntempX[j, k], sumpnx[i, j, k], ord)
                            TaylorSeries.add!(tmp1521[i, j, k], V_t_pn2[i, j], pNY_t_pn3[i, j, k], ord)
                            TaylorSeries.add!(termpny[i, j, k], Y_t_pn1[i, j, k], tmp1521[i, j, k], ord)
                            TaylorSeries.add!(sumpny[i, j, k], pntempY[j, k], termpny[i, j, k], ord)
                            TaylorSeries.identity!(pntempY[j, k], sumpny[i, j, k], ord)
                            TaylorSeries.add!(tmp1524[i, j, k], W_t_pn2[i, j], pNZ_t_pn3[i, j, k], ord)
                            TaylorSeries.add!(termpnz[i, j, k], Z_t_pn1[i, j, k], tmp1524[i, j, k], ord)
                            TaylorSeries.add!(sumpnz[i, j, k], pntempZ[j, k], termpnz[i, j, k], ord)
                            TaylorSeries.identity!(pntempZ[j, k], sumpnz[i, j, k], ord)
                        end
                    end
                    TaylorSeries.mul!(postNewtonX[j, k + 1], pntempX[j, k], c_m2, ord)
                    TaylorSeries.mul!(postNewtonY[j, k + 1], pntempY[j, k], c_m2, ord)
                    TaylorSeries.mul!(postNewtonZ[j, k + 1], pntempZ[j, k], c_m2, ord)
                end
        end
        TaylorSeries.identity!(X_E_τ_0, q_del_τ_0[3ea - 2], ord)
        TaylorSeries.identity!(Y_E_τ_0, q_del_τ_0[3ea - 1], ord)
        TaylorSeries.identity!(Z_E_τ_0, q_del_τ_0[3ea], ord)
        TaylorSeries.identity!(X_E_τ_1, q_del_τ_1[3ea - 2], ord)
        TaylorSeries.identity!(Y_E_τ_1, q_del_τ_1[3ea - 1], ord)
        TaylorSeries.identity!(Z_E_τ_1, q_del_τ_1[3ea], ord)
        TaylorSeries.identity!(X_E_τ_2, q_del_τ_2[3ea - 2], ord)
        TaylorSeries.identity!(Y_E_τ_2, q_del_τ_2[3ea - 1], ord)
        TaylorSeries.identity!(Z_E_τ_2, q_del_τ_2[3ea], ord)
        TaylorSeries.subst!(X_ME_τ_0, q_del_τ_0[3mo - 2], X_E_τ_0, ord)
        TaylorSeries.subst!(Y_ME_τ_0, q_del_τ_0[3mo - 1], Y_E_τ_0, ord)
        TaylorSeries.subst!(Z_ME_τ_0, q_del_τ_0[3mo], Z_E_τ_0, ord)
        TaylorSeries.subst!(X_ME_τ_1, q_del_τ_1[3mo - 2], X_E_τ_1, ord)
        TaylorSeries.subst!(Y_ME_τ_1, q_del_τ_1[3mo - 1], Y_E_τ_1, ord)
        TaylorSeries.subst!(Z_ME_τ_1, q_del_τ_1[3mo], Z_E_τ_1, ord)
        TaylorSeries.subst!(X_ME_τ_2, q_del_τ_2[3mo - 2], X_E_τ_2, ord)
        TaylorSeries.subst!(Y_ME_τ_2, q_del_τ_2[3mo - 1], Y_E_τ_2, ord)
        TaylorSeries.subst!(Z_ME_τ_2, q_del_τ_2[3mo], Z_E_τ_2, ord)
        TaylorSeries.subst!(X_SE_τ_0, q_del_τ_0[3su - 2], X_E_τ_0, ord)
        TaylorSeries.subst!(Y_SE_τ_0, q_del_τ_0[3su - 1], Y_E_τ_0, ord)
        TaylorSeries.subst!(Z_SE_τ_0, q_del_τ_0[3su], Z_E_τ_0, ord)
        TaylorSeries.subst!(X_SE_τ_1, q_del_τ_1[3su - 2], X_E_τ_1, ord)
        TaylorSeries.subst!(Y_SE_τ_1, q_del_τ_1[3su - 1], Y_E_τ_1, ord)
        TaylorSeries.subst!(Z_SE_τ_1, q_del_τ_1[3su], Z_E_τ_1, ord)
        TaylorSeries.subst!(X_SE_τ_2, q_del_τ_2[3su - 2], X_E_τ_2, ord)
        TaylorSeries.subst!(Y_SE_τ_2, q_del_τ_2[3su - 1], Y_E_τ_2, ord)
        TaylorSeries.subst!(Z_SE_τ_2, q_del_τ_2[3su], Z_E_τ_2, ord)
        TaylorSeries.mul!(tmp1548, R30[1, 1], X_ME_τ_0, ord)
        TaylorSeries.mul!(tmp1549, R30[1, 2], Y_ME_τ_0, ord)
        TaylorSeries.add!(tmp1550, tmp1548, tmp1549, ord)
        TaylorSeries.mul!(tmp1551, R30[1, 3], Z_ME_τ_0, ord)
        TaylorSeries.add!(r_star_M_0[1], tmp1550, tmp1551, ord)
        TaylorSeries.mul!(tmp1553, R30[2, 1], X_ME_τ_0, ord)
        TaylorSeries.mul!(tmp1554, R30[2, 2], Y_ME_τ_0, ord)
        TaylorSeries.add!(tmp1555, tmp1553, tmp1554, ord)
        TaylorSeries.mul!(tmp1556, R30[2, 3], Z_ME_τ_0, ord)
        TaylorSeries.add!(r_star_M_0[2], tmp1555, tmp1556, ord)
        TaylorSeries.mul!(tmp1558, R30[3, 1], X_ME_τ_0, ord)
        TaylorSeries.mul!(tmp1559, R30[3, 2], Y_ME_τ_0, ord)
        TaylorSeries.add!(tmp1560, tmp1558, tmp1559, ord)
        TaylorSeries.mul!(tmp1561, R30[3, 3], Z_ME_τ_0, ord)
        TaylorSeries.add!(r_star_M_0[3], tmp1560, tmp1561, ord)
        TaylorSeries.mul!(tmp1563, R31[1, 1], X_ME_τ_1, ord)
        TaylorSeries.mul!(tmp1564, R31[1, 2], Y_ME_τ_1, ord)
        TaylorSeries.add!(tmp1565, tmp1563, tmp1564, ord)
        TaylorSeries.mul!(tmp1566, R31[1, 3], Z_ME_τ_1, ord)
        TaylorSeries.add!(r_star_M_1[1], tmp1565, tmp1566, ord)
        TaylorSeries.mul!(tmp1568, R31[2, 1], X_ME_τ_1, ord)
        TaylorSeries.mul!(tmp1569, R31[2, 2], Y_ME_τ_1, ord)
        TaylorSeries.add!(tmp1570, tmp1568, tmp1569, ord)
        TaylorSeries.mul!(tmp1571, R31[2, 3], Z_ME_τ_1, ord)
        TaylorSeries.add!(r_star_M_1[2], tmp1570, tmp1571, ord)
        TaylorSeries.mul!(tmp1573, R31[3, 1], X_ME_τ_1, ord)
        TaylorSeries.mul!(tmp1574, R31[3, 2], Y_ME_τ_1, ord)
        TaylorSeries.add!(tmp1575, tmp1573, tmp1574, ord)
        TaylorSeries.mul!(tmp1576, R31[3, 3], Z_ME_τ_1, ord)
        TaylorSeries.add!(r_star_M_1[3], tmp1575, tmp1576, ord)
        TaylorSeries.mul!(tmp1578, R32[1, 1], X_ME_τ_2, ord)
        TaylorSeries.mul!(tmp1579, R32[1, 2], Y_ME_τ_2, ord)
        TaylorSeries.add!(tmp1580, tmp1578, tmp1579, ord)
        TaylorSeries.mul!(tmp1581, R32[1, 3], Z_ME_τ_2, ord)
        TaylorSeries.add!(r_star_M_2[1], tmp1580, tmp1581, ord)
        TaylorSeries.mul!(tmp1583, R32[2, 1], X_ME_τ_2, ord)
        TaylorSeries.mul!(tmp1584, R32[2, 2], Y_ME_τ_2, ord)
        TaylorSeries.add!(tmp1585, tmp1583, tmp1584, ord)
        TaylorSeries.mul!(tmp1586, R32[2, 3], Z_ME_τ_2, ord)
        TaylorSeries.add!(r_star_M_2[2], tmp1585, tmp1586, ord)
        TaylorSeries.mul!(tmp1588, R32[3, 1], X_ME_τ_2, ord)
        TaylorSeries.mul!(tmp1589, R32[3, 2], Y_ME_τ_2, ord)
        TaylorSeries.add!(tmp1590, tmp1588, tmp1589, ord)
        TaylorSeries.mul!(tmp1591, R32[3, 3], Z_ME_τ_2, ord)
        TaylorSeries.add!(r_star_M_2[3], tmp1590, tmp1591, ord)
        TaylorSeries.mul!(tmp1593, R30[1, 1], X_SE_τ_0, ord)
        TaylorSeries.mul!(tmp1594, R30[1, 2], Y_SE_τ_0, ord)
        TaylorSeries.add!(tmp1595, tmp1593, tmp1594, ord)
        TaylorSeries.mul!(tmp1596, R30[1, 3], Z_SE_τ_0, ord)
        TaylorSeries.add!(r_star_S_0[1], tmp1595, tmp1596, ord)
        TaylorSeries.mul!(tmp1598, R30[2, 1], X_SE_τ_0, ord)
        TaylorSeries.mul!(tmp1599, R30[2, 2], Y_SE_τ_0, ord)
        TaylorSeries.add!(tmp1600, tmp1598, tmp1599, ord)
        TaylorSeries.mul!(tmp1601, R30[2, 3], Z_SE_τ_0, ord)
        TaylorSeries.add!(r_star_S_0[2], tmp1600, tmp1601, ord)
        TaylorSeries.mul!(tmp1603, R30[3, 1], X_SE_τ_0, ord)
        TaylorSeries.mul!(tmp1604, R30[3, 2], Y_SE_τ_0, ord)
        TaylorSeries.add!(tmp1605, tmp1603, tmp1604, ord)
        TaylorSeries.mul!(tmp1606, R30[3, 3], Z_SE_τ_0, ord)
        TaylorSeries.add!(r_star_S_0[3], tmp1605, tmp1606, ord)
        TaylorSeries.mul!(tmp1608, R31[1, 1], X_SE_τ_1, ord)
        TaylorSeries.mul!(tmp1609, R31[1, 2], Y_SE_τ_1, ord)
        TaylorSeries.add!(tmp1610, tmp1608, tmp1609, ord)
        TaylorSeries.mul!(tmp1611, R31[1, 3], Z_SE_τ_1, ord)
        TaylorSeries.add!(r_star_S_1[1], tmp1610, tmp1611, ord)
        TaylorSeries.mul!(tmp1613, R31[2, 1], X_SE_τ_1, ord)
        TaylorSeries.mul!(tmp1614, R31[2, 2], Y_SE_τ_1, ord)
        TaylorSeries.add!(tmp1615, tmp1613, tmp1614, ord)
        TaylorSeries.mul!(tmp1616, R31[2, 3], Z_SE_τ_1, ord)
        TaylorSeries.add!(r_star_S_1[2], tmp1615, tmp1616, ord)
        TaylorSeries.mul!(tmp1618, R31[3, 1], X_SE_τ_1, ord)
        TaylorSeries.mul!(tmp1619, R31[3, 2], Y_SE_τ_1, ord)
        TaylorSeries.add!(tmp1620, tmp1618, tmp1619, ord)
        TaylorSeries.mul!(tmp1621, R31[3, 3], Z_SE_τ_1, ord)
        TaylorSeries.add!(r_star_S_1[3], tmp1620, tmp1621, ord)
        TaylorSeries.mul!(tmp1623, R32[1, 1], X_SE_τ_2, ord)
        TaylorSeries.mul!(tmp1624, R32[1, 2], Y_SE_τ_2, ord)
        TaylorSeries.add!(tmp1625, tmp1623, tmp1624, ord)
        TaylorSeries.mul!(tmp1626, R32[1, 3], Z_SE_τ_2, ord)
        TaylorSeries.add!(r_star_S_2[1], tmp1625, tmp1626, ord)
        TaylorSeries.mul!(tmp1628, R32[2, 1], X_SE_τ_2, ord)
        TaylorSeries.mul!(tmp1629, R32[2, 2], Y_SE_τ_2, ord)
        TaylorSeries.add!(tmp1630, tmp1628, tmp1629, ord)
        TaylorSeries.mul!(tmp1631, R32[2, 3], Z_SE_τ_2, ord)
        TaylorSeries.add!(r_star_S_2[2], tmp1630, tmp1631, ord)
        TaylorSeries.mul!(tmp1633, R32[3, 1], X_SE_τ_2, ord)
        TaylorSeries.mul!(tmp1634, R32[3, 2], Y_SE_τ_2, ord)
        TaylorSeries.add!(tmp1635, tmp1633, tmp1634, ord)
        TaylorSeries.mul!(tmp1636, R32[3, 3], Z_SE_τ_2, ord)
        TaylorSeries.add!(r_star_S_2[3], tmp1635, tmp1636, ord)
        TaylorSeries.pow!(tmp1639, r_star_M_0[1], 2, ord)
        TaylorSeries.pow!(tmp1641, r_star_M_0[2], 2, ord)
        TaylorSeries.add!(ρ0s2_M, tmp1639, tmp1641, ord)
        TaylorSeries.sqrt!(ρ0s_M, ρ0s2_M, ord)
        TaylorSeries.pow!(z0s2_M, r_star_M_0[3], 2, ord)
        TaylorSeries.add!(r0s2_M, ρ0s2_M, z0s2_M, ord)
        TaylorSeries.sqrt!(r0s_M, r0s2_M, ord)
        TaylorSeries.pow!(r0s5_M, r0s_M, 5, ord)
        TaylorSeries.pow!(tmp1651, r_star_S_0[1], 2, ord)
        TaylorSeries.pow!(tmp1653, r_star_S_0[2], 2, ord)
        TaylorSeries.add!(ρ0s2_S, tmp1651, tmp1653, ord)
        TaylorSeries.sqrt!(ρ0s_S, ρ0s2_S, ord)
        TaylorSeries.pow!(z0s2_S, r_star_S_0[3], 2, ord)
        TaylorSeries.add!(r0s2_S, ρ0s2_S, z0s2_S, ord)
        TaylorSeries.sqrt!(r0s_S, r0s2_S, ord)
        TaylorSeries.pow!(r0s5_S, r0s_S, 5, ord)
        TaylorSeries.mul!(tmp1663, Z_bf[mo, ea], r_star_M_0[3], ord)
        TaylorSeries.pow!(tmp1665, tmp1663, 2, ord)
        TaylorSeries.mul!(tmp1667, r_xy[mo, ea], ρ0s_M, ord)
        TaylorSeries.pow!(tmp1669, tmp1667, 2, ord)
        TaylorSeries.mul!(tmp1670, 0.5, tmp1669, ord)
        TaylorSeries.add!(tmp1671, tmp1665, tmp1670, ord)
        TaylorSeries.div!(tmp1672, tmp1671, r_p2[mo, ea], ord)
        TaylorSeries.mul!(tmp1673, 5, tmp1672, ord)
        TaylorSeries.subst!(coeff0_M, r0s2_M, tmp1673, ord)
        TaylorSeries.mul!(tmp1676, Z_bf[mo, ea], r_star_S_0[3], ord)
        TaylorSeries.pow!(tmp1678, tmp1676, 2, ord)
        TaylorSeries.mul!(tmp1680, r_xy[mo, ea], ρ0s_S, ord)
        TaylorSeries.pow!(tmp1682, tmp1680, 2, ord)
        TaylorSeries.mul!(tmp1683, 0.5, tmp1682, ord)
        TaylorSeries.add!(tmp1684, tmp1678, tmp1683, ord)
        TaylorSeries.div!(tmp1685, tmp1684, r_p2[mo, ea], ord)
        TaylorSeries.mul!(tmp1686, 5, tmp1685, ord)
        TaylorSeries.subst!(coeff0_S, r0s2_S, tmp1686, ord)
        TaylorSeries.div!(k_20E_div_r0s5_M, k_20E, r0s5_M, ord)
        TaylorSeries.div!(k_20E_div_r0s5_S, k_20E, r0s5_S, ord)
        TaylorSeries.add!(tmp1690, ρ0s2_M, coeff0_M, ord)
        TaylorSeries.mul!(tmp1691, k_20E_div_r0s5_M, tmp1690, ord)
        TaylorSeries.mul!(aux0_M_x, tmp1691, X_bf[mo, ea], ord)
        TaylorSeries.add!(tmp1693, ρ0s2_M, coeff0_M, ord)
        TaylorSeries.mul!(tmp1694, k_20E_div_r0s5_M, tmp1693, ord)
        TaylorSeries.mul!(aux0_M_y, tmp1694, Y_bf[mo, ea], ord)
        TaylorSeries.mul!(tmp1697, 2, z0s2_M, ord)
        TaylorSeries.add!(tmp1698, tmp1697, coeff0_M, ord)
        TaylorSeries.mul!(tmp1699, k_20E_div_r0s5_M, tmp1698, ord)
        TaylorSeries.mul!(aux0_M_z, tmp1699, Z_bf[mo, ea], ord)
        TaylorSeries.add!(tmp1701, ρ0s2_S, coeff0_S, ord)
        TaylorSeries.mul!(tmp1702, k_20E_div_r0s5_S, tmp1701, ord)
        TaylorSeries.mul!(aux0_S_x, tmp1702, X_bf[mo, ea], ord)
        TaylorSeries.add!(tmp1704, ρ0s2_S, coeff0_S, ord)
        TaylorSeries.mul!(tmp1705, k_20E_div_r0s5_S, tmp1704, ord)
        TaylorSeries.mul!(aux0_S_y, tmp1705, Y_bf[mo, ea], ord)
        TaylorSeries.mul!(tmp1708, 2, z0s2_S, ord)
        TaylorSeries.add!(tmp1709, tmp1708, coeff0_S, ord)
        TaylorSeries.mul!(tmp1710, k_20E_div_r0s5_S, tmp1709, ord)
        TaylorSeries.mul!(aux0_S_z, tmp1710, Z_bf[mo, ea], ord)
        TaylorSeries.pow!(tmp1713, r_star_M_1[1], 2, ord)
        TaylorSeries.pow!(tmp1715, r_star_M_1[2], 2, ord)
        TaylorSeries.add!(ρ1s2_M, tmp1713, tmp1715, ord)
        TaylorSeries.sqrt!(ρ1s_M, ρ1s2_M, ord)
        TaylorSeries.pow!(z1s2_M, r_star_M_1[3], 2, ord)
        TaylorSeries.add!(r1s2_M, ρ1s2_M, z1s2_M, ord)
        TaylorSeries.sqrt!(r1s_M, r1s2_M, ord)
        TaylorSeries.pow!(r1s5_M, r1s_M, 5, ord)
        TaylorSeries.pow!(tmp1725, r_star_S_1[1], 2, ord)
        TaylorSeries.pow!(tmp1727, r_star_S_1[2], 2, ord)
        TaylorSeries.add!(ρ1s2_S, tmp1725, tmp1727, ord)
        TaylorSeries.sqrt!(ρ1s_S, ρ1s2_S, ord)
        TaylorSeries.pow!(z1s2_S, r_star_S_1[3], 2, ord)
        TaylorSeries.add!(r1s2_S, ρ1s2_S, z1s2_S, ord)
        TaylorSeries.sqrt!(r1s_S, r1s2_S, ord)
        TaylorSeries.pow!(r1s5_S, r1s_S, 5, ord)
        TaylorSeries.mul!(tmp1736, X_bf[mo, ea], r_star_M_1[1], ord)
        TaylorSeries.mul!(tmp1737, Y_bf[mo, ea], r_star_M_1[2], ord)
        TaylorSeries.add!(coeff1_1_M, tmp1736, tmp1737, ord)
        TaylorSeries.mul!(tmp1739, X_bf[mo, ea], r_star_S_1[1], ord)
        TaylorSeries.mul!(tmp1740, Y_bf[mo, ea], r_star_S_1[2], ord)
        TaylorSeries.add!(coeff1_1_S, tmp1739, tmp1740, ord)
        TaylorSeries.mul!(coeff2_1_M, Z_bf[mo, ea], r_star_M_1[3], ord)
        TaylorSeries.mul!(coeff2_1_S, Z_bf[mo, ea], r_star_S_1[3], ord)
        TaylorSeries.mul!(tmp1745, 10, coeff1_1_M, ord)
        TaylorSeries.mul!(tmp1746, tmp1745, coeff2_1_M, ord)
        TaylorSeries.div!(coeff3_1_M, tmp1746, r_p2[mo, ea], ord)
        TaylorSeries.mul!(tmp1749, 10, coeff1_1_S, ord)
        TaylorSeries.mul!(tmp1750, tmp1749, coeff2_1_S, ord)
        TaylorSeries.div!(coeff3_1_S, tmp1750, r_p2[mo, ea], ord)
        TaylorSeries.div!(k_21E_div_r1s5_M, k_21E, r1s5_M, ord)
        TaylorSeries.div!(k_21E_div_r1s5_S, k_21E, r1s5_S, ord)
        TaylorSeries.mul!(tmp1755, 2, coeff2_1_M, ord)
        TaylorSeries.mul!(tmp1756, tmp1755, r_star_M_1[1], ord)
        TaylorSeries.mul!(tmp1757, coeff3_1_M, X_bf[mo, ea], ord)
        TaylorSeries.subst!(tmp1758, tmp1756, tmp1757, ord)
        TaylorSeries.mul!(aux1_M_x, k_21E_div_r1s5_M, tmp1758, ord)
        TaylorSeries.mul!(tmp1761, 2, coeff2_1_M, ord)
        TaylorSeries.mul!(tmp1762, tmp1761, r_star_M_1[2], ord)
        TaylorSeries.mul!(tmp1763, coeff3_1_M, Y_bf[mo, ea], ord)
        TaylorSeries.subst!(tmp1764, tmp1762, tmp1763, ord)
        TaylorSeries.mul!(aux1_M_y, k_21E_div_r1s5_M, tmp1764, ord)
        TaylorSeries.mul!(tmp1767, 2, coeff1_1_M, ord)
        TaylorSeries.mul!(tmp1768, tmp1767, r_star_M_1[3], ord)
        TaylorSeries.mul!(tmp1769, coeff3_1_M, Z_bf[mo, ea], ord)
        TaylorSeries.subst!(tmp1770, tmp1768, tmp1769, ord)
        TaylorSeries.mul!(aux1_M_z, k_21E_div_r1s5_M, tmp1770, ord)
        TaylorSeries.mul!(tmp1773, 2, coeff2_1_S, ord)
        TaylorSeries.mul!(tmp1774, tmp1773, r_star_S_1[1], ord)
        TaylorSeries.mul!(tmp1775, coeff3_1_S, X_bf[mo, ea], ord)
        TaylorSeries.subst!(tmp1776, tmp1774, tmp1775, ord)
        TaylorSeries.mul!(aux1_S_x, k_21E_div_r1s5_S, tmp1776, ord)
        TaylorSeries.mul!(tmp1779, 2, coeff2_1_S, ord)
        TaylorSeries.mul!(tmp1780, tmp1779, r_star_S_1[2], ord)
        TaylorSeries.mul!(tmp1781, coeff3_1_S, Y_bf[mo, ea], ord)
        TaylorSeries.subst!(tmp1782, tmp1780, tmp1781, ord)
        TaylorSeries.mul!(aux1_S_y, k_21E_div_r1s5_S, tmp1782, ord)
        TaylorSeries.mul!(tmp1785, 2, coeff1_1_S, ord)
        TaylorSeries.mul!(tmp1786, tmp1785, r_star_S_1[3], ord)
        TaylorSeries.mul!(tmp1787, coeff3_1_S, Z_bf[mo, ea], ord)
        TaylorSeries.subst!(tmp1788, tmp1786, tmp1787, ord)
        TaylorSeries.mul!(aux1_S_z, k_21E_div_r1s5_S, tmp1788, ord)
        TaylorSeries.pow!(tmp1791, r_star_M_2[1], 2, ord)
        TaylorSeries.pow!(tmp1793, r_star_M_2[2], 2, ord)
        TaylorSeries.add!(ρ2s2_M, tmp1791, tmp1793, ord)
        TaylorSeries.sqrt!(ρ2s_M, ρ2s2_M, ord)
        TaylorSeries.pow!(z2s2_M, r_star_M_2[3], 2, ord)
        TaylorSeries.add!(r2s2_M, ρ2s2_M, z2s2_M, ord)
        TaylorSeries.sqrt!(r2s_M, r2s2_M, ord)
        TaylorSeries.pow!(r2s5_M, r2s_M, 5, ord)
        TaylorSeries.pow!(tmp1803, r_star_S_2[1], 2, ord)
        TaylorSeries.pow!(tmp1805, r_star_S_2[2], 2, ord)
        TaylorSeries.add!(ρ2s2_S, tmp1803, tmp1805, ord)
        TaylorSeries.sqrt!(ρ2s_S, ρ2s2_S, ord)
        TaylorSeries.pow!(z2s2_S, r_star_S_2[3], 2, ord)
        TaylorSeries.add!(r2s2_S, ρ2s2_S, z2s2_S, ord)
        TaylorSeries.sqrt!(r2s_S, r2s2_S, ord)
        TaylorSeries.pow!(r2s5_S, r2s_S, 5, ord)
        TaylorSeries.mul!(tmp1814, X_bf[mo, ea], r_star_M_2[1], ord)
        TaylorSeries.mul!(tmp1815, Y_bf[mo, ea], r_star_M_2[2], ord)
        TaylorSeries.add!(coeff1_2_M, tmp1814, tmp1815, ord)
        TaylorSeries.mul!(tmp1817, X_bf[mo, ea], r_star_S_2[1], ord)
        TaylorSeries.mul!(tmp1818, Y_bf[mo, ea], r_star_S_2[2], ord)
        TaylorSeries.add!(coeff1_2_S, tmp1817, tmp1818, ord)
        TaylorSeries.pow!(tmp1822, coeff1_2_M, 2, ord)
        TaylorSeries.pow!(tmp1825, r_xy[mo, ea], 2, ord)
        TaylorSeries.mul!(tmp1826, 0.5, tmp1825, ord)
        TaylorSeries.mul!(tmp1827, tmp1826, ρ2s2_M, ord)
        TaylorSeries.subst!(tmp1828, tmp1822, tmp1827, ord)
        TaylorSeries.mul!(tmp1829, 5, tmp1828, ord)
        TaylorSeries.div!(coeff3_2_M, tmp1829, r_p2[mo, ea], ord)
        TaylorSeries.pow!(tmp1833, coeff1_2_S, 2, ord)
        TaylorSeries.pow!(tmp1836, r_xy[mo, ea], 2, ord)
        TaylorSeries.mul!(tmp1837, 0.5, tmp1836, ord)
        TaylorSeries.mul!(tmp1838, tmp1837, ρ2s2_S, ord)
        TaylorSeries.subst!(tmp1839, tmp1833, tmp1838, ord)
        TaylorSeries.mul!(tmp1840, 5, tmp1839, ord)
        TaylorSeries.div!(coeff3_2_S, tmp1840, r_p2[mo, ea], ord)
        TaylorSeries.div!(k_22E_div_r2s5_M, k_22E, r2s5_M, ord)
        TaylorSeries.div!(k_22E_div_r2s5_S, k_22E, r2s5_S, ord)
        TaylorSeries.mul!(tmp1845, 2, coeff1_2_M, ord)
        TaylorSeries.mul!(tmp1846, tmp1845, r_star_M_2[1], ord)
        TaylorSeries.add!(tmp1847, ρ2s2_M, coeff3_2_M, ord)
        TaylorSeries.mul!(tmp1848, tmp1847, X_bf[mo, ea], ord)
        TaylorSeries.subst!(tmp1849, tmp1846, tmp1848, ord)
        TaylorSeries.mul!(aux2_M_x, k_22E_div_r2s5_M, tmp1849, ord)
        TaylorSeries.mul!(tmp1852, 2, coeff1_2_M, ord)
        TaylorSeries.mul!(tmp1853, tmp1852, r_star_M_2[2], ord)
        TaylorSeries.add!(tmp1854, ρ2s2_M, coeff3_2_M, ord)
        TaylorSeries.mul!(tmp1855, tmp1854, Y_bf[mo, ea], ord)
        TaylorSeries.subst!(tmp1856, tmp1853, tmp1855, ord)
        TaylorSeries.mul!(aux2_M_y, k_22E_div_r2s5_M, tmp1856, ord)
        TaylorSeries.subst!(tmp1858, coeff3_2_M, ord)
        TaylorSeries.mul!(tmp1859, k_22E_div_r2s5_M, tmp1858, ord)
        TaylorSeries.mul!(aux2_M_z, tmp1859, Z_bf[mo, ea], ord)
        TaylorSeries.mul!(tmp1862, 2, coeff1_2_S, ord)
        TaylorSeries.mul!(tmp1863, tmp1862, r_star_S_2[1], ord)
        TaylorSeries.add!(tmp1864, ρ2s2_S, coeff3_2_S, ord)
        TaylorSeries.mul!(tmp1865, tmp1864, X_bf[mo, ea], ord)
        TaylorSeries.subst!(tmp1866, tmp1863, tmp1865, ord)
        TaylorSeries.mul!(aux2_S_x, k_22E_div_r2s5_S, tmp1866, ord)
        TaylorSeries.mul!(tmp1869, 2, coeff1_2_S, ord)
        TaylorSeries.mul!(tmp1870, tmp1869, r_star_S_2[2], ord)
        TaylorSeries.add!(tmp1871, ρ2s2_S, coeff3_2_S, ord)
        TaylorSeries.mul!(tmp1872, tmp1871, Y_bf[mo, ea], ord)
        TaylorSeries.subst!(tmp1873, tmp1870, tmp1872, ord)
        TaylorSeries.mul!(aux2_S_y, k_22E_div_r2s5_S, tmp1873, ord)
        TaylorSeries.subst!(tmp1875, coeff3_2_S, ord)
        TaylorSeries.mul!(tmp1876, k_22E_div_r2s5_S, tmp1875, ord)
        TaylorSeries.mul!(aux2_S_z, tmp1876, Z_bf[mo, ea], ord)
        TaylorSeries.div!(tmp1878, RE_au, r_p1d2[mo, ea], ord)
        TaylorSeries.pow!(RE_div_r_p5, tmp1878, 5, ord)
        TaylorSeries.mul!(aux_tidacc, tid_num_coeff, RE_div_r_p5, ord)
        TaylorSeries.mul!(tide_acc_coeff_M, μ[mo], aux_tidacc, ord)
        TaylorSeries.mul!(tide_acc_coeff_S, μ[su], aux_tidacc, ord)
        TaylorSeries.add!(tmp1884, aux0_M_x, aux1_M_x, ord)
        TaylorSeries.add!(tmp1885, tmp1884, aux2_M_x, ord)
        TaylorSeries.mul!(tmp1886, tide_acc_coeff_M, tmp1885, ord)
        TaylorSeries.add!(tmp1887, aux0_S_x, aux1_S_x, ord)
        TaylorSeries.add!(tmp1888, tmp1887, aux2_S_x, ord)
        TaylorSeries.mul!(tmp1889, tide_acc_coeff_S, tmp1888, ord)
        TaylorSeries.add!(tidal_bf_x, tmp1886, tmp1889, ord)
        TaylorSeries.add!(tmp1891, aux0_M_y, aux1_M_y, ord)
        TaylorSeries.add!(tmp1892, tmp1891, aux2_M_y, ord)
        TaylorSeries.mul!(tmp1893, tide_acc_coeff_M, tmp1892, ord)
        TaylorSeries.add!(tmp1894, aux0_S_y, aux1_S_y, ord)
        TaylorSeries.add!(tmp1895, tmp1894, aux2_S_y, ord)
        TaylorSeries.mul!(tmp1896, tide_acc_coeff_S, tmp1895, ord)
        TaylorSeries.add!(tidal_bf_y, tmp1893, tmp1896, ord)
        TaylorSeries.add!(tmp1898, aux0_M_z, aux1_M_z, ord)
        TaylorSeries.add!(tmp1899, tmp1898, aux2_M_z, ord)
        TaylorSeries.mul!(tmp1900, tide_acc_coeff_M, tmp1899, ord)
        TaylorSeries.add!(tmp1901, aux0_S_z, aux1_S_z, ord)
        TaylorSeries.add!(tmp1902, tmp1901, aux2_S_z, ord)
        TaylorSeries.mul!(tmp1903, tide_acc_coeff_S, tmp1902, ord)
        TaylorSeries.add!(tidal_bf_z, tmp1900, tmp1903, ord)
        TaylorSeries.mul!(tmp1905, M_[1, 1, ea], tidal_bf_x, ord)
        TaylorSeries.mul!(tmp1906, M_[2, 1, ea], tidal_bf_y, ord)
        TaylorSeries.add!(tmp1907, tmp1905, tmp1906, ord)
        TaylorSeries.mul!(tmp1908, M_[3, 1, ea], tidal_bf_z, ord)
        TaylorSeries.add!(tidal_x, tmp1907, tmp1908, ord)
        TaylorSeries.mul!(tmp1910, M_[1, 2, ea], tidal_bf_x, ord)
        TaylorSeries.mul!(tmp1911, M_[2, 2, ea], tidal_bf_y, ord)
        TaylorSeries.add!(tmp1912, tmp1910, tmp1911, ord)
        TaylorSeries.mul!(tmp1913, M_[3, 2, ea], tidal_bf_z, ord)
        TaylorSeries.add!(tidal_y, tmp1912, tmp1913, ord)
        TaylorSeries.mul!(tmp1915, M_[1, 3, ea], tidal_bf_x, ord)
        TaylorSeries.mul!(tmp1916, M_[2, 3, ea], tidal_bf_y, ord)
        TaylorSeries.add!(tmp1917, tmp1915, tmp1916, ord)
        TaylorSeries.mul!(tmp1918, M_[3, 3, ea], tidal_bf_z, ord)
        TaylorSeries.add!(tidal_z, tmp1917, tmp1918, ord)
        TaylorSeries.add!(accX_mo_tides, accX[mo], tidal_x, ord)
        TaylorSeries.add!(accY_mo_tides, accY[mo], tidal_y, ord)
        TaylorSeries.add!(accZ_mo_tides, accZ[mo], tidal_z, ord)
        TaylorSeries.identity!(accX[mo], accX_mo_tides, ord)
        TaylorSeries.identity!(accY[mo], accY_mo_tides, ord)
        TaylorSeries.identity!(accZ[mo], accZ_mo_tides, ord)
        #= REPL[2]:744 =# Threads.@threads for i = 1:N_ext
                TaylorSeries.add!(dq[3 * (N + i) - 2], postNewtonX[i, postnewton_iter + 1], accX[i], ord)
                TaylorSeries.add!(dq[3 * (N + i) - 1], postNewtonY[i, postnewton_iter + 1], accY[i], ord)
                TaylorSeries.add!(dq[3 * (N + i)], postNewtonZ[i, postnewton_iter + 1], accZ[i], ord)
            end
        #= REPL[2]:749 =# Threads.@threads for i = N_ext + 1:N
                TaylorSeries.identity!(dq[3 * (N + i) - 2], postNewtonX[i, postnewton_iter + 1], ord)
                TaylorSeries.identity!(dq[3 * (N + i) - 1], postNewtonY[i, postnewton_iter + 1], ord)
                TaylorSeries.identity!(dq[3 * (N + i)], postNewtonZ[i, postnewton_iter + 1], ord)
            end
        for __idx = eachindex(q)
            (q[__idx]).coeffs[ordnext + 1] = (dq[__idx]).coeffs[ordnext] / ordnext
        end
    end
    return nothing
end