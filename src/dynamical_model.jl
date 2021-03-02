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
@taylorize function NBP_pN_A_J23E_J23M_J2S!(dq, q, params, t)
    # N: number of bodies
    # jd0: initial Julian date
    local N, jd0 = params
    local S = eltype(q)
    local N_ext = 11 # number of bodies in extended-body accelerations

    local zero_q_1 = zero(q[1])
    local one_t = one(t)
    local dsj2k = t+(jd0-J2000) # days since J2000.0 (TDB)
    local ϕ0, θ0, ψ0, ωx0, ωy0, ωz0 = q[6N+1:6N+6]
    local dϕ0 = ((ωx0*sin(ψ0)) + (ωy0*cos(ψ0)) )/sin(θ0)
    local dθ0 = (ωx0*cos(ψ0)) - (ωy0*sin(ψ0))
    local dψ0 = ωz0 - (dϕ0*cos(θ0))
    local __t = Taylor1(t.order)
    local eulang_t = [ϕ0 + __t*dϕ0, θ0 + __t*dθ0, ψ0 + __t*dψ0]
    local I_m_t = ITM_und.*one_t # ITM(q_del_τ_M, eulang_t_del) # matrix elements of lunar moment of inertia (Folkner et al. 2014, eq. 41)
    local dI_m_t = derivative.(I_m_t) # time-derivative of lunar mantle I at time t
    local inv_I_m_t = inv(I_m_t) # inverse of lunar mantle I matrix at time t

    # parameters related to speed of light, c
    local c_p2 = 29979.063823897606 # c^2 = 29979.063823897606 au^2/d^2
    local c_m2 = 3.3356611996764786e-5 # c^-2 = 3.3356611996764786e-5 d^2/au^2

    X = Array{S}(undef, N, N)
    Y = Array{S}(undef, N, N)
    Z = Array{S}(undef, N, N)

    r_p2 = Array{S}(undef, N, N)
    r_p3d2 = Array{S}(undef, N, N)
    r_p7d2 = Array{S}(undef, N, N)

    #Newtonian accelerations
    newtonX = Array{S}(undef, N)
    newtonY = Array{S}(undef, N)
    newtonZ = Array{S}(undef, N)

    newtonianCoeff = Array{S}(undef, N, N)

    #post-Newtonian stuff
    U = Array{S}(undef, N, N)
    V = Array{S}(undef, N, N)
    W = Array{S}(undef, N, N)

    _4U_m_3X = Array{S}(undef, N, N)
    _4V_m_3Y = Array{S}(undef, N, N)
    _4W_m_3Z = Array{S}(undef, N, N)

    UU = Array{S}(undef, N, N)
    VV = Array{S}(undef, N, N)
    WW = Array{S}(undef, N, N)

    r_p1d2 = Array{S}(undef, N, N)

    newtonianNb_Potential = Array{S}(undef, N)
    newtonian1b_Potential = Array{S}(undef, N, N)
    newtonianCoeff = Array{S}(undef, N, N)
    newton_acc_X = Array{S}(undef, N, N)
    newton_acc_Y = Array{S}(undef, N, N)
    newton_acc_Z = Array{S}(undef, N, N)

    v2 = Array{S}(undef, N)
    vi_dot_vj = Array{S}(undef, N, N)
    pn2 = Array{S}(undef, N, N)
    pn3 = Array{S}(undef, N, N)
    _4ϕj = Array{S}(undef, N, N)
    ϕi_plus_4ϕj = Array{S}(undef, N, N)
    sj2_plus_2si2_minus_4vivj = Array{S}(undef, N, N)
    ϕs_and_vs = Array{S}(undef, N, N)
    U_t_pn2 = Array{S}(undef, N, N)
    V_t_pn2 = Array{S}(undef, N, N)
    W_t_pn2 = Array{S}(undef, N, N)
    pn1t1_7 = Array{S}(undef, N, N)

    pntempX = Array{S}(undef, N)
    pntempY = Array{S}(undef, N)
    pntempZ = Array{S}(undef, N)
    pn1 = Array{S}(undef, N, N)
    X_t_pn1 = Array{S}(undef, N, N)
    Y_t_pn1 = Array{S}(undef, N, N)
    Z_t_pn1 = Array{S}(undef, N, N)
    pNX_t_pn3 = Array{S}(undef, N, N)
    pNY_t_pn3 = Array{S}(undef, N, N)
    pNZ_t_pn3 = Array{S}(undef, N, N)
    pNX_t_X = Array{S}(undef, N, N)
    pNY_t_Y = Array{S}(undef, N, N)
    pNZ_t_Z = Array{S}(undef, N, N)
    postNewtonX = Array{S}(undef, N)
    postNewtonY = Array{S}(undef, N)
    postNewtonZ = Array{S}(undef, N)

    # (Jn, Cmn, Smn) acceleration auxiliaries
    X_bf_1 = Array{S}(undef, N_ext, N_ext)
    Y_bf_1 = Array{S}(undef, N_ext, N_ext)
    Z_bf_1 = Array{S}(undef, N_ext, N_ext)
    X_bf_2 = Array{S}(undef, N_ext, N_ext)
    Y_bf_2 = Array{S}(undef, N_ext, N_ext)
    Z_bf_2 = Array{S}(undef, N_ext, N_ext)
    X_bf_3 = Array{S}(undef, N_ext, N_ext)
    Y_bf_3 = Array{S}(undef, N_ext, N_ext)
    Z_bf_3 = Array{S}(undef, N_ext, N_ext)
    X_bf = Array{S}(undef, N_ext, N_ext)
    Y_bf = Array{S}(undef, N_ext, N_ext)
    Z_bf = Array{S}(undef, N_ext, N_ext)
    F_JCS_x = Array{S}(undef, N_ext, N_ext)
    F_JCS_y = Array{S}(undef, N_ext, N_ext)
    F_JCS_z = Array{S}(undef, N_ext, N_ext)
    temp_accX_j = Array{S}(undef, N_ext, N_ext)
    temp_accY_j = Array{S}(undef, N_ext, N_ext)
    temp_accZ_j = Array{S}(undef, N_ext, N_ext)
    temp_accX_i = Array{S}(undef, N_ext, N_ext)
    temp_accY_i = Array{S}(undef, N_ext, N_ext)
    temp_accZ_i = Array{S}(undef, N_ext, N_ext)
    sin_ϕ = Array{S}(undef, N_ext, N_ext)
    cos_ϕ = Array{S}(undef, N_ext, N_ext)
    sin_λ = Array{S}(undef, N_ext, N_ext)
    cos_λ = Array{S}(undef, N_ext, N_ext)
    r_xy = Array{S}(undef, N_ext, N_ext)
    r_p4 = Array{S}(undef, N_ext, N_ext)
    P_n = Array{S}(undef, N_ext, N_ext, maximum(n1SEM)+1)
    dP_n = Array{S}(undef, N_ext, N_ext, maximum(n1SEM)+1)
    temp_fjξ = Array{S}(undef, N_ext, N_ext, maximum(n1SEM)+1)
    temp_fjζ = Array{S}(undef, N_ext, N_ext, maximum(n1SEM)+1)
    temp_rn = Array{S}(undef, N_ext, N_ext, maximum(n1SEM)+1)
    F_CS_ξ_36 = Array{S}(undef, N_ext, N_ext)
    F_CS_η_36 = Array{S}(undef, N_ext, N_ext)
    F_CS_ζ_36 = Array{S}(undef, N_ext, N_ext)
    F_J_ξ_36 = Array{S}(undef, N_ext, N_ext)
    F_J_ζ_36 = Array{S}(undef, N_ext, N_ext)
    sin_mλ = Array{S}(undef, N_ext, N_ext, n1SEM[mo])
    cos_mλ = Array{S}(undef, N_ext, N_ext, n1SEM[mo])
    secϕ_P_nm = Array{S}(undef, N_ext, N_ext, n1SEM[mo]+1, n1SEM[mo]+1)
    P_nm = Array{S}(undef, N_ext, N_ext, n1SEM[mo]+1, n1SEM[mo]+1)
    cosϕ_dP_nm = Array{S}(undef, N_ext, N_ext, n1SEM[mo]+1, n1SEM[mo]+1)
    F_J_ξ = Array{S}(undef, N_ext, N_ext)
    F_J_η = Array{S}(undef, N_ext, N_ext)
    F_J_ζ = Array{S}(undef, N_ext, N_ext)
    F_CS_ξ = Array{S}(undef, N_ext, N_ext)
    F_CS_η = Array{S}(undef, N_ext, N_ext)
    F_CS_ζ = Array{S}(undef, N_ext, N_ext)
    F_JCS_ξ = Array{S}(undef, N_ext, N_ext)
    F_JCS_η = Array{S}(undef, N_ext, N_ext)
    F_JCS_ζ = Array{S}(undef, N_ext, N_ext)
    Rb2p = Array{S}(undef, N_ext, N_ext, 3, 3) #R matrix body-fixed to "primed" ξηζ frame (Moyer, 1971, eq. 161)
    Gc2p = Array{S}(undef, N_ext, N_ext, 3, 3) #G matrix "space-fixed" to "primed" ξηζ frame (Moyer, 1971, eq. 163)

    # extended-body accelerations
    accX = Array{S}(undef, N_ext)
    accY = Array{S}(undef, N_ext)
    accZ = Array{S}(undef, N_ext)

    # lunar torques: Folkner et al. (2014), Eq. (43)
    N_MfigM_pmA_x = Array{S}(undef, N_ext)
    N_MfigM_pmA_y = Array{S}(undef, N_ext)
    N_MfigM_pmA_z = Array{S}(undef, N_ext)
    temp_N_M_x = Array{S}(undef, N_ext)
    temp_N_M_y = Array{S}(undef, N_ext)
    temp_N_M_z = Array{S}(undef, N_ext)
    N_MfigM = Array{S}(undef, 3)
    N_MfigM[1] = zero_q_1
    N_MfigM[2] = zero_q_1
    N_MfigM[3] = zero_q_1

    # rotations to and from Earth, Sun and Moon pole-oriented frames
    local αs = deg2rad(α_p_sun*one_t)
    local δs = deg2rad(δ_p_sun*one_t)
    local αm = eulang_t[1] - (pi/2)
    local δm = (pi/2) - eulang_t[2]
    local Wm = eulang_t[3]
    RotM = Array{S}(undef, 3, 3, 5) # space-fixed -> body-fixed coordinate transformations
    local RotM[:,:,ea] = c2t_jpl_de430(dsj2k)
    local RotM[:,:,su] = pole_rotation(αs, δs)
    local RotM[:,:,mo] = pole_rotation(αm, δm, Wm)
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
    J2_t = Array{S}(undef, 5)
    J2_t[su] = J2S_t
    J2_t[ea] = J2E_t
    # lunar torques: overall numerical factor in Eq. (44) of Folkner et al. (2014)
    local N_MfigM_figE_factor = 7.5*μ[ea]*J2E_t

    for j in 1:N
        newtonX[j] = zero_q_1
        newtonY[j] = zero_q_1
        newtonZ[j] = zero_q_1
        newtonianNb_Potential[j] = zero_q_1
        dq[3j-2] = q[3(N+j)-2]
        dq[3j-1] = q[3(N+j)-1]
        dq[3j  ] = q[3(N+j)  ]
    end

    for j in 1:N_ext
        accX[j] = zero_q_1
        accY[j] = zero_q_1
        accZ[j] = zero_q_1
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

    # 2nd-order lunar zonal and tesseral harmonics
    J2M_t = ( I_m_t[3,3] - ((I_m_t[1,1]+I_m_t[2,2])/2) )/(μ[mo]) # J_{2,M}*R_M^2
    C22M_t = ((I_m_t[2,2] - I_m_t[1,1])/(μ[mo]))/4 # C_{22,M}*R_M^2
    C21M_t = (-I_m_t[1,3])/(μ[mo]) # C_{21,M}*R_M^2
    S21M_t = (-I_m_t[3,2])/(μ[mo]) # S_{21,M}*R_M^2
    S22M_t = ((-I_m_t[2,1])/(μ[mo]))/2 # S_{22,M}*R_M^2
    J2_t[mo] = J2M_t

    for j in 1:N_ext
        for i in 1:N_ext
            # i == j && continue
            if i == j
                continue
            else
                #Jn, Cnm, Snm accelerations, if j-th body is flattened
                if UJ_interaction[i,j]
                    # rotate from inertial frame to extended-body frame
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
                                P_nm[i,j,1,1] = cos_ϕ[i,j]
                                cosϕ_dP_nm[i,j,1,1] = sin_ϕ[i,j]*lnm3[1] #+0 (second term in Eq. 183 from Moyer, 1971, vanishes when n=m)
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
                    # G_{i,j} = \sum_k R_{i,k} RotM{k,j}
                    Gc2p[i,j,1,1] = ((Rb2p[i,j,1,1]*RotM[1,1,j]) + (Rb2p[i,j,1,2]*RotM[2,1,j])) + (Rb2p[i,j,1,3]*RotM[3,1,j])
                    Gc2p[i,j,2,1] = ((Rb2p[i,j,2,1]*RotM[1,1,j]) + (Rb2p[i,j,2,2]*RotM[2,1,j])) + (Rb2p[i,j,2,3]*RotM[3,1,j])
                    Gc2p[i,j,3,1] = ((Rb2p[i,j,3,1]*RotM[1,1,j]) + (Rb2p[i,j,3,2]*RotM[2,1,j])) + (Rb2p[i,j,3,3]*RotM[3,1,j])
                    Gc2p[i,j,1,2] = ((Rb2p[i,j,1,1]*RotM[1,2,j]) + (Rb2p[i,j,1,2]*RotM[2,2,j])) + (Rb2p[i,j,1,3]*RotM[3,2,j])
                    Gc2p[i,j,2,2] = ((Rb2p[i,j,2,1]*RotM[1,2,j]) + (Rb2p[i,j,2,2]*RotM[2,2,j])) + (Rb2p[i,j,2,3]*RotM[3,2,j])
                    Gc2p[i,j,3,2] = ((Rb2p[i,j,3,1]*RotM[1,2,j]) + (Rb2p[i,j,3,2]*RotM[2,2,j])) + (Rb2p[i,j,3,3]*RotM[3,2,j])
                    Gc2p[i,j,1,3] = ((Rb2p[i,j,1,1]*RotM[1,3,j]) + (Rb2p[i,j,1,2]*RotM[2,3,j])) + (Rb2p[i,j,1,3]*RotM[3,3,j])
                    Gc2p[i,j,2,3] = ((Rb2p[i,j,2,1]*RotM[1,3,j]) + (Rb2p[i,j,2,2]*RotM[2,3,j])) + (Rb2p[i,j,2,3]*RotM[3,3,j])
                    Gc2p[i,j,3,3] = ((Rb2p[i,j,3,1]*RotM[1,3,j]) + (Rb2p[i,j,3,2]*RotM[2,3,j])) + (Rb2p[i,j,3,3]*RotM[3,3,j])
                    # compute cartesian coordinates of acceleration due to body figure in inertial frame (Moyer, 1971, eq. (169))
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

                    if j == mo
                        # compute torques acting upon the body-figure of the Moon due to external point masses
                        # Folkner et al. (2014), Eq. (43)
                        N_MfigM_pmA_x[i] = μ[i]*( (Y[i,j]*F_JCS_z[i,j]) - (Z[i,j]*F_JCS_y[i,j]) )
                        N_MfigM_pmA_y[i] = μ[i]*( (Z[i,j]*F_JCS_x[i,j]) - (X[i,j]*F_JCS_z[i,j]) )
                        N_MfigM_pmA_z[i] = μ[i]*( (X[i,j]*F_JCS_y[i,j]) - (Y[i,j]*F_JCS_x[i,j]) )
                        # expressions below have minus sign since N_MfigM_pmA_{x,y,z} have inverted signs in cross product
                        temp_N_M_x[i] = N_MfigM[1] - (μ[j]*N_MfigM_pmA_x[i])
                        N_MfigM[1] = temp_N_M_x[i]
                        temp_N_M_y[i] = N_MfigM[2] - (μ[j]*N_MfigM_pmA_y[i])
                        N_MfigM[2] = temp_N_M_y[i]
                        temp_N_M_z[i] = N_MfigM[3] - (μ[j]*N_MfigM_pmA_z[i])
                        N_MfigM[3] = temp_N_M_z[i]
                    end
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
                ###
                pn1[i,j] = zero_q_1
                X_t_pn1[i,j] = zero_q_1
                Y_t_pn1[i,j] = zero_q_1
                Z_t_pn1[i,j] = zero_q_1
                pNX_t_pn3[i,j] = zero_q_1
                pNY_t_pn3[i,j] = zero_q_1
                pNZ_t_pn3[i,j] = zero_q_1
                pNX_t_X[i,j] = zero_q_1
                pNY_t_Y[i,j] = zero_q_1
                pNZ_t_Z[i,j] = zero_q_1
            end # else (i != j)
        end
        pntempX[j] = zero_q_1
        pntempY[j] = zero_q_1
        pntempZ[j] = zero_q_1
    end

    for j in 1:N
        for i in 1:N
            # i == j && continue
            if i == j
                continue
            else
                pNX_t_X[i,j] = newtonX[i]*X[i,j]
                pNY_t_Y[i,j] = newtonY[i]*Y[i,j]
                pNZ_t_Z[i,j] = newtonZ[i]*Z[i,j]
                pn1[i,j] = (  pn1t1_7[i,j]  +  0.5*( (pNX_t_X[i,j]+pNY_t_Y[i,j]) + pNZ_t_Z[i,j] )  )

                X_t_pn1[i,j] = newton_acc_X[i,j]*pn1[i,j]
                Y_t_pn1[i,j] = newton_acc_Y[i,j]*pn1[i,j]
                Z_t_pn1[i,j] = newton_acc_Z[i,j]*pn1[i,j]

                pNX_t_pn3[i,j] = newtonX[i]*pn3[i,j]
                pNY_t_pn3[i,j] = newtonY[i]*pn3[i,j]
                pNZ_t_pn3[i,j] = newtonZ[i]*pn3[i,j]

                termpnx = ( X_t_pn1[i,j] + (U_t_pn2[i,j]+pNX_t_pn3[i,j]) )
                sumpnx = pntempX[j] + termpnx
                pntempX[j] = sumpnx
                termpny = ( Y_t_pn1[i,j] + (V_t_pn2[i,j]+pNY_t_pn3[i,j]) )
                sumpny = pntempY[j] + termpny
                pntempY[j] = sumpny
                termpnz = ( Z_t_pn1[i,j] + (W_t_pn2[i,j]+pNZ_t_pn3[i,j]) )
                sumpnz = pntempZ[j] + termpnz
                pntempZ[j] = sumpnz
            end # else (i != j)
        end
        postNewtonX[j] = pntempX[j]*c_m2
        postNewtonY[j] = pntempY[j]*c_m2
        postNewtonZ[j] = pntempZ[j]*c_m2
    end

    #fill accelerations (post-Newtonian and extended body accelerations)
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

    Iω_x = (I_m_t[1,1]*q[6N+4]) + ((I_m_t[1,2]*q[6N+5]) + (I_m_t[1,3]*q[6N+6]))
    Iω_y = (I_m_t[2,1]*q[6N+4]) + ((I_m_t[2,2]*q[6N+5]) + (I_m_t[2,3]*q[6N+6]))
    Iω_z = (I_m_t[3,1]*q[6N+4]) + ((I_m_t[3,2]*q[6N+5]) + (I_m_t[3,3]*q[6N+6]))

    # ω × (I*ω)
    ωxIω_x = (q[6N+5]*Iω_z) - (q[6N+6]*Iω_y)
    ωxIω_y = (q[6N+6]*Iω_x) - (q[6N+4]*Iω_z)
    ωxIω_z = (q[6N+4]*Iω_y) - (q[6N+5]*Iω_x)

    # (dI/dt)*ω
    dIω_x = (dI_m_t[1,1]*q[6N+4]) + ((dI_m_t[1,2]*q[6N+5]) + (dI_m_t[1,3]*q[6N+6]))
    dIω_y = (dI_m_t[2,1]*q[6N+4]) + ((dI_m_t[2,2]*q[6N+5]) + (dI_m_t[2,3]*q[6N+6]))
    dIω_z = (dI_m_t[3,1]*q[6N+4]) + ((dI_m_t[3,2]*q[6N+5]) + (dI_m_t[3,3]*q[6N+6]))

    # Moon->Earth radial unit vector (inertial coordinates)
    er_EM_I_1 = X[ea,mo]/r_p1d2[ea,mo]
    er_EM_I_2 = Y[ea,mo]/r_p1d2[ea,mo]
    er_EM_I_3 = Z[ea,mo]/r_p1d2[ea,mo]

    # Earth pole unit vector (inertial coordinates)
    p_E_I_1 = RotM[3,1,ea]
    p_E_I_2 = RotM[3,2,ea]
    p_E_I_3 = RotM[3,3,ea]

    # transform er_EM_I_i, p_E_I_i to lunar mantle frame coordinates
    er_EM_1 = (RotM[1,1,mo]*er_EM_I_1) + ((RotM[1,2,mo]*er_EM_I_2) + (RotM[1,3,mo]*er_EM_I_3))
    er_EM_2 = (RotM[2,1,mo]*er_EM_I_1) + ((RotM[2,2,mo]*er_EM_I_2) + (RotM[2,3,mo]*er_EM_I_3))
    er_EM_3 = (RotM[3,1,mo]*er_EM_I_1) + ((RotM[3,2,mo]*er_EM_I_2) + (RotM[3,3,mo]*er_EM_I_3))
    p_E_1 = (RotM[1,1,mo]*p_E_I_1) + ((RotM[1,2,mo]*p_E_I_2) + (RotM[1,3,mo]*p_E_I_3))
    p_E_2 = (RotM[2,1,mo]*p_E_I_1) + ((RotM[2,2,mo]*p_E_I_2) + (RotM[2,3,mo]*p_E_I_3))
    p_E_3 = (RotM[3,1,mo]*p_E_I_1) + ((RotM[3,2,mo]*p_E_I_2) + (RotM[3,3,mo]*p_E_I_3))

    # evaluate Eq. (44) of Folkner et al. (2014) in lunar mantle frame coords
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

    one_minus_7sin2ϕEM = one_t - (7((sin_ϕ[mo,ea])^2))
    two_sinϕEM = 2sin_ϕ[mo,ea]

    N_MfigM_figE_factor_div_rEMp5 = (N_MfigM_figE_factor/(r_p1d2[mo,ea]^5))
    N_MfigM_figE_1 = N_MfigM_figE_factor_div_rEMp5*( (one_minus_7sin2ϕEM*er_EM_cross_I_er_EM_1) + (two_sinϕEM*(er_EM_cross_I_p_E_1+p_E_cross_I_er_EM_1)) - (0.4p_E_cross_I_p_E_1))
    N_MfigM_figE_2 = N_MfigM_figE_factor_div_rEMp5*( (one_minus_7sin2ϕEM*er_EM_cross_I_er_EM_2) + (two_sinϕEM*(er_EM_cross_I_p_E_2+p_E_cross_I_er_EM_2)) - (0.4p_E_cross_I_p_E_2))
    N_MfigM_figE_3 = N_MfigM_figE_factor_div_rEMp5*( (one_minus_7sin2ϕEM*er_EM_cross_I_er_EM_3) + (two_sinϕEM*(er_EM_cross_I_p_E_3+p_E_cross_I_er_EM_3)) - (0.4p_E_cross_I_p_E_3))

    # torques acting upon lunar body-figure due to external point masses: transform coordinates from inertial frame to lunar mantle frame
    N_1_LMF = (RotM[1,1,mo]*N_MfigM[1]) + ((RotM[1,2,mo]*N_MfigM[2]) + (RotM[1,3,mo]*N_MfigM[3]))
    N_2_LMF = (RotM[2,1,mo]*N_MfigM[1]) + ((RotM[2,2,mo]*N_MfigM[2]) + (RotM[2,3,mo]*N_MfigM[3]))
    N_3_LMF = (RotM[3,1,mo]*N_MfigM[1]) + ((RotM[3,2,mo]*N_MfigM[2]) + (RotM[3,3,mo]*N_MfigM[3]))

    # I*(dω/dt)
    I_dω_1 = (N_1_LMF + N_MfigM_figE_1) - (dIω_x + ωxIω_x)
    I_dω_2 = (N_2_LMF + N_MfigM_figE_2) - (dIω_y + ωxIω_y)
    I_dω_3 = (N_3_LMF + N_MfigM_figE_3) - (dIω_z + ωxIω_z)

    # lunar physical librations: Folkner et al. (2014), Eq. (14)
    dq[6N+1] = ((q[6N+4]*sin(q[6N+3])) + (q[6N+5]*cos(q[6N+3])) )/sin(q[6N+2])
    dq[6N+2] = (q[6N+4]*cos(q[6N+3])) - (q[6N+5]*sin(q[6N+3]))
    dq[6N+3] = q[6N+6] - (dq[6N+1]*cos(q[6N+2]))

    # lunar physical librations: Folkner et al. (2014), Eq. (34)
    # TODO: add lunar mantle-core boundary interaction
    dq[6N+4] = (inv_I_m_t[1,1]*I_dω_1) + ( (inv_I_m_t[1,2]*I_dω_2) + (inv_I_m_t[1,3]*I_dω_3) )
    dq[6N+5] = (inv_I_m_t[2,1]*I_dω_1) + ( (inv_I_m_t[2,2]*I_dω_2) + (inv_I_m_t[2,3]*I_dω_3) )
    dq[6N+6] = (inv_I_m_t[3,1]*I_dω_1) + ( (inv_I_m_t[3,2]*I_dω_2) + (inv_I_m_t[3,3]*I_dω_3) )

    # TT-TDB
    # TODO: implement TT-TDB integration
    dq[6N+7] = zero_q_1

    nothing
end

@taylorize function NBP_pN_A_J23E_J23M_J2S_threads!(dq, q, params, t)
    # N: number of bodies
    # jd0: initial Julian date
    local N, jd0 = params
    local S = eltype(q)
    local N_ext = 11 # number of bodies in extended-body accelerations

    local zero_q_1 = zero(q[1])
    local one_t = one(t)
    local dsj2k = t+(jd0-J2000) # days since J2000.0 (TDB)
    local ϕ0, θ0, ψ0, ωx0, ωy0, ωz0 = q[6N+1:6N+6]
    local dϕ0 = ((ωx0*sin(ψ0)) + (ωy0*cos(ψ0)) )/sin(θ0)
    local dθ0 = (ωx0*cos(ψ0)) - (ωy0*sin(ψ0))
    local dψ0 = ωz0 - (dϕ0*cos(θ0))
    local __t = Taylor1(t.order)
    local eulang_t = [ϕ0 + __t*dϕ0, θ0 + __t*dθ0, ψ0 + __t*dψ0]
    local I_m_t = ITM_und.*one_t # ITM(q_del_τ_M, eulang_t_del) # matrix elements of lunar moment of inertia (Folkner et al. 2014, eq. 41)
    local dI_m_t = derivative.(I_m_t) # time-derivative of lunar mantle I at time t
    local inv_I_m_t = inv(I_m_t) # inverse of lunar mantle I matrix at time t

    # parameters related to speed of light, c
    local c_p2 = 29979.063823897606 # c^2 = 29979.063823897606 au^2/d^2
    local c_m2 = 3.3356611996764786e-5 # c^-2 = 3.3356611996764786e-5 d^2/au^2

    X = Array{S}(undef, N, N)
    Y = Array{S}(undef, N, N)
    Z = Array{S}(undef, N, N)

    r_p2 = Array{S}(undef, N, N)
    r_p3d2 = Array{S}(undef, N, N)
    r_p7d2 = Array{S}(undef, N, N)

    #Newtonian accelerations
    newtonX = Array{S}(undef, N)
    newtonY = Array{S}(undef, N)
    newtonZ = Array{S}(undef, N)

    newtonianCoeff = Array{S}(undef, N, N)

    #post-Newtonian stuff
    U = Array{S}(undef, N, N)
    V = Array{S}(undef, N, N)
    W = Array{S}(undef, N, N)

    _4U_m_3X = Array{S}(undef, N, N)
    _4V_m_3Y = Array{S}(undef, N, N)
    _4W_m_3Z = Array{S}(undef, N, N)

    UU = Array{S}(undef, N, N)
    VV = Array{S}(undef, N, N)
    WW = Array{S}(undef, N, N)

    r_p1d2 = Array{S}(undef, N, N)

    newtonianNb_Potential = Array{S}(undef, N)
    newtonian1b_Potential = Array{S}(undef, N, N)
    newtonianCoeff = Array{S}(undef, N, N)
    newton_acc_X = Array{S}(undef, N, N)
    newton_acc_Y = Array{S}(undef, N, N)
    newton_acc_Z = Array{S}(undef, N, N)

    v2 = Array{S}(undef, N)
    vi_dot_vj = Array{S}(undef, N, N)
    pn2 = Array{S}(undef, N, N)
    pn3 = Array{S}(undef, N, N)
    _4ϕj = Array{S}(undef, N, N)
    ϕi_plus_4ϕj = Array{S}(undef, N, N)
    sj2_plus_2si2_minus_4vivj = Array{S}(undef, N, N)
    ϕs_and_vs = Array{S}(undef, N, N)
    U_t_pn2 = Array{S}(undef, N, N)
    V_t_pn2 = Array{S}(undef, N, N)
    W_t_pn2 = Array{S}(undef, N, N)
    pn1t1_7 = Array{S}(undef, N, N)

    pntempX = Array{S}(undef, N)
    pntempY = Array{S}(undef, N)
    pntempZ = Array{S}(undef, N)
    pn1 = Array{S}(undef, N, N)
    X_t_pn1 = Array{S}(undef, N, N)
    Y_t_pn1 = Array{S}(undef, N, N)
    Z_t_pn1 = Array{S}(undef, N, N)
    pNX_t_pn3 = Array{S}(undef, N, N)
    pNY_t_pn3 = Array{S}(undef, N, N)
    pNZ_t_pn3 = Array{S}(undef, N, N)
    pNX_t_X = Array{S}(undef, N, N)
    pNY_t_Y = Array{S}(undef, N, N)
    pNZ_t_Z = Array{S}(undef, N, N)
    postNewtonX = Array{S}(undef, N)
    postNewtonY = Array{S}(undef, N)
    postNewtonZ = Array{S}(undef, N)

    # (Jn, Cmn, Smn) acceleration auxiliaries
    X_bf_1 = Array{S}(undef, N_ext, N_ext)
    Y_bf_1 = Array{S}(undef, N_ext, N_ext)
    Z_bf_1 = Array{S}(undef, N_ext, N_ext)
    X_bf_2 = Array{S}(undef, N_ext, N_ext)
    Y_bf_2 = Array{S}(undef, N_ext, N_ext)
    Z_bf_2 = Array{S}(undef, N_ext, N_ext)
    X_bf_3 = Array{S}(undef, N_ext, N_ext)
    Y_bf_3 = Array{S}(undef, N_ext, N_ext)
    Z_bf_3 = Array{S}(undef, N_ext, N_ext)
    X_bf = Array{S}(undef, N_ext, N_ext)
    Y_bf = Array{S}(undef, N_ext, N_ext)
    Z_bf = Array{S}(undef, N_ext, N_ext)
    F_JCS_x = Array{S}(undef, N_ext, N_ext)
    F_JCS_y = Array{S}(undef, N_ext, N_ext)
    F_JCS_z = Array{S}(undef, N_ext, N_ext)
    temp_accX_j = Array{S}(undef, N_ext, N_ext)
    temp_accY_j = Array{S}(undef, N_ext, N_ext)
    temp_accZ_j = Array{S}(undef, N_ext, N_ext)
    temp_accX_i = Array{S}(undef, N_ext, N_ext)
    temp_accY_i = Array{S}(undef, N_ext, N_ext)
    temp_accZ_i = Array{S}(undef, N_ext, N_ext)
    sin_ϕ = Array{S}(undef, N_ext, N_ext)
    cos_ϕ = Array{S}(undef, N_ext, N_ext)
    sin_λ = Array{S}(undef, N_ext, N_ext)
    cos_λ = Array{S}(undef, N_ext, N_ext)
    r_xy = Array{S}(undef, N_ext, N_ext)
    r_p4 = Array{S}(undef, N_ext, N_ext)
    P_n = Array{S}(undef, N_ext, N_ext, maximum(n1SEM)+1)
    dP_n = Array{S}(undef, N_ext, N_ext, maximum(n1SEM)+1)
    temp_fjξ = Array{S}(undef, N_ext, N_ext, maximum(n1SEM)+1)
    temp_fjζ = Array{S}(undef, N_ext, N_ext, maximum(n1SEM)+1)
    temp_rn = Array{S}(undef, N_ext, N_ext, maximum(n1SEM)+1)
    F_CS_ξ_36 = Array{S}(undef, N_ext, N_ext)
    F_CS_η_36 = Array{S}(undef, N_ext, N_ext)
    F_CS_ζ_36 = Array{S}(undef, N_ext, N_ext)
    F_J_ξ_36 = Array{S}(undef, N_ext, N_ext)
    F_J_ζ_36 = Array{S}(undef, N_ext, N_ext)
    sin_mλ = Array{S}(undef, N_ext, N_ext, n1SEM[mo])
    cos_mλ = Array{S}(undef, N_ext, N_ext, n1SEM[mo])
    secϕ_P_nm = Array{S}(undef, N_ext, N_ext, n1SEM[mo]+1, n1SEM[mo]+1)
    P_nm = Array{S}(undef, N_ext, N_ext, n1SEM[mo]+1, n1SEM[mo]+1)
    cosϕ_dP_nm = Array{S}(undef, N_ext, N_ext, n1SEM[mo]+1, n1SEM[mo]+1)
    F_J_ξ = Array{S}(undef, N_ext, N_ext)
    F_J_η = Array{S}(undef, N_ext, N_ext)
    F_J_ζ = Array{S}(undef, N_ext, N_ext)
    F_CS_ξ = Array{S}(undef, N_ext, N_ext)
    F_CS_η = Array{S}(undef, N_ext, N_ext)
    F_CS_ζ = Array{S}(undef, N_ext, N_ext)
    F_JCS_ξ = Array{S}(undef, N_ext, N_ext)
    F_JCS_η = Array{S}(undef, N_ext, N_ext)
    F_JCS_ζ = Array{S}(undef, N_ext, N_ext)
    Rb2p = Array{S}(undef, N_ext, N_ext, 3, 3) #R matrix body-fixed to "primed" ξηζ frame (Moyer, 1971, eq. 161)
    Gc2p = Array{S}(undef, N_ext, N_ext, 3, 3) #G matrix "space-fixed" to "primed" ξηζ frame (Moyer, 1971, eq. 163)

    # extended-body accelerations
    accX = Array{S}(undef, N_ext)
    accY = Array{S}(undef, N_ext)
    accZ = Array{S}(undef, N_ext)

    # lunar torques: Folkner et al. (2014), Eq. (43)
    N_MfigM_pmA_x = Array{S}(undef, N_ext)
    N_MfigM_pmA_y = Array{S}(undef, N_ext)
    N_MfigM_pmA_z = Array{S}(undef, N_ext)
    temp_N_M_x = Array{S}(undef, N_ext)
    temp_N_M_y = Array{S}(undef, N_ext)
    temp_N_M_z = Array{S}(undef, N_ext)
    N_MfigM = Array{S}(undef, 3)
    N_MfigM[1] = zero_q_1
    N_MfigM[2] = zero_q_1
    N_MfigM[3] = zero_q_1

    # rotations to and from Earth, Sun and Moon pole-oriented frames
    local αs = deg2rad(α_p_sun*one_t)
    local δs = deg2rad(δ_p_sun*one_t)
    local αm = eulang_t[1] - (pi/2)
    local δm = (pi/2) - eulang_t[2]
    local Wm = eulang_t[3]
    RotM = Array{S}(undef, 3, 3, 5) # space-fixed -> body-fixed coordinate transformations
    local RotM[:,:,ea] = c2t_jpl_de430(dsj2k)
    local RotM[:,:,su] = pole_rotation(αs, δs)
    local RotM[:,:,mo] = pole_rotation(αm, δm, Wm)
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
    J2_t = Array{S}(undef, 5)
    J2_t[su] = J2S_t
    J2_t[ea] = J2E_t
    # lunar torques: overall numerical factor in Eq. (44) of Folkner et al. (2014)
    local N_MfigM_figE_factor = 7.5*μ[ea]*J2E_t

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

    # 2nd-order lunar zonal and tesseral harmonics
    J2M_t = ( I_m_t[3,3] - ((I_m_t[1,1]+I_m_t[2,2])/2) )/(μ[mo]) # J_{2,M}*R_M^2
    C22M_t = ((I_m_t[2,2] - I_m_t[1,1])/(μ[mo]))/4 # C_{22,M}*R_M^2
    C21M_t = (-I_m_t[1,3])/(μ[mo]) # C_{21,M}*R_M^2
    S21M_t = (-I_m_t[3,2])/(μ[mo]) # S_{21,M}*R_M^2
    S22M_t = ((-I_m_t[2,1])/(μ[mo]))/2 # S_{22,M}*R_M^2
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
                                P_nm[i,j,1,1] = cos_ϕ[i,j]
                                cosϕ_dP_nm[i,j,1,1] = sin_ϕ[i,j]*lnm3[1] #+0 (second term in Eq. 183 from Moyer, 1971, vanishes when n=m)
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
                    # G_{i,j} = \sum_k R_{i,k} RotM{k,j}
                    Gc2p[i,j,1,1] = ((Rb2p[i,j,1,1]*RotM[1,1,j]) + (Rb2p[i,j,1,2]*RotM[2,1,j])) + (Rb2p[i,j,1,3]*RotM[3,1,j])
                    Gc2p[i,j,2,1] = ((Rb2p[i,j,2,1]*RotM[1,1,j]) + (Rb2p[i,j,2,2]*RotM[2,1,j])) + (Rb2p[i,j,2,3]*RotM[3,1,j])
                    Gc2p[i,j,3,1] = ((Rb2p[i,j,3,1]*RotM[1,1,j]) + (Rb2p[i,j,3,2]*RotM[2,1,j])) + (Rb2p[i,j,3,3]*RotM[3,1,j])
                    Gc2p[i,j,1,2] = ((Rb2p[i,j,1,1]*RotM[1,2,j]) + (Rb2p[i,j,1,2]*RotM[2,2,j])) + (Rb2p[i,j,1,3]*RotM[3,2,j])
                    Gc2p[i,j,2,2] = ((Rb2p[i,j,2,1]*RotM[1,2,j]) + (Rb2p[i,j,2,2]*RotM[2,2,j])) + (Rb2p[i,j,2,3]*RotM[3,2,j])
                    Gc2p[i,j,3,2] = ((Rb2p[i,j,3,1]*RotM[1,2,j]) + (Rb2p[i,j,3,2]*RotM[2,2,j])) + (Rb2p[i,j,3,3]*RotM[3,2,j])
                    Gc2p[i,j,1,3] = ((Rb2p[i,j,1,1]*RotM[1,3,j]) + (Rb2p[i,j,1,2]*RotM[2,3,j])) + (Rb2p[i,j,1,3]*RotM[3,3,j])
                    Gc2p[i,j,2,3] = ((Rb2p[i,j,2,1]*RotM[1,3,j]) + (Rb2p[i,j,2,2]*RotM[2,3,j])) + (Rb2p[i,j,2,3]*RotM[3,3,j])
                    Gc2p[i,j,3,3] = ((Rb2p[i,j,3,1]*RotM[1,3,j]) + (Rb2p[i,j,3,2]*RotM[2,3,j])) + (Rb2p[i,j,3,3]*RotM[3,3,j])
                    # compute cartesian coordinates of acceleration due to body figure in inertial frame (Moyer, 1971, eq. (169))
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

                    if j == mo
                        # compute torques acting upon the body-figure of the Moon due to external point masses
                        # Folkner et al. (2014), Eq. (43)
                        N_MfigM_pmA_x[i] = μ[i]*( (Y[i,j]*F_JCS_z[i,j]) - (Z[i,j]*F_JCS_y[i,j]) )
                        N_MfigM_pmA_y[i] = μ[i]*( (Z[i,j]*F_JCS_x[i,j]) - (X[i,j]*F_JCS_z[i,j]) )
                        N_MfigM_pmA_z[i] = μ[i]*( (X[i,j]*F_JCS_y[i,j]) - (Y[i,j]*F_JCS_x[i,j]) )
                        # expressions below have minus sign since N_MfigM_pmA_{x,y,z} have inverted signs in cross product
                        temp_N_M_x[i] = N_MfigM[1] - (μ[j]*N_MfigM_pmA_x[i])
                        N_MfigM[1] = temp_N_M_x[i]
                        temp_N_M_y[i] = N_MfigM[2] - (μ[j]*N_MfigM_pmA_y[i])
                        N_MfigM[2] = temp_N_M_y[i]
                        temp_N_M_z[i] = N_MfigM[3] - (μ[j]*N_MfigM_pmA_z[i])
                        N_MfigM[3] = temp_N_M_z[i]
                    end
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
                ###
                pn1[i,j] = zero_q_1
                X_t_pn1[i,j] = zero_q_1
                Y_t_pn1[i,j] = zero_q_1
                Z_t_pn1[i,j] = zero_q_1
                pNX_t_pn3[i,j] = zero_q_1
                pNY_t_pn3[i,j] = zero_q_1
                pNZ_t_pn3[i,j] = zero_q_1
                pNX_t_X[i,j] = zero_q_1
                pNY_t_Y[i,j] = zero_q_1
                pNZ_t_Z[i,j] = zero_q_1
            end # else (i != j)
        end
        pntempX[j] = zero_q_1
        pntempY[j] = zero_q_1
        pntempZ[j] = zero_q_1
    end

    Threads.@threads for j in 1:N
        for i in 1:N
            # i == j && continue
            if i == j
                continue
            else
                pNX_t_X[i,j] = newtonX[i]*X[i,j]
                pNY_t_Y[i,j] = newtonY[i]*Y[i,j]
                pNZ_t_Z[i,j] = newtonZ[i]*Z[i,j]
                pn1[i,j] = (  pn1t1_7[i,j]  +  0.5*( (pNX_t_X[i,j]+pNY_t_Y[i,j]) + pNZ_t_Z[i,j] )  )

                X_t_pn1[i,j] = newton_acc_X[i,j]*pn1[i,j]
                Y_t_pn1[i,j] = newton_acc_Y[i,j]*pn1[i,j]
                Z_t_pn1[i,j] = newton_acc_Z[i,j]*pn1[i,j]

                pNX_t_pn3[i,j] = newtonX[i]*pn3[i,j]
                pNY_t_pn3[i,j] = newtonY[i]*pn3[i,j]
                pNZ_t_pn3[i,j] = newtonZ[i]*pn3[i,j]

                termpnx = ( X_t_pn1[i,j] + (U_t_pn2[i,j]+pNX_t_pn3[i,j]) )
                sumpnx = pntempX[j] + termpnx
                pntempX[j] = sumpnx
                termpny = ( Y_t_pn1[i,j] + (V_t_pn2[i,j]+pNY_t_pn3[i,j]) )
                sumpny = pntempY[j] + termpny
                pntempY[j] = sumpny
                termpnz = ( Z_t_pn1[i,j] + (W_t_pn2[i,j]+pNZ_t_pn3[i,j]) )
                sumpnz = pntempZ[j] + termpnz
                pntempZ[j] = sumpnz
            end # else (i != j)
        end
        postNewtonX[j] = pntempX[j]*c_m2
        postNewtonY[j] = pntempY[j]*c_m2
        postNewtonZ[j] = pntempZ[j]*c_m2
    end

    #fill accelerations (post-Newtonian and extended body accelerations)
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

    Iω_x = (I_m_t[1,1]*q[6N+4]) + ((I_m_t[1,2]*q[6N+5]) + (I_m_t[1,3]*q[6N+6]))
    Iω_y = (I_m_t[2,1]*q[6N+4]) + ((I_m_t[2,2]*q[6N+5]) + (I_m_t[2,3]*q[6N+6]))
    Iω_z = (I_m_t[3,1]*q[6N+4]) + ((I_m_t[3,2]*q[6N+5]) + (I_m_t[3,3]*q[6N+6]))

    # ω × (I*ω)
    ωxIω_x = (q[6N+5]*Iω_z) - (q[6N+6]*Iω_y)
    ωxIω_y = (q[6N+6]*Iω_x) - (q[6N+4]*Iω_z)
    ωxIω_z = (q[6N+4]*Iω_y) - (q[6N+5]*Iω_x)

    # (dI/dt)*ω
    dIω_x = (dI_m_t[1,1]*q[6N+4]) + ((dI_m_t[1,2]*q[6N+5]) + (dI_m_t[1,3]*q[6N+6]))
    dIω_y = (dI_m_t[2,1]*q[6N+4]) + ((dI_m_t[2,2]*q[6N+5]) + (dI_m_t[2,3]*q[6N+6]))
    dIω_z = (dI_m_t[3,1]*q[6N+4]) + ((dI_m_t[3,2]*q[6N+5]) + (dI_m_t[3,3]*q[6N+6]))

    # Moon->Earth radial unit vector (inertial coordinates)
    er_EM_I_1 = X[ea,mo]/r_p1d2[ea,mo]
    er_EM_I_2 = Y[ea,mo]/r_p1d2[ea,mo]
    er_EM_I_3 = Z[ea,mo]/r_p1d2[ea,mo]

    # Earth pole unit vector (inertial coordinates)
    p_E_I_1 = RotM[3,1,ea]
    p_E_I_2 = RotM[3,2,ea]
    p_E_I_3 = RotM[3,3,ea]

    # transform er_EM_I_i, p_E_I_i to lunar mantle frame coordinates
    er_EM_1 = (RotM[1,1,mo]*er_EM_I_1) + ((RotM[1,2,mo]*er_EM_I_2) + (RotM[1,3,mo]*er_EM_I_3))
    er_EM_2 = (RotM[2,1,mo]*er_EM_I_1) + ((RotM[2,2,mo]*er_EM_I_2) + (RotM[2,3,mo]*er_EM_I_3))
    er_EM_3 = (RotM[3,1,mo]*er_EM_I_1) + ((RotM[3,2,mo]*er_EM_I_2) + (RotM[3,3,mo]*er_EM_I_3))
    p_E_1 = (RotM[1,1,mo]*p_E_I_1) + ((RotM[1,2,mo]*p_E_I_2) + (RotM[1,3,mo]*p_E_I_3))
    p_E_2 = (RotM[2,1,mo]*p_E_I_1) + ((RotM[2,2,mo]*p_E_I_2) + (RotM[2,3,mo]*p_E_I_3))
    p_E_3 = (RotM[3,1,mo]*p_E_I_1) + ((RotM[3,2,mo]*p_E_I_2) + (RotM[3,3,mo]*p_E_I_3))

    # evaluate Eq. (44) of Folkner et al. (2014) in lunar mantle frame coords
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

    one_minus_7sin2ϕEM = one_t - (7((sin_ϕ[mo,ea])^2))
    two_sinϕEM = 2sin_ϕ[mo,ea]

    N_MfigM_figE_factor_div_rEMp5 = (N_MfigM_figE_factor/(r_p1d2[mo,ea]^5))
    N_MfigM_figE_1 = N_MfigM_figE_factor_div_rEMp5*( (one_minus_7sin2ϕEM*er_EM_cross_I_er_EM_1) + (two_sinϕEM*(er_EM_cross_I_p_E_1+p_E_cross_I_er_EM_1)) - (0.4p_E_cross_I_p_E_1))
    N_MfigM_figE_2 = N_MfigM_figE_factor_div_rEMp5*( (one_minus_7sin2ϕEM*er_EM_cross_I_er_EM_2) + (two_sinϕEM*(er_EM_cross_I_p_E_2+p_E_cross_I_er_EM_2)) - (0.4p_E_cross_I_p_E_2))
    N_MfigM_figE_3 = N_MfigM_figE_factor_div_rEMp5*( (one_minus_7sin2ϕEM*er_EM_cross_I_er_EM_3) + (two_sinϕEM*(er_EM_cross_I_p_E_3+p_E_cross_I_er_EM_3)) - (0.4p_E_cross_I_p_E_3))

    # torques acting upon lunar body-figure due to external point masses: transform coordinates from inertial frame to lunar mantle frame
    N_1_LMF = (RotM[1,1,mo]*N_MfigM[1]) + ((RotM[1,2,mo]*N_MfigM[2]) + (RotM[1,3,mo]*N_MfigM[3]))
    N_2_LMF = (RotM[2,1,mo]*N_MfigM[1]) + ((RotM[2,2,mo]*N_MfigM[2]) + (RotM[2,3,mo]*N_MfigM[3]))
    N_3_LMF = (RotM[3,1,mo]*N_MfigM[1]) + ((RotM[3,2,mo]*N_MfigM[2]) + (RotM[3,3,mo]*N_MfigM[3]))

    # I*(dω/dt)
    I_dω_1 = (N_1_LMF + N_MfigM_figE_1) - (dIω_x + ωxIω_x)
    I_dω_2 = (N_2_LMF + N_MfigM_figE_2) - (dIω_y + ωxIω_y)
    I_dω_3 = (N_3_LMF + N_MfigM_figE_3) - (dIω_z + ωxIω_z)

    # lunar physical librations: Folkner et al. (2014), Eq. (14)
    dq[6N+1] = ((q[6N+4]*sin(q[6N+3])) + (q[6N+5]*cos(q[6N+3])) )/sin(q[6N+2])
    dq[6N+2] = (q[6N+4]*cos(q[6N+3])) - (q[6N+5]*sin(q[6N+3]))
    dq[6N+3] = q[6N+6] - (dq[6N+1]*cos(q[6N+2]))

    # lunar physical librations: Folkner et al. (2014), Eq. (34)
    # TODO: add lunar mantle-core boundary interaction
    dq[6N+4] = (inv_I_m_t[1,1]*I_dω_1) + ( (inv_I_m_t[1,2]*I_dω_2) + (inv_I_m_t[1,3]*I_dω_3) )
    dq[6N+5] = (inv_I_m_t[2,1]*I_dω_1) + ( (inv_I_m_t[2,2]*I_dω_2) + (inv_I_m_t[2,3]*I_dω_3) )
    dq[6N+6] = (inv_I_m_t[3,1]*I_dω_1) + ( (inv_I_m_t[3,2]*I_dω_2) + (inv_I_m_t[3,3]*I_dω_3) )

    # TT-TDB
    # TODO: implement TT-TDB integration
    dq[6N+7] = zero_q_1

    nothing
end

@taylorize function DE430!(dq, q, params, t)
    # N: number of bodies
    # jd0: initial Julian date
    local N, jd0 = params
    local S = eltype(q)
    local N_ext = 11 # number of bodies in extended-body accelerations
    local N_bwd = 11 # number of bodies in backward integration
    local params_bwd = (N_bwd, jd0)
    local qq_bwd = Taylor1.(constant_term.(  q[ union(nbodyind(N,1:N_bwd),lastindex(q)-6:lastindex(q)) ]), t.order )
    local dqq_bwd = similar(qq_bwd)
    # local xaux_bwd = similar(qq_bwd)
    # local jc = TaylorIntegration.__jetcoeffs!(Val(false), NBP_pN_A_J23E_J23M_J2S_threads!, t, qq_bwd, dqq_bwd, xaux_bwd, params_bwd)
    # local jc = TaylorIntegration.__jetcoeffs!(Val(true), NBP_pN_A_J23E_J23M_J2S_threads!, t, qq_bwd, dqq_bwd, xaux_bwd, params_bwd)
    # local jc = TaylorIntegration.jetcoeffs!(NBP_pN_A_J23E_J23M_J2S_threads!, t, qq_bwd, dqq_bwd, xaux_bwd, params_bwd)
    local jc = TaylorIntegration.jetcoeffs!(Val(NBP_pN_A_J23E_J23M_J2S_threads!), t, qq_bwd, dqq_bwd, params_bwd)
    local __t = Taylor1(t.order)
    local q_del_τ_M = qq_bwd(__t-τ_M)
    local q_del_τ_0 = qq_bwd(__t-τ_0p)
    local q_del_τ_1 = qq_bwd(__t-τ_1p)
    local q_del_τ_2 = qq_bwd(__t-τ_2p)

    local zero_q_1 = zero(q[1])
    local one_t = one(t)
    local dsj2k = t+(jd0-J2000) # days since J2000.0 (TDB)
    local eulang_t = qq_bwd[6N_bwd+1:6N_bwd+3]
    local eulang_t_del = q_del_τ_M[6N_bwd+1:6N_bwd+3]
    local I_m_t = ITM(q_del_τ_M, eulang_t_del) # matrix elements of lunar moment of inertia (Folkner et al. 2014, eq. 41)
    local dI_m_t = derivative.(I_m_t) # time-derivative of lunar mantle I at time t
    local inv_I_m_t = inv(I_m_t) # inverse of lunar mantle I matrix at time t

    # parameters related to speed of light, c
    local c_p2 = 29979.063823897606 # c^2 = 29979.063823897606 au^2/d^2
    local c_m2 = 3.3356611996764786e-5 # c^-2 = 3.3356611996764786e-5 d^2/au^2

    X = Array{S}(undef, N, N)
    Y = Array{S}(undef, N, N)
    Z = Array{S}(undef, N, N)

    r_p2 = Array{S}(undef, N, N)
    r_p3d2 = Array{S}(undef, N, N)
    r_p7d2 = Array{S}(undef, N, N)

    #Newtonian accelerations
    newtonX = Array{S}(undef, N)
    newtonY = Array{S}(undef, N)
    newtonZ = Array{S}(undef, N)

    newtonianCoeff = Array{S}(undef, N, N)

    #post-Newtonian stuff
    U = Array{S}(undef, N, N)
    V = Array{S}(undef, N, N)
    W = Array{S}(undef, N, N)

    _4U_m_3X = Array{S}(undef, N, N)
    _4V_m_3Y = Array{S}(undef, N, N)
    _4W_m_3Z = Array{S}(undef, N, N)

    UU = Array{S}(undef, N, N)
    VV = Array{S}(undef, N, N)
    WW = Array{S}(undef, N, N)

    r_p1d2 = Array{S}(undef, N, N)

    newtonianNb_Potential = Array{S}(undef, N)
    newtonian1b_Potential = Array{S}(undef, N, N)
    newtonianCoeff = Array{S}(undef, N, N)
    newton_acc_X = Array{S}(undef, N, N)
    newton_acc_Y = Array{S}(undef, N, N)
    newton_acc_Z = Array{S}(undef, N, N)

    v2 = Array{S}(undef, N)
    vi_dot_vj = Array{S}(undef, N, N)
    pn2 = Array{S}(undef, N, N)
    pn3 = Array{S}(undef, N, N)
    _4ϕj = Array{S}(undef, N, N)
    ϕi_plus_4ϕj = Array{S}(undef, N, N)
    sj2_plus_2si2_minus_4vivj = Array{S}(undef, N, N)
    ϕs_and_vs = Array{S}(undef, N, N)
    U_t_pn2 = Array{S}(undef, N, N)
    V_t_pn2 = Array{S}(undef, N, N)
    W_t_pn2 = Array{S}(undef, N, N)
    pn1t1_7 = Array{S}(undef, N, N)

    pntempX = Array{S}(undef, N)
    pntempY = Array{S}(undef, N)
    pntempZ = Array{S}(undef, N)
    pn1 = Array{S}(undef, N, N)
    X_t_pn1 = Array{S}(undef, N, N)
    Y_t_pn1 = Array{S}(undef, N, N)
    Z_t_pn1 = Array{S}(undef, N, N)
    pNX_t_pn3 = Array{S}(undef, N, N)
    pNY_t_pn3 = Array{S}(undef, N, N)
    pNZ_t_pn3 = Array{S}(undef, N, N)
    pNX_t_X = Array{S}(undef, N, N)
    pNY_t_Y = Array{S}(undef, N, N)
    pNZ_t_Z = Array{S}(undef, N, N)
    postNewtonX = Array{S}(undef, N)
    postNewtonY = Array{S}(undef, N)
    postNewtonZ = Array{S}(undef, N)

    # (Jn, Cmn, Smn) acceleration auxiliaries
    X_bf_1 = Array{S}(undef, N_ext, N_ext)
    Y_bf_1 = Array{S}(undef, N_ext, N_ext)
    Z_bf_1 = Array{S}(undef, N_ext, N_ext)
    X_bf_2 = Array{S}(undef, N_ext, N_ext)
    Y_bf_2 = Array{S}(undef, N_ext, N_ext)
    Z_bf_2 = Array{S}(undef, N_ext, N_ext)
    X_bf_3 = Array{S}(undef, N_ext, N_ext)
    Y_bf_3 = Array{S}(undef, N_ext, N_ext)
    Z_bf_3 = Array{S}(undef, N_ext, N_ext)
    X_bf = Array{S}(undef, N_ext, N_ext)
    Y_bf = Array{S}(undef, N_ext, N_ext)
    Z_bf = Array{S}(undef, N_ext, N_ext)
    F_JCS_x = Array{S}(undef, N_ext, N_ext)
    F_JCS_y = Array{S}(undef, N_ext, N_ext)
    F_JCS_z = Array{S}(undef, N_ext, N_ext)
    temp_accX_j = Array{S}(undef, N_ext, N_ext)
    temp_accY_j = Array{S}(undef, N_ext, N_ext)
    temp_accZ_j = Array{S}(undef, N_ext, N_ext)
    temp_accX_i = Array{S}(undef, N_ext, N_ext)
    temp_accY_i = Array{S}(undef, N_ext, N_ext)
    temp_accZ_i = Array{S}(undef, N_ext, N_ext)
    sin_ϕ = Array{S}(undef, N_ext, N_ext)
    cos_ϕ = Array{S}(undef, N_ext, N_ext)
    sin_λ = Array{S}(undef, N_ext, N_ext)
    cos_λ = Array{S}(undef, N_ext, N_ext)
    r_xy = Array{S}(undef, N_ext, N_ext)
    r_p4 = Array{S}(undef, N_ext, N_ext)
    P_n = Array{S}(undef, N_ext, N_ext, maximum(n1SEM)+1)
    dP_n = Array{S}(undef, N_ext, N_ext, maximum(n1SEM)+1)
    temp_fjξ = Array{S}(undef, N_ext, N_ext, maximum(n1SEM)+1)
    temp_fjζ = Array{S}(undef, N_ext, N_ext, maximum(n1SEM)+1)
    temp_rn = Array{S}(undef, N_ext, N_ext, maximum(n1SEM)+1)
    F_CS_ξ_36 = Array{S}(undef, N_ext, N_ext)
    F_CS_η_36 = Array{S}(undef, N_ext, N_ext)
    F_CS_ζ_36 = Array{S}(undef, N_ext, N_ext)
    F_J_ξ_36 = Array{S}(undef, N_ext, N_ext)
    F_J_ζ_36 = Array{S}(undef, N_ext, N_ext)
    sin_mλ = Array{S}(undef, N_ext, N_ext, n1SEM[mo])
    cos_mλ = Array{S}(undef, N_ext, N_ext, n1SEM[mo])
    secϕ_P_nm = Array{S}(undef, N_ext, N_ext, n1SEM[mo]+1, n1SEM[mo]+1)
    P_nm = Array{S}(undef, N_ext, N_ext, n1SEM[mo]+1, n1SEM[mo]+1)
    cosϕ_dP_nm = Array{S}(undef, N_ext, N_ext, n1SEM[mo]+1, n1SEM[mo]+1)
    F_J_ξ = Array{S}(undef, N_ext, N_ext)
    F_J_η = Array{S}(undef, N_ext, N_ext)
    F_J_ζ = Array{S}(undef, N_ext, N_ext)
    F_CS_ξ = Array{S}(undef, N_ext, N_ext)
    F_CS_η = Array{S}(undef, N_ext, N_ext)
    F_CS_ζ = Array{S}(undef, N_ext, N_ext)
    F_JCS_ξ = Array{S}(undef, N_ext, N_ext)
    F_JCS_η = Array{S}(undef, N_ext, N_ext)
    F_JCS_ζ = Array{S}(undef, N_ext, N_ext)
    Rb2p = Array{S}(undef, N_ext, N_ext, 3, 3) #R matrix body-fixed to "primed" ξηζ frame (Moyer, 1971, eq. 161)
    Gc2p = Array{S}(undef, N_ext, N_ext, 3, 3) #G matrix "space-fixed" to "primed" ξηζ frame (Moyer, 1971, eq. 163)

    # extended-body accelerations
    accX = Array{S}(undef, N_ext)
    accY = Array{S}(undef, N_ext)
    accZ = Array{S}(undef, N_ext)

    # lunar torques: Folkner et al. (2014), Eq. (43)
    N_MfigM_pmA_x = Array{S}(undef, N_ext)
    N_MfigM_pmA_y = Array{S}(undef, N_ext)
    N_MfigM_pmA_z = Array{S}(undef, N_ext)
    temp_N_M_x = Array{S}(undef, N_ext)
    temp_N_M_y = Array{S}(undef, N_ext)
    temp_N_M_z = Array{S}(undef, N_ext)
    N_MfigM = Array{S}(undef, 3)
    N_MfigM[1] = zero_q_1
    N_MfigM[2] = zero_q_1
    N_MfigM[3] = zero_q_1

    # rotations to and from Earth, Sun and Moon pole-oriented frames
    local αs = deg2rad(α_p_sun*one_t)
    local δs = deg2rad(δ_p_sun*one_t)
    local αm = eulang_t[1] - (pi/2)
    local δm = (pi/2) - eulang_t[2]
    local Wm = eulang_t[3]
    RotM = Array{S}(undef, 3, 3, 5) # space-fixed -> body-fixed coordinate transformations
    local RotM[:,:,ea] = c2t_jpl_de430(dsj2k)
    local RotM[:,:,su] = pole_rotation(αs, δs)
    local RotM[:,:,mo] = pole_rotation(αm, δm, Wm)
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
    J2_t = Array{S}(undef, 5)
    J2_t[su] = J2S_t
    J2_t[ea] = J2E_t
    # tidal acceleration auxiliaries
    local μ_mo_div_μ_ea = μ[mo]/μ[ea]
    local tid_num_coeff = 1.5*(1.0 + μ_mo_div_μ_ea)
    # Folkner et al. (2014), Eq. (31)
    # time-delayed geocentric Moon position
    local q_ME_τ_0 = q_del_τ_0[3mo-2:3mo] .- q_del_τ_0[3ea-2:3ea]
    local q_ME_τ_1 = q_del_τ_1[3mo-2:3mo] .- q_del_τ_1[3ea-2:3ea]
    local q_ME_τ_2 = q_del_τ_2[3mo-2:3mo] .- q_del_τ_2[3ea-2:3ea]
    # time-delayed geocentric Sun position
    local q_SE_τ_0 = q_del_τ_0[3su-2:3su] .- q_del_τ_0[3ea-2:3ea]
    local q_SE_τ_1 = q_del_τ_1[3su-2:3su] .- q_del_τ_1[3ea-2:3ea]
    local q_SE_τ_2 = q_del_τ_2[3su-2:3su] .- q_del_τ_2[3ea-2:3ea]
    # R3X: geocentric space-fixed -> rotational time-delay -> geocentric Earth true-equator-of-date frame
    local R30 = RotM[:,:,ea] # ###*Rz(-ω_E*τ_0) == Id(3x3), since τ_0=0
    local R31 = Rz(-ω_E*τ_1)*RotM[:,:,ea]
    local R32 = Rz(-ω_E*τ_2)*RotM[:,:,ea]
    # tidal accelerations
    local r_star_M_0 = R30*q_ME_τ_0 # r-star 0, Moon
    local r_star_M_1 = R31*q_ME_τ_1 # r-star 1, Moon
    local r_star_M_2 = R32*q_ME_τ_2 # r-star 2, Moon
    local r_star_S_0 = R30*q_SE_τ_0 # r-star 0, Sun
    local r_star_S_1 = R31*q_SE_τ_1 # r-star 1, Sun
    local r_star_S_2 = R32*q_SE_τ_2 # r-star 2, Sun
    # lunar torques: overall numerical factor in Eq. (44) of Folkner et al. (2014)
    local N_MfigM_figE_factor = 7.5*μ[ea]*J2E_t

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

    # 2nd-order lunar zonal and tesseral harmonics
    J2M_t = ( I_m_t[3,3] - ((I_m_t[1,1]+I_m_t[2,2])/2) )/(μ[mo]) # J_{2,M}*R_M^2
    C22M_t = ((I_m_t[2,2] - I_m_t[1,1])/(μ[mo]))/4 # C_{22,M}*R_M^2
    C21M_t = (-I_m_t[1,3])/(μ[mo]) # C_{21,M}*R_M^2
    S21M_t = (-I_m_t[3,2])/(μ[mo]) # S_{21,M}*R_M^2
    S22M_t = ((-I_m_t[2,1])/(μ[mo]))/2 # S_{22,M}*R_M^2
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
                                P_nm[i,j,1,1] = cos_ϕ[i,j]
                                cosϕ_dP_nm[i,j,1,1] = sin_ϕ[i,j]*lnm3[1] #+0 (second term in Eq. 183 from Moyer, 1971, vanishes when n=m)
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
                    # G_{i,j} = \sum_k R_{i,k} RotM{k,j}
                    Gc2p[i,j,1,1] = ((Rb2p[i,j,1,1]*RotM[1,1,j]) + (Rb2p[i,j,1,2]*RotM[2,1,j])) + (Rb2p[i,j,1,3]*RotM[3,1,j])
                    Gc2p[i,j,2,1] = ((Rb2p[i,j,2,1]*RotM[1,1,j]) + (Rb2p[i,j,2,2]*RotM[2,1,j])) + (Rb2p[i,j,2,3]*RotM[3,1,j])
                    Gc2p[i,j,3,1] = ((Rb2p[i,j,3,1]*RotM[1,1,j]) + (Rb2p[i,j,3,2]*RotM[2,1,j])) + (Rb2p[i,j,3,3]*RotM[3,1,j])
                    Gc2p[i,j,1,2] = ((Rb2p[i,j,1,1]*RotM[1,2,j]) + (Rb2p[i,j,1,2]*RotM[2,2,j])) + (Rb2p[i,j,1,3]*RotM[3,2,j])
                    Gc2p[i,j,2,2] = ((Rb2p[i,j,2,1]*RotM[1,2,j]) + (Rb2p[i,j,2,2]*RotM[2,2,j])) + (Rb2p[i,j,2,3]*RotM[3,2,j])
                    Gc2p[i,j,3,2] = ((Rb2p[i,j,3,1]*RotM[1,2,j]) + (Rb2p[i,j,3,2]*RotM[2,2,j])) + (Rb2p[i,j,3,3]*RotM[3,2,j])
                    Gc2p[i,j,1,3] = ((Rb2p[i,j,1,1]*RotM[1,3,j]) + (Rb2p[i,j,1,2]*RotM[2,3,j])) + (Rb2p[i,j,1,3]*RotM[3,3,j])
                    Gc2p[i,j,2,3] = ((Rb2p[i,j,2,1]*RotM[1,3,j]) + (Rb2p[i,j,2,2]*RotM[2,3,j])) + (Rb2p[i,j,2,3]*RotM[3,3,j])
                    Gc2p[i,j,3,3] = ((Rb2p[i,j,3,1]*RotM[1,3,j]) + (Rb2p[i,j,3,2]*RotM[2,3,j])) + (Rb2p[i,j,3,3]*RotM[3,3,j])
                    # compute cartesian coordinates of acceleration due to body figure in inertial frame (Moyer, 1971, eq. (169))
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

                    if j == mo
                        # compute torques acting upon the body-figure of the Moon due to external point masses
                        # Folkner et al. (2014), Eq. (43)
                        N_MfigM_pmA_x[i] = μ[i]*( (Y[i,j]*F_JCS_z[i,j]) - (Z[i,j]*F_JCS_y[i,j]) )
                        N_MfigM_pmA_y[i] = μ[i]*( (Z[i,j]*F_JCS_x[i,j]) - (X[i,j]*F_JCS_z[i,j]) )
                        N_MfigM_pmA_z[i] = μ[i]*( (X[i,j]*F_JCS_y[i,j]) - (Y[i,j]*F_JCS_x[i,j]) )
                        # expressions below have minus sign since N_MfigM_pmA_{x,y,z} have inverted signs in cross product
                        temp_N_M_x[i] = N_MfigM[1] - (μ[j]*N_MfigM_pmA_x[i])
                        N_MfigM[1] = temp_N_M_x[i]
                        temp_N_M_y[i] = N_MfigM[2] - (μ[j]*N_MfigM_pmA_y[i])
                        N_MfigM[2] = temp_N_M_y[i]
                        temp_N_M_z[i] = N_MfigM[3] - (μ[j]*N_MfigM_pmA_z[i])
                        N_MfigM[3] = temp_N_M_z[i]
                    end
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
                ###
                pn1[i,j] = zero_q_1
                X_t_pn1[i,j] = zero_q_1
                Y_t_pn1[i,j] = zero_q_1
                Z_t_pn1[i,j] = zero_q_1
                pNX_t_pn3[i,j] = zero_q_1
                pNY_t_pn3[i,j] = zero_q_1
                pNZ_t_pn3[i,j] = zero_q_1
                pNX_t_X[i,j] = zero_q_1
                pNY_t_Y[i,j] = zero_q_1
                pNZ_t_Z[i,j] = zero_q_1
            end # else (i != j)
        end
        pntempX[j] = zero_q_1
        pntempY[j] = zero_q_1
        pntempZ[j] = zero_q_1
    end

    Threads.@threads for j in 1:N
        for i in 1:N
            # i == j && continue
            if i == j
                continue
            else
                pNX_t_X[i,j] = newtonX[i]*X[i,j]
                pNY_t_Y[i,j] = newtonY[i]*Y[i,j]
                pNZ_t_Z[i,j] = newtonZ[i]*Z[i,j]
                pn1[i,j] = (  pn1t1_7[i,j]  +  0.5*( (pNX_t_X[i,j]+pNY_t_Y[i,j]) + pNZ_t_Z[i,j] )  )

                X_t_pn1[i,j] = newton_acc_X[i,j]*pn1[i,j]
                Y_t_pn1[i,j] = newton_acc_Y[i,j]*pn1[i,j]
                Z_t_pn1[i,j] = newton_acc_Z[i,j]*pn1[i,j]

                pNX_t_pn3[i,j] = newtonX[i]*pn3[i,j]
                pNY_t_pn3[i,j] = newtonY[i]*pn3[i,j]
                pNZ_t_pn3[i,j] = newtonZ[i]*pn3[i,j]

                termpnx = ( X_t_pn1[i,j] + (U_t_pn2[i,j]+pNX_t_pn3[i,j]) )
                sumpnx = pntempX[j] + termpnx
                pntempX[j] = sumpnx
                termpny = ( Y_t_pn1[i,j] + (V_t_pn2[i,j]+pNY_t_pn3[i,j]) )
                sumpny = pntempY[j] + termpny
                pntempY[j] = sumpny
                termpnz = ( Z_t_pn1[i,j] + (W_t_pn2[i,j]+pNZ_t_pn3[i,j]) )
                sumpnz = pntempZ[j] + termpnz
                pntempZ[j] = sumpnz
            end # else (i != j)
        end
        postNewtonX[j] = pntempX[j]*c_m2
        postNewtonY[j] = pntempY[j]*c_m2
        postNewtonZ[j] = pntempZ[j]*c_m2
    end

    # compute tidal acceleration of the Moon due to tides raised on Earth by the Sun and Moon
    # Folkner et al. (2014) Eq. (32)
    # X_bf[mo,ea] are geocentric, Earth-fixed position coordinates in "unprimed" (inertial) frame of perturbed body (Moon)
    x0s_M = r_star_M_0[1]
    y0s_M = r_star_M_0[2]
    z0s_M = r_star_M_0[3]
    ρ0s2_M = (x0s_M^2) + (y0s_M^2)
    ρ0s_M = sqrt(ρ0s2_M)
    z0s2_M = z0s_M^2
    r0s2_M = ρ0s2_M + z0s2_M
    r0s_M = sqrt(r0s2_M)
    r0s5_M = r0s_M^5

    x0s_S = r_star_S_0[1]
    y0s_S = r_star_S_0[2]
    z0s_S = r_star_S_0[3]
    ρ0s2_S = (x0s_S^2) + (y0s_S^2)
    ρ0s_S = sqrt(ρ0s2_S)
    z0s2_S = z0s_S^2
    r0s2_S = ρ0s2_S + z0s2_S
    r0s_S = sqrt(r0s2_S)
    r0s5_S = r0s_S^5

    coeff0_M = r0s2_M - 5( ( ((Z_bf[mo,ea]*r_star_M_0[3])^2) + 0.5((r_xy[mo,ea]*ρ0s_M)^2) )/r_p2[mo,ea] )
    coeff0_S = r0s2_S - 5( ( ((Z_bf[mo,ea]*r_star_S_0[3])^2) + 0.5((r_xy[mo,ea]*ρ0s_S)^2) )/r_p2[mo,ea] )

    k_20E_div_r0s5_M = k_20E/r0s5_M
    k_20E_div_r0s5_S = k_20E/r0s5_S

    a_tid_0_M_x = k_20E_div_r0s5_M*((ρ0s2_M + coeff0_M)*X_bf[mo,ea])
    a_tid_0_M_y = k_20E_div_r0s5_M*((ρ0s2_M + coeff0_M)*Y_bf[mo,ea])
    a_tid_0_M_z = k_20E_div_r0s5_M*(((2z0s2_M) + coeff0_M)*Z_bf[mo,ea])
    a_tid_0_S_x = k_20E_div_r0s5_S*((ρ0s2_S + coeff0_S)*X_bf[mo,ea])
    a_tid_0_S_y = k_20E_div_r0s5_S*((ρ0s2_S + coeff0_S)*Y_bf[mo,ea])
    a_tid_0_S_z = k_20E_div_r0s5_S*(((2z0s2_S) + coeff0_S)*Z_bf[mo,ea])

    x1s_M = r_star_M_1[1]
    y1s_M = r_star_M_1[2]
    z1s_M = r_star_M_1[3]
    ρ1s2_M = (x1s_M^2) + (y1s_M^2)
    ρ1s_M = sqrt(ρ1s2_M)
    z1s2_M = z1s_M^2
    r1s2_M = ρ1s2_M + z1s2_M
    r1s_M = sqrt(r1s2_M)
    r1s5_M = r1s_M^5

    x1s_S = r_star_S_1[1]
    y1s_S = r_star_S_1[2]
    z1s_S = r_star_S_1[3]
    ρ1s2_S = (x1s_S^2) + (y1s_S^2)
    ρ1s_S = sqrt(ρ1s2_S)
    z1s2_S = z1s_S^2
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

    a_tid_1_M_x = k_21E_div_r1s5_M*((2coeff2_1_M*r_star_M_1[1]) - (coeff3_1_M*X_bf[mo,ea]))
    a_tid_1_M_y = k_21E_div_r1s5_M*((2coeff2_1_M*r_star_M_1[2]) - (coeff3_1_M*Y_bf[mo,ea]))
    a_tid_1_M_z = k_21E_div_r1s5_M*((2coeff1_1_M*r_star_M_1[3]) - (coeff3_1_M*Z_bf[mo,ea]))
    a_tid_1_S_x = k_21E_div_r1s5_S*((2coeff2_1_S*r_star_S_1[1]) - (coeff3_1_S*X_bf[mo,ea]))
    a_tid_1_S_y = k_21E_div_r1s5_S*((2coeff2_1_S*r_star_S_1[2]) - (coeff3_1_S*Y_bf[mo,ea]))
    a_tid_1_S_z = k_21E_div_r1s5_S*((2coeff1_1_S*r_star_S_1[3]) - (coeff3_1_S*Z_bf[mo,ea]))

    x2s_M = r_star_M_2[1]
    y2s_M = r_star_M_2[2]
    z2s_M = r_star_M_2[3]
    ρ2s2_M = (x2s_M^2) + (y2s_M^2)
    ρ2s_M = sqrt(ρ2s2_M)
    z2s2_M = z2s_M^2
    r2s2_M = ρ2s2_M + z2s2_M
    r2s_M = sqrt(r2s2_M)
    r2s5_M = r2s_M^5

    x2s_S = r_star_S_2[1]
    y2s_S = r_star_S_2[2]
    z2s_S = r_star_S_2[3]
    ρ2s2_S = (x2s_S^2) + (y2s_S^2)
    ρ2s_S = sqrt(ρ2s2_S)
    z2s2_S = z2s_S^2
    r2s2_S = ρ2s2_S + z2s2_S
    r2s_S = sqrt(r2s2_S)
    r2s5_S = r2s_S^5

    coeff1_2_M = (X_bf[mo,ea]*r_star_M_2[1]) + (Y_bf[mo,ea]*r_star_M_2[2])
    coeff1_2_S = (X_bf[mo,ea]*r_star_S_2[1]) + (Y_bf[mo,ea]*r_star_S_2[2])

    coeff3_2_M = 5( (coeff1_2_M^2) - 0.5((r_xy[mo,ea]^2)*ρ2s2_M) )/r_p2[mo,ea]
    coeff3_2_S = 5( (coeff1_2_S^2) - 0.5((r_xy[mo,ea]^2)*ρ2s2_S) )/r_p2[mo,ea]

    k_22E_div_r2s5_M = k_22E/r2s5_M
    k_22E_div_r2s5_S = k_22E/r2s5_S

    a_tid_2_M_x = k_22E_div_r2s5_M*(  (2coeff1_2_M*r_star_M_2[1])-(ρ2s2_M+coeff3_2_M)*X_bf[mo,ea]  )
    a_tid_2_M_y = k_22E_div_r2s5_M*(  (2coeff1_2_M*r_star_M_2[2])-(ρ2s2_M+coeff3_2_M)*Y_bf[mo,ea]  )
    a_tid_2_M_z = k_22E_div_r2s5_M*(  -coeff3_2_M*Z_bf[mo,ea]  )
    a_tid_2_S_x = k_22E_div_r2s5_S*(  (2coeff1_2_S*r_star_S_2[1])-(ρ2s2_S+coeff3_2_S)*X_bf[mo,ea]  )
    a_tid_2_S_y = k_22E_div_r2s5_S*(  (2coeff1_2_S*r_star_S_2[2])-(ρ2s2_S+coeff3_2_S)*Y_bf[mo,ea]  )
    a_tid_2_S_z = k_22E_div_r2s5_S*(  -coeff3_2_S*Z_bf[mo,ea]  )

    RE_div_r_p5 = (RE_au/r_p1d2[mo,ea])^5
    aux_tidacc = tid_num_coeff*RE_div_r_p5
    a_tidal_coeff_M = μ[mo]*aux_tidacc
    a_tidal_coeff_S = μ[su]*aux_tidacc

    # add contributions from long-period, diurnal and semi-diurnal tides (true-of-date coordinates)
    a_tidal_tod_x = (a_tidal_coeff_M*((a_tid_0_M_x+a_tid_1_M_x)+a_tid_2_M_x)) + (a_tidal_coeff_S*((a_tid_0_S_x+a_tid_1_S_x)+a_tid_2_S_x))
    a_tidal_tod_y = (a_tidal_coeff_M*((a_tid_0_M_y+a_tid_1_M_y)+a_tid_2_M_y)) + (a_tidal_coeff_S*((a_tid_0_S_y+a_tid_1_S_y)+a_tid_2_S_y))
    a_tidal_tod_z = (a_tidal_coeff_M*((a_tid_0_M_z+a_tid_1_M_z)+a_tid_2_M_z)) + (a_tidal_coeff_S*((a_tid_0_S_z+a_tid_1_S_z)+a_tid_2_S_z))

    # transform from geocentric Earth-true-equator-of-date coordinates to geocentric mean equator of J2000.0 coordinates
    a_tidal_x = ((RotM[1,1,ea]*a_tidal_tod_x)+(RotM[2,1,ea]*a_tidal_tod_y)) + (RotM[3,1,ea]*a_tidal_tod_z)
    a_tidal_y = ((RotM[1,2,ea]*a_tidal_tod_x)+(RotM[2,2,ea]*a_tidal_tod_y)) + (RotM[3,2,ea]*a_tidal_tod_z)
    a_tidal_z = ((RotM[1,3,ea]*a_tidal_tod_x)+(RotM[2,3,ea]*a_tidal_tod_y)) + (RotM[3,3,ea]*a_tidal_tod_z)

    # add tidal acceleration to Moon's acceleration due to extended-body effects
    accX_mo_tides = accX[mo] + a_tidal_x
    accY_mo_tides = accY[mo] + a_tidal_y
    accZ_mo_tides = accZ[mo] + a_tidal_z
    accX[mo] = accX_mo_tides
    accY[mo] = accY_mo_tides
    accZ[mo] = accZ_mo_tides

    #fill accelerations (post-Newtonian and extended body accelerations)
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

    Iω_x = (I_m_t[1,1]*q[6N+4]) + ((I_m_t[1,2]*q[6N+5]) + (I_m_t[1,3]*q[6N+6]))
    Iω_y = (I_m_t[2,1]*q[6N+4]) + ((I_m_t[2,2]*q[6N+5]) + (I_m_t[2,3]*q[6N+6]))
    Iω_z = (I_m_t[3,1]*q[6N+4]) + ((I_m_t[3,2]*q[6N+5]) + (I_m_t[3,3]*q[6N+6]))

    # ω × (I*ω)
    ωxIω_x = (q[6N+5]*Iω_z) - (q[6N+6]*Iω_y)
    ωxIω_y = (q[6N+6]*Iω_x) - (q[6N+4]*Iω_z)
    ωxIω_z = (q[6N+4]*Iω_y) - (q[6N+5]*Iω_x)

    # (dI/dt)*ω
    dIω_x = (dI_m_t[1,1]*q[6N+4]) + ((dI_m_t[1,2]*q[6N+5]) + (dI_m_t[1,3]*q[6N+6]))
    dIω_y = (dI_m_t[2,1]*q[6N+4]) + ((dI_m_t[2,2]*q[6N+5]) + (dI_m_t[2,3]*q[6N+6]))
    dIω_z = (dI_m_t[3,1]*q[6N+4]) + ((dI_m_t[3,2]*q[6N+5]) + (dI_m_t[3,3]*q[6N+6]))

    # Moon->Earth radial unit vector (inertial coordinates)
    er_EM_I_1 = X[ea,mo]/r_p1d2[ea,mo]
    er_EM_I_2 = Y[ea,mo]/r_p1d2[ea,mo]
    er_EM_I_3 = Z[ea,mo]/r_p1d2[ea,mo]

    # Earth pole unit vector (inertial coordinates)
    p_E_I_1 = RotM[3,1,ea]
    p_E_I_2 = RotM[3,2,ea]
    p_E_I_3 = RotM[3,3,ea]

    # transform er_EM_I_i, p_E_I_i to lunar mantle frame coordinates
    er_EM_1 = (RotM[1,1,mo]*er_EM_I_1) + ((RotM[1,2,mo]*er_EM_I_2) + (RotM[1,3,mo]*er_EM_I_3))
    er_EM_2 = (RotM[2,1,mo]*er_EM_I_1) + ((RotM[2,2,mo]*er_EM_I_2) + (RotM[2,3,mo]*er_EM_I_3))
    er_EM_3 = (RotM[3,1,mo]*er_EM_I_1) + ((RotM[3,2,mo]*er_EM_I_2) + (RotM[3,3,mo]*er_EM_I_3))
    p_E_1 = (RotM[1,1,mo]*p_E_I_1) + ((RotM[1,2,mo]*p_E_I_2) + (RotM[1,3,mo]*p_E_I_3))
    p_E_2 = (RotM[2,1,mo]*p_E_I_1) + ((RotM[2,2,mo]*p_E_I_2) + (RotM[2,3,mo]*p_E_I_3))
    p_E_3 = (RotM[3,1,mo]*p_E_I_1) + ((RotM[3,2,mo]*p_E_I_2) + (RotM[3,3,mo]*p_E_I_3))

    # evaluate Eq. (44) of Folkner et al. (2014) in lunar mantle frame coords
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

    one_minus_7sin2ϕEM = one_t - (7((sin_ϕ[mo,ea])^2))
    two_sinϕEM = 2sin_ϕ[mo,ea]

    N_MfigM_figE_factor_div_rEMp5 = (N_MfigM_figE_factor/(r_p1d2[mo,ea]^5))
    N_MfigM_figE_1 = N_MfigM_figE_factor_div_rEMp5*( (one_minus_7sin2ϕEM*er_EM_cross_I_er_EM_1) + (two_sinϕEM*(er_EM_cross_I_p_E_1+p_E_cross_I_er_EM_1)) - (0.4p_E_cross_I_p_E_1))
    N_MfigM_figE_2 = N_MfigM_figE_factor_div_rEMp5*( (one_minus_7sin2ϕEM*er_EM_cross_I_er_EM_2) + (two_sinϕEM*(er_EM_cross_I_p_E_2+p_E_cross_I_er_EM_2)) - (0.4p_E_cross_I_p_E_2))
    N_MfigM_figE_3 = N_MfigM_figE_factor_div_rEMp5*( (one_minus_7sin2ϕEM*er_EM_cross_I_er_EM_3) + (two_sinϕEM*(er_EM_cross_I_p_E_3+p_E_cross_I_er_EM_3)) - (0.4p_E_cross_I_p_E_3))

    # torques acting upon lunar body-figure due to external point masses: transform coordinates from inertial frame to lunar mantle frame
    N_1_LMF = (RotM[1,1,mo]*N_MfigM[1]) + ((RotM[1,2,mo]*N_MfigM[2]) + (RotM[1,3,mo]*N_MfigM[3]))
    N_2_LMF = (RotM[2,1,mo]*N_MfigM[1]) + ((RotM[2,2,mo]*N_MfigM[2]) + (RotM[2,3,mo]*N_MfigM[3]))
    N_3_LMF = (RotM[3,1,mo]*N_MfigM[1]) + ((RotM[3,2,mo]*N_MfigM[2]) + (RotM[3,3,mo]*N_MfigM[3]))

    # I*(dω/dt)
    I_dω_1 = (N_1_LMF + N_MfigM_figE_1) - (dIω_x + ωxIω_x)
    I_dω_2 = (N_2_LMF + N_MfigM_figE_2) - (dIω_y + ωxIω_y)
    I_dω_3 = (N_3_LMF + N_MfigM_figE_3) - (dIω_z + ωxIω_z)

    # lunar physical librations: Folkner et al. (2014), Eq. (14)
    dq[6N+1] = ((q[6N+4]*sin(q[6N+3])) + (q[6N+5]*cos(q[6N+3])) )/sin(q[6N+2])
    dq[6N+2] = (q[6N+4]*cos(q[6N+3])) - (q[6N+5]*sin(q[6N+3]))
    dq[6N+3] = q[6N+6] - (dq[6N+1]*cos(q[6N+2]))

    # lunar physical librations: Folkner et al. (2014), Eq. (34)
    # TODO: add lunar mantle-core boundary interaction
    dq[6N+4] = (inv_I_m_t[1,1]*I_dω_1) + ( (inv_I_m_t[1,2]*I_dω_2) + (inv_I_m_t[1,3]*I_dω_3) )
    dq[6N+5] = (inv_I_m_t[2,1]*I_dω_1) + ( (inv_I_m_t[2,2]*I_dω_2) + (inv_I_m_t[2,3]*I_dω_3) )
    dq[6N+6] = (inv_I_m_t[3,1]*I_dω_1) + ( (inv_I_m_t[3,2]*I_dω_2) + (inv_I_m_t[3,3]*I_dω_3) )

    # TT-TDB
    # TODO: implement TT-TDB integration
    dq[6N+7] = zero_q_1

    nothing
end
