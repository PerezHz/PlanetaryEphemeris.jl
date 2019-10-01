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
function NBP_pN_A_J23E_J23M_J2S!(dq, q, params, t)
    local S = eltype(q[1])
    local N = Int((length(q))/6) # number of bodies, including NEA
    local _1_to_N = Base.OneTo(N) # iterator over all bodies

    #TODO: handle appropiately @taylorize'd version with postnewton_iter>1
    local postnewton_iter = 1 # number of iterations of post-Newtonian subroutine
    local j2_body_index = [su, ea, mo] # indices of bodies with J2 flattening (note: Earth and Moon also have J3)

    # parameters related to speed of light, c
    local c_p2 = 29979.063823897606 # c^2 = 29979.063823897606 au^2/d^2
    local c_m2 = 3.3356611996764786e-5 # c^-2 = 3.3356611996764786e-5 d^2/au^2

    local zero_q_1 = zero(q[1])

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

    for j in _1_to_N
        for i in _1_to_N
            # i == j && continue
            if i == j
            else
                pn2[i,j] = zero_q_1
                pn3[i,j] = zero_q_1
                _4ϕj[i,j] = zero_q_1
                ϕi_plus_4ϕj[i,j] = zero_q_1
                sj2_plus_2si2_minus_4vivj[i,j] = zero_q_1
                ϕs_and_vs[i,j] = zero_q_1
                UU[i,j] = zero_q_1
                VV[i,j] = zero_q_1
                WW[i,j] = zero_q_1
                U_t_pn2[i,j] = zero_q_1
                V_t_pn2[i,j] = zero_q_1
                W_t_pn2[i,j] = zero_q_1
                vi_dot_vj[i,j] = zero_q_1
                newton_acc_X[i,j] = zero_q_1
                newton_acc_Y[i,j] = zero_q_1
                newton_acc_Z[i,j] = zero_q_1
                pn1t1_7[i,j] = zero_q_1
            end
        end
    end

    # J2 acceleration auxiliaries
    t31 = Array{Taylor1{S}}(undef, N, N)
    t32 = Array{Taylor1{S}}(undef, N, N)
    t33 = Array{Taylor1{S}}(undef, N, N)
    r_sin_ϕ = Array{Taylor1{S}}(undef, N, N)
    F_J2_x = Array{Taylor1{S}}(undef, N, N)
    F_J2_y = Array{Taylor1{S}}(undef, N, N)
    F_J2_z = Array{Taylor1{S}}(undef, N, N)
    F_J2_x1 = Array{Taylor1{S}}(undef, N, N)
    F_J2_y1 = Array{Taylor1{S}}(undef, N, N)
    F_J2_z1 = Array{Taylor1{S}}(undef, N, N)
    F_J2_x2 = Array{Taylor1{S}}(undef, N, N)
    F_J2_y2 = Array{Taylor1{S}}(undef, N, N)
    F_J2_z2 = Array{Taylor1{S}}(undef, N, N)
    sin_ϕ = Array{Taylor1{S}}(undef, N, N)
    sin2_ϕ = Array{Taylor1{S}}(undef, N, N)
    sin3_ϕ = Array{Taylor1{S}}(undef, N, N)
    sin4_ϕ = Array{Taylor1{S}}(undef, N, N)
    ϕ = Array{Taylor1{S}}(undef, N, N)
    cos_ϕ = Array{Taylor1{S}}(undef, N, N)
    P_2_sin_ϕ = Array{Taylor1{S}}(undef, N, N)
    ∂P_2_sin_ϕ = Array{Taylor1{S}}(undef, N, N)
    P_3_sin_ϕ = Array{Taylor1{S}}(undef, N, N)
    ∂P_3_sin_ϕ = Array{Taylor1{S}}(undef, N, N)
    m_c_ϕ_∂P_2 = Array{Taylor1{S}}(undef, N, N)
    m_c_ϕ_∂P_3 = Array{Taylor1{S}}(undef, N, N)
    Λ2j_div_r4 = Array{Taylor1{S}}(undef, N, N)
    Λ3j_div_r5 = Array{Taylor1{S}}(undef, N, N)
    F_J_ξ = Array{Taylor1{S}}(undef, N, N)
    F_J_η = Array{Taylor1{S}}(undef, N, N)
    F_J_ζ = Array{Taylor1{S}}(undef, N, N)
    F_J2_ξ = Array{Taylor1{S}}(undef, N, N)
    F_J2_η = Array{Taylor1{S}}(undef, N, N)
    F_J2_ζ = Array{Taylor1{S}}(undef, N, N)
    F_J3_ξ = Array{Taylor1{S}}(undef, N, N)
    F_J3_η = Array{Taylor1{S}}(undef, N, N)
    F_J3_ζ = Array{Taylor1{S}}(undef, N, N)
    ξx = Array{Taylor1{S}}(undef, N, N)
    ξy = Array{Taylor1{S}}(undef, N, N)
    ξz = Array{Taylor1{S}}(undef, N, N)
    ηx = Array{Taylor1{S}}(undef, N, N)
    ηy = Array{Taylor1{S}}(undef, N, N)
    ηz = Array{Taylor1{S}}(undef, N, N)
    ηx1 = Array{Taylor1{S}}(undef, N, N)
    ηy1 = Array{Taylor1{S}}(undef, N, N)
    ηz1 = Array{Taylor1{S}}(undef, N, N)
    ηx2 = Array{Taylor1{S}}(undef, N, N)
    ηy2 = Array{Taylor1{S}}(undef, N, N)
    ηz2 = Array{Taylor1{S}}(undef, N, N)
    ζx = Array{Taylor1{S}}(undef, N, N)
    ζy = Array{Taylor1{S}}(undef, N, N)
    ζz = Array{Taylor1{S}}(undef, N, N)
    ζx1 = Array{Taylor1{S}}(undef, N, N)
    ζy1 = Array{Taylor1{S}}(undef, N, N)
    ζz1 = Array{Taylor1{S}}(undef, N, N)
    ζx2 = Array{Taylor1{S}}(undef, N, N)
    ζy2 = Array{Taylor1{S}}(undef, N, N)
    ζz2 = Array{Taylor1{S}}(undef, N, N)

    # extended-body accelerations
    accX = Array{Taylor1{S}}(undef, N)
    accY = Array{Taylor1{S}}(undef, N)
    accZ = Array{Taylor1{S}}(undef, N)

    # rotations to and from Earth, Sun and Moon pole-oriented frames
    local dsj2k = t-2.451545e6 # J2000.0 = 2.451545e6
    local αs = deg2rad(α_p_sun*one(t))
    local δs = deg2rad(δ_p_sun*one(t))
    local αm = moon_pole_ra(dsj2k)
    local δm = moon_pole_dec(dsj2k)
    local M_ = Array{Taylor1{S}}(undef, 3, 3, N)
    local M_[:,:,ea] = t2c_jpl_de430(dsj2k)
    local M_[:,:,su] = pole_rotation( αs, δs )
    local M_[:,:,mo] = pole_rotation( αm, δm )

    for j in _1_to_N
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

    for j in j2_body_index
        for i in _1_to_N
            if i == j
            else
                t31[i,j] = zero_q_1
                t32[i,j] = zero_q_1
                t33[i,j] = zero_q_1
                F_J2_x[i,j] = zero_q_1
                F_J2_y[i,j] = zero_q_1
                F_J2_z[i,j] = zero_q_1
                F_J2_x1[i,j] = zero_q_1
                F_J2_y1[i,j] = zero_q_1
                F_J2_z1[i,j] = zero_q_1
                F_J2_x2[i,j] = zero_q_1
                F_J2_y2[i,j] = zero_q_1
                F_J2_z2[i,j] = zero_q_1
                sin_ϕ[i,j] = zero_q_1
                sin2_ϕ[i,j] = zero_q_1
                sin3_ϕ[i,j] = zero_q_1
                sin4_ϕ[i,j] = zero_q_1
                ϕ[i,j] = zero_q_1
                cos_ϕ[i,j] = zero_q_1
                P_2_sin_ϕ[i,j] = zero_q_1
                ∂P_2_sin_ϕ[i,j] = zero_q_1
                P_3_sin_ϕ[i,j] = zero_q_1
                ∂P_3_sin_ϕ[i,j] = zero_q_1
                m_c_ϕ_∂P_2[i,j] = zero_q_1
                m_c_ϕ_∂P_3[i,j] = zero_q_1
                Λ2j_div_r4[i,j] = zero_q_1
                Λ3j_div_r5[i,j] = zero_q_1
                F_J_ξ[i,j] = zero_q_1
                F_J_η[i,j] = zero_q_1
                F_J_ζ[i,j] = zero_q_1
                F_J2_ξ[i,j] = zero_q_1
                F_J2_η[i,j] = zero_q_1
                F_J2_ζ[i,j] = zero_q_1
                F_J3_ξ[i,j] = zero_q_1
                F_J3_η[i,j] = zero_q_1
                F_J3_ζ[i,j] = zero_q_1
                ξx[i,j] = zero_q_1
                ξy[i,j] = zero_q_1
                ξz[i,j] = zero_q_1
                ηx[i,j] = zero_q_1
                ηy[i,j] = zero_q_1
                ηz[i,j] = zero_q_1
                ηx1[i,j] = zero_q_1
                ηy1[i,j] = zero_q_1
                ηz1[i,j] = zero_q_1
                ηx2[i,j] = zero_q_1
                ηy2[i,j] = zero_q_1
                ηz2[i,j] = zero_q_1
                ζx[i,j] = zero_q_1
                ζy[i,j] = zero_q_1
                ζz[i,j] = zero_q_1
                ζx1[i,j] = zero_q_1
                ζy1[i,j] = zero_q_1
                ζz1[i,j] = zero_q_1
                ζx2[i,j] = zero_q_1
                ζy2[i,j] = zero_q_1
                ζz2[i,j] = zero_q_1
            end #if i == j
        end #for i in _1_to_N
    end #for j in j2_body_index

    #compute point-mass Newtonian accelerations, all bodies
    for j in _1_to_N
        for i in _1_to_N
            # i == j && continue
            if i == j
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

                #J2 accelerations, if j-th body is flattened
                if UJ_interaction[i,j]
                    # # rotate from inertial frame to extended-body frame
                    t31[i,j] = X[i,j]*M_[1,3,j]
                    t32[i,j] = Y[i,j]*M_[2,3,j]
                    t33[i,j] = Z[i,j]*M_[3,3,j]
                    r_sin_ϕ[i,j] = (t31[i,j]+t32[i,j])+t33[i,j]

                    # compute cartesian coordinates of acceleration due to body figure in body frame
                    sin_ϕ[i,j] = r_sin_ϕ[i,j]/r_p1d2[i,j] # z/r = latitude of point mass wrt body frame
                    ϕ[i,j] = asin(sin_ϕ[i,j])
                    cos_ϕ[i,j] = cos(ϕ[i,j])
                    sin2_ϕ[i,j] = sin_ϕ[i,j]^2
                    sin3_ϕ[i,j] = sin_ϕ[i,j]^3
                    P_2_sin_ϕ[i,j] = 1.5sin2_ϕ[i,j] - 0.5
                    ∂P_2_sin_ϕ[i,j] = 3sin_ϕ[i,j]
                    P_3_sin_ϕ[i,j] = (-1.5sin_ϕ[i,j]) + (2.5sin3_ϕ[i,j])
                    ∂P_3_sin_ϕ[i,j] = -1.5 + 7.5sin2_ϕ[i,j]
                    Λ2j_div_r4[i,j] = (-Λ2[j])/(r_p2[i,j]^2)
                    Λ3j_div_r5[i,j] = (-Λ3[j])/(r_p1d2[i,j]^5)
                    m_c_ϕ_∂P_2[i,j] = (-cos_ϕ[i,j])*∂P_2_sin_ϕ[i,j]
                    m_c_ϕ_∂P_3[i,j] = (-cos_ϕ[i,j])*∂P_3_sin_ϕ[i,j]
                    F_J2_ξ[i,j] = ( Λ2j_div_r4[i,j]*(3P_2_sin_ϕ[i,j]) )
                    #F_J2_η[i,j] = zero_q_1
                    F_J2_ζ[i,j] = Λ2j_div_r4[i,j]*m_c_ϕ_∂P_2[i,j]
                    F_J3_ξ[i,j] = ( Λ3j_div_r5[i,j]*(4P_3_sin_ϕ[i,j]) )
                    #F_J3_η[i,j] = zero_q_1
                    F_J3_ζ[i,j] = Λ3j_div_r5[i,j]*m_c_ϕ_∂P_3[i,j]
                    F_J_ξ[i,j] = F_J2_ξ[i,j] + F_J3_ξ[i,j]
                    #F_J_η[i,j] = zero_q_1
                    F_J_ζ[i,j] = F_J2_ζ[i,j] + F_J3_ζ[i,j]
                    #Compute unit vectors ξ,η,ζ
                    ξx[i,j] = X[i,j]/r_p1d2[i,j]
                    ξy[i,j] = Y[i,j]/r_p1d2[i,j]
                    ξz[i,j] = Z[i,j]/r_p1d2[i,j]
                    #Compute η = p x ξ
                    ηx1[i,j] = M_[2,3,j]*ξz[i,j]
                    ηy1[i,j] = M_[3,3,j]*ξx[i,j]
                    ηz1[i,j] = M_[1,3,j]*ξy[i,j]
                    ηx2[i,j] = M_[3,3,j]*ξy[i,j]
                    ηy2[i,j] = M_[1,3,j]*ξz[i,j]
                    ηz2[i,j] = M_[2,3,j]*ξx[i,j]
                    ηx[i,j] = ηx1[i,j] - ηx2[i,j]
                    ηy[i,j] = ηy1[i,j] - ηy2[i,j]
                    ηz[i,j] = ηz1[i,j] - ηz2[i,j]
                    #Compute ζ = ξ x η
                    ζx1[i,j] = ξy[i,j]*ηz[i,j]
                    ζy1[i,j] = ξz[i,j]*ηx[i,j]
                    ζz1[i,j] = ξx[i,j]*ηy[i,j]
                    ζx2[i,j] = ξz[i,j]*ηy[i,j]
                    ζy2[i,j] = ξx[i,j]*ηz[i,j]
                    ζz2[i,j] = ξy[i,j]*ηx[i,j]
                    ζx[i,j] = ζx1[i,j] - ζx2[i,j]
                    ζy[i,j] = ζy1[i,j] - ζy2[i,j]
                    ζz[i,j] = ζz1[i,j] - ζz2[i,j]
                    # compute cartesian coordinates of acceleration due to body figure in inertial frame
                    F_J2_x1[i,j] = F_J_ξ[i,j]*ξx[i,j]
                    F_J2_y1[i,j] = F_J_ξ[i,j]*ξy[i,j]
                    F_J2_z1[i,j] = F_J_ξ[i,j]*ξz[i,j]
                    F_J2_x2[i,j] = F_J_ζ[i,j]*ζx[i,j]
                    F_J2_y2[i,j] = F_J_ζ[i,j]*ζy[i,j]
                    F_J2_z2[i,j] = F_J_ζ[i,j]*ζz[i,j]
                    F_J2_x[i,j] = F_J2_x1[i,j] + F_J2_x2[i,j]
                    F_J2_y[i,j] = F_J2_y1[i,j] + F_J2_y2[i,j]
                    F_J2_z[i,j] = F_J2_z1[i,j] + F_J2_z2[i,j]
                end # if UJ_interaction[i,j]
            end #if i != j
        end #for, i
        v2[j] = ( (dq[3j-2]^2)+(dq[3j-1]^2) ) + (dq[3j]^2)
    end #for, j

    for j in _1_to_N
        for i in _1_to_N
            # i == j && continue
            if i == j
            else
                newtonianNb_Potential[j] = newtonianNb_Potential[j] + newtonian1b_Potential[i, j]
                newtonX[j] = newtonX[j] + newton_acc_X[i,j]
                newtonY[j] = newtonY[j] + newton_acc_Y[i,j]
                newtonZ[j] = newtonZ[j] + newton_acc_Z[i,j]
                if UJ_interaction[i,j]
                    # # add result to total acceleration on upon j-th body figure due to i-th point mass
                    # @show "acc",j,"+μ",i,"Λ2",j
                    accX[j] = accX[j] + (μ[i]*F_J2_x[i,j])
                    accY[j] = accY[j] + (μ[i]*F_J2_y[i,j])
                    accZ[j] = accZ[j] + (μ[i]*F_J2_z[i,j])
                    # # reaction force on i-th body
                    # @show "acc",i,"-μ",j,"Λ2",j
                    accX[i] = accX[i] - (μ[j]*F_J2_x[i,j])
                    accY[i] = accY[i] - (μ[j]*F_J2_y[i,j])
                    accZ[i] = accZ[i] - (μ[j]*F_J2_z[i,j])
                end
            end
        end
    end

    #post-Newtonian corrections to gravitational acceleration
    #Moyer, 1971, page 7 eq. 35
    # post-Newtonian iterative procedure setup and initialization
    for j in _1_to_N
        for i in _1_to_N
            # i == j && continue
            if i == j
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
                for k in Base.OneTo(postnewton_iter)
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
            end
        end
        postNewtonX[j,1] = newtonX[j]
        postNewtonY[j,1] = newtonY[j]
        postNewtonZ[j,1] = newtonZ[j]
        for k in Base.OneTo(postnewton_iter)
            pntempX[j,k] = zero_q_1
            pntempY[j,k] = zero_q_1
            pntempZ[j,k] = zero_q_1
        end
    end

    # post-Newtonian iterations
    for k in Base.OneTo(postnewton_iter)
        for j in _1_to_N
            for i in _1_to_N
                # i == j && continue
                if i == j
                else
                    pNX_t_X[i,j,k] = postNewtonX[i,k]*X[i,j]
                    pNY_t_Y[i,j,k] = postNewtonY[i,k]*Y[i,j]
                    pNZ_t_Z[i,j,k] = postNewtonZ[i,k]*Z[i,j]
                    pn1[i,j,k] = (  pn1t1_7[i,j]  +  ( (pNX_t_X[i,j,k]+pNY_t_Y[i,j,k]) + pNZ_t_Z[i,j,k] )  )

                    X_t_pn1[i,j,k] = newton_acc_X[i,j]*pn1[i,j,k]
                    Y_t_pn1[i,j,k] = newton_acc_Y[i,j]*pn1[i,j,k]
                    Z_t_pn1[i,j,k] = newton_acc_Z[i,j]*pn1[i,j,k]

                    pNX_t_pn3[i,j,k] = postNewtonX[i,k]*pn3[i,j]
                    pNY_t_pn3[i,j,k] = postNewtonY[i,k]*pn3[i,j]
                    pNZ_t_pn3[i,j,k] = postNewtonZ[i,k]*pn3[i,j]
                end # if i != j
            end #for i
        end #for j
        for j in _1_to_N
            for i in _1_to_N
                # i == j && continue
                if i == j
                else
                    pntempX[j,k] = pntempX[j,k] + ( X_t_pn1[i,j,k] + (U_t_pn2[i,j]+pNX_t_pn3[i,j,k]) )
                    pntempY[j,k] = pntempY[j,k] + ( Y_t_pn1[i,j,k] + (V_t_pn2[i,j]+pNY_t_pn3[i,j,k]) )
                    pntempZ[j,k] = pntempZ[j,k] + ( Z_t_pn1[i,j,k] + (W_t_pn2[i,j]+pNZ_t_pn3[i,j,k]) )
                end
            end
            postNewtonX[j,k+1] = pntempX[j,k]*c_m2
            postNewtonY[j,k+1] = pntempY[j,k]*c_m2
            postNewtonZ[j,k+1] = pntempZ[j,k]*c_m2
        end
    end #for k in Base.OneTo(postnewton_iter) # (post-Newtonian iterations)

    #fill accelerations (post-Newtonian and extended body accelerations)
    for i in _1_to_N
        dq[3(N+i)-2] = postNewtonX[i,postnewton_iter+1] + accX[i]
        dq[3(N+i)-1] = postNewtonY[i,postnewton_iter+1] + accY[i]
        dq[3(N+i)  ] = postNewtonZ[i,postnewton_iter+1] + accZ[i]
    end

    nothing
end

function TaylorIntegration.jetcoeffs!(::Val{NBP_pN_A_J23E_J23M_J2S!}, t::Taylor1{_T}, q::AbstractVector{Taylor1{_S}}, dq::AbstractVector{Taylor1{_S}}, params) where {_T <: Real, _S <: Number}
    order = t.order
    local S = eltype(q[1])
    local N = Int(length(q) / 6)
    local _1_to_N = Base.OneTo(N)
    local postnewton_iter = 1
    local j2_body_index = [su, ea, mo]
    local c_p2 = 29979.063823897606
    local c_m2 = 3.3356611996764786e-5
    local zero_q_1 = zero(q[1])
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
    for j = _1_to_N
        for i = _1_to_N
            if i == j
            else
                pn2[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                pn3[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                _4ϕj[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                ϕi_plus_4ϕj[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                sj2_plus_2si2_minus_4vivj[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                ϕs_and_vs[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                UU[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                VV[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                WW[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                U_t_pn2[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                V_t_pn2[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                W_t_pn2[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                vi_dot_vj[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                newton_acc_X[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                newton_acc_Y[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                newton_acc_Z[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                pn1t1_7[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
            end
        end
    end
    t31 = Array{Taylor1{S}}(undef, N, N)
    t32 = Array{Taylor1{S}}(undef, N, N)
    t33 = Array{Taylor1{S}}(undef, N, N)
    r_sin_ϕ = Array{Taylor1{S}}(undef, N, N)
    F_J2_x = Array{Taylor1{S}}(undef, N, N)
    F_J2_y = Array{Taylor1{S}}(undef, N, N)
    F_J2_z = Array{Taylor1{S}}(undef, N, N)
    F_J2_x1 = Array{Taylor1{S}}(undef, N, N)
    F_J2_y1 = Array{Taylor1{S}}(undef, N, N)
    F_J2_z1 = Array{Taylor1{S}}(undef, N, N)
    F_J2_x2 = Array{Taylor1{S}}(undef, N, N)
    F_J2_y2 = Array{Taylor1{S}}(undef, N, N)
    F_J2_z2 = Array{Taylor1{S}}(undef, N, N)
    sin_ϕ = Array{Taylor1{S}}(undef, N, N)
    sin2_ϕ = Array{Taylor1{S}}(undef, N, N)
    sin3_ϕ = Array{Taylor1{S}}(undef, N, N)
    sin4_ϕ = Array{Taylor1{S}}(undef, N, N)
    ϕ = Array{Taylor1{S}}(undef, N, N)
    cos_ϕ = Array{Taylor1{S}}(undef, N, N)
    P_2_sin_ϕ = Array{Taylor1{S}}(undef, N, N)
    ∂P_2_sin_ϕ = Array{Taylor1{S}}(undef, N, N)
    P_3_sin_ϕ = Array{Taylor1{S}}(undef, N, N)
    ∂P_3_sin_ϕ = Array{Taylor1{S}}(undef, N, N)
    m_c_ϕ_∂P_2 = Array{Taylor1{S}}(undef, N, N)
    m_c_ϕ_∂P_3 = Array{Taylor1{S}}(undef, N, N)
    Λ2j_div_r4 = Array{Taylor1{S}}(undef, N, N)
    Λ3j_div_r5 = Array{Taylor1{S}}(undef, N, N)
    F_J_ξ = Array{Taylor1{S}}(undef, N, N)
    F_J_η = Array{Taylor1{S}}(undef, N, N)
    F_J_ζ = Array{Taylor1{S}}(undef, N, N)
    F_J2_ξ = Array{Taylor1{S}}(undef, N, N)
    F_J2_η = Array{Taylor1{S}}(undef, N, N)
    F_J2_ζ = Array{Taylor1{S}}(undef, N, N)
    F_J3_ξ = Array{Taylor1{S}}(undef, N, N)
    F_J3_η = Array{Taylor1{S}}(undef, N, N)
    F_J3_ζ = Array{Taylor1{S}}(undef, N, N)
    ξx = Array{Taylor1{S}}(undef, N, N)
    ξy = Array{Taylor1{S}}(undef, N, N)
    ξz = Array{Taylor1{S}}(undef, N, N)
    ηx = Array{Taylor1{S}}(undef, N, N)
    ηy = Array{Taylor1{S}}(undef, N, N)
    ηz = Array{Taylor1{S}}(undef, N, N)
    ηx1 = Array{Taylor1{S}}(undef, N, N)
    ηy1 = Array{Taylor1{S}}(undef, N, N)
    ηz1 = Array{Taylor1{S}}(undef, N, N)
    ηx2 = Array{Taylor1{S}}(undef, N, N)
    ηy2 = Array{Taylor1{S}}(undef, N, N)
    ηz2 = Array{Taylor1{S}}(undef, N, N)
    ζx = Array{Taylor1{S}}(undef, N, N)
    ζy = Array{Taylor1{S}}(undef, N, N)
    ζz = Array{Taylor1{S}}(undef, N, N)
    ζx1 = Array{Taylor1{S}}(undef, N, N)
    ζy1 = Array{Taylor1{S}}(undef, N, N)
    ζz1 = Array{Taylor1{S}}(undef, N, N)
    ζx2 = Array{Taylor1{S}}(undef, N, N)
    ζy2 = Array{Taylor1{S}}(undef, N, N)
    ζz2 = Array{Taylor1{S}}(undef, N, N)
    accX = Array{Taylor1{S}}(undef, N)
    accY = Array{Taylor1{S}}(undef, N)
    accZ = Array{Taylor1{S}}(undef, N)
    local dsj2k = t - 2.451545e6
    local αs = deg2rad(α_p_sun * one(t))
    local δs = deg2rad(δ_p_sun * one(t))
    local αm = moon_pole_ra(dsj2k)
    local δm = moon_pole_dec(dsj2k)
    local M_ = Array{Taylor1{S}}(undef, 3, 3, N)
    local M_[:, :, ea] = t2c_jpl_de430(dsj2k)
    local M_[:, :, su] = pole_rotation(αs, δs)
    local M_[:, :, mo] = pole_rotation(αm, δm)
    for j = _1_to_N
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
    for j = j2_body_index
        for i = _1_to_N
            if i == j
            else
                t31[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                t32[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                t33[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                F_J2_x[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                F_J2_y[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                F_J2_z[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                F_J2_x1[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                F_J2_y1[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                F_J2_z1[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                F_J2_x2[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                F_J2_y2[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                F_J2_z2[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                sin_ϕ[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                sin2_ϕ[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                sin3_ϕ[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                sin4_ϕ[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                ϕ[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                cos_ϕ[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                P_2_sin_ϕ[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                ∂P_2_sin_ϕ[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                P_3_sin_ϕ[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                ∂P_3_sin_ϕ[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                m_c_ϕ_∂P_2[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                m_c_ϕ_∂P_3[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                Λ2j_div_r4[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                Λ3j_div_r5[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                F_J_ξ[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                F_J_η[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                F_J_ζ[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                F_J2_ξ[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                F_J2_η[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                F_J2_ζ[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                F_J3_ξ[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                F_J3_η[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                F_J3_ζ[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                ξx[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                ξy[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                ξz[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                ηx[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                ηy[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                ηz[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                ηx1[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                ηy1[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                ηz1[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                ηx2[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                ηy2[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                ηz2[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                ζx[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                ζy[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                ζz[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                ζx1[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                ζy1[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                ζz1[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                ζx2[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                ζy2[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                ζz2[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
            end
        end
    end
    tmp849 = Array{Taylor1{_S}}(undef, size(dq))
    tmp849 .= Taylor1(zero(_S), order)
    tmp851 = Array{Taylor1{_S}}(undef, size(dq))
    tmp851 .= Taylor1(zero(_S), order)
    tmp854 = Array{Taylor1{_S}}(undef, size(dq))
    tmp854 .= Taylor1(zero(_S), order)
    tmp856 = Array{Taylor1{_S}}(undef, size(dq))
    tmp856 .= Taylor1(zero(_S), order)
    tmp859 = Array{Taylor1{_S}}(undef, size(dq))
    tmp859 .= Taylor1(zero(_S), order)
    tmp861 = Array{Taylor1{_S}}(undef, size(dq))
    tmp861 .= Taylor1(zero(_S), order)
    pn2x = Array{Taylor1{_S}}(undef, size(X))
    pn2x .= Taylor1(zero(_S), order)
    pn2y = Array{Taylor1{_S}}(undef, size(Y))
    pn2y .= Taylor1(zero(_S), order)
    pn2z = Array{Taylor1{_S}}(undef, size(Z))
    pn2z .= Taylor1(zero(_S), order)
    tmp869 = Array{Taylor1{_S}}(undef, size(UU))
    tmp869 .= Taylor1(zero(_S), order)
    tmp872 = Array{Taylor1{_S}}(undef, size(X))
    tmp872 .= Taylor1(zero(_S), order)
    tmp874 = Array{Taylor1{_S}}(undef, size(Y))
    tmp874 .= Taylor1(zero(_S), order)
    tmp875 = Array{Taylor1{_S}}(undef, size(tmp872))
    tmp875 .= Taylor1(zero(_S), order)
    tmp877 = Array{Taylor1{_S}}(undef, size(Z))
    tmp877 .= Taylor1(zero(_S), order)
    tmp885 = Array{Taylor1{_S}}(undef, size(pn2x))
    tmp885 .= Taylor1(zero(_S), order)
    tmp886 = Array{Taylor1{_S}}(undef, size(tmp885))
    tmp886 .= Taylor1(zero(_S), order)
    tmp900 = Array{Taylor1{_S}}(undef, size(t31))
    tmp900 .= Taylor1(zero(_S), order)
    tmp1049 = Array{Taylor1{_S}}(undef, size(sin_ϕ))
    tmp1049 .= Taylor1(zero(_S), order)
    tmp1050 = Array{Taylor1{_S}}(undef, size(ϕ))
    tmp1050 .= Taylor1(zero(_S), order)
    tmp910 = Array{Taylor1{_S}}(undef, size(sin2_ϕ))
    tmp910 .= Taylor1(zero(_S), order)
    tmp916 = Array{Taylor1{_S}}(undef, size(sin_ϕ))
    tmp916 .= Taylor1(zero(_S), order)
    tmp918 = Array{Taylor1{_S}}(undef, size(sin3_ϕ))
    tmp918 .= Taylor1(zero(_S), order)
    tmp922 = Array{Taylor1{_S}}(undef, size(sin2_ϕ))
    tmp922 .= Taylor1(zero(_S), order)
    tmp924 = Array{Taylor1{_S}}(undef, size(Λ2))
    tmp924 .= Taylor1(zero(_S), order)
    tmp926 = Array{Taylor1{_S}}(undef, size(r_p2))
    tmp926 .= Taylor1(zero(_S), order)
    tmp928 = Array{Taylor1{_S}}(undef, size(Λ3))
    tmp928 .= Taylor1(zero(_S), order)
    tmp930 = Array{Taylor1{_S}}(undef, size(r_p1d2))
    tmp930 .= Taylor1(zero(_S), order)
    tmp932 = Array{Taylor1{_S}}(undef, size(cos_ϕ))
    tmp932 .= Taylor1(zero(_S), order)
    tmp934 = Array{Taylor1{_S}}(undef, size(cos_ϕ))
    tmp934 .= Taylor1(zero(_S), order)
    tmp937 = Array{Taylor1{_S}}(undef, size(Λ2j_div_r4))
    tmp937 .= Taylor1(zero(_S), order)
    tmp941 = Array{Taylor1{_S}}(undef, size(Λ3j_div_r5))
    tmp941 .= Taylor1(zero(_S), order)
    tmp977 = Array{Taylor1{_S}}(undef, size(dq))
    tmp977 .= Taylor1(zero(_S), order)
    tmp979 = Array{Taylor1{_S}}(undef, size(dq))
    tmp979 .= Taylor1(zero(_S), order)
    tmp980 = Array{Taylor1{_S}}(undef, size(tmp977))
    tmp980 .= Taylor1(zero(_S), order)
    tmp982 = Array{Taylor1{_S}}(undef, size(dq))
    tmp982 .= Taylor1(zero(_S), order)
    for j = _1_to_N
        for i = _1_to_N
            if i == j
            else
                X[i, j] = Taylor1(constant_term(q[3i - 2]) - constant_term(q[3j - 2]), order)
                Y[i, j] = Taylor1(constant_term(q[3i - 1]) - constant_term(q[3j - 1]), order)
                Z[i, j] = Taylor1(constant_term(q[3i]) - constant_term(q[3j]), order)
                U[i, j] = Taylor1(constant_term(dq[3i - 2]) - constant_term(dq[3j - 2]), order)
                V[i, j] = Taylor1(constant_term(dq[3i - 1]) - constant_term(dq[3j - 1]), order)
                W[i, j] = Taylor1(constant_term(dq[3i]) - constant_term(dq[3j]), order)
                tmp849[3j - 2] = Taylor1(constant_term(4) * constant_term(dq[3j - 2]), order)
                tmp851[3i - 2] = Taylor1(constant_term(3) * constant_term(dq[3i - 2]), order)
                _4U_m_3X[i, j] = Taylor1(constant_term(tmp849[3j - 2]) - constant_term(tmp851[3i - 2]), order)
                tmp854[3j - 1] = Taylor1(constant_term(4) * constant_term(dq[3j - 1]), order)
                tmp856[3i - 1] = Taylor1(constant_term(3) * constant_term(dq[3i - 1]), order)
                _4V_m_3Y[i, j] = Taylor1(constant_term(tmp854[3j - 1]) - constant_term(tmp856[3i - 1]), order)
                tmp859[3j] = Taylor1(constant_term(4) * constant_term(dq[3j]), order)
                tmp861[3i] = Taylor1(constant_term(3) * constant_term(dq[3i]), order)
                _4W_m_3Z[i, j] = Taylor1(constant_term(tmp859[3j]) - constant_term(tmp861[3i]), order)
                pn2x[i, j] = Taylor1(constant_term(X[i, j]) * constant_term(_4U_m_3X[i, j]), order)
                pn2y[i, j] = Taylor1(constant_term(Y[i, j]) * constant_term(_4V_m_3Y[i, j]), order)
                pn2z[i, j] = Taylor1(constant_term(Z[i, j]) * constant_term(_4W_m_3Z[i, j]), order)
                UU[i, j] = Taylor1(constant_term(dq[3i - 2]) * constant_term(dq[3j - 2]), order)
                VV[i, j] = Taylor1(constant_term(dq[3i - 1]) * constant_term(dq[3j - 1]), order)
                WW[i, j] = Taylor1(constant_term(dq[3i]) * constant_term(dq[3j]), order)
                tmp869[i, j] = Taylor1(constant_term(UU[i, j]) + constant_term(VV[i, j]), order)
                vi_dot_vj[i, j] = Taylor1(constant_term(tmp869[i, j]) + constant_term(WW[i, j]), order)
                tmp872[i, j] = Taylor1(constant_term(X[i, j]) ^ constant_term(2), order)
                tmp874[i, j] = Taylor1(constant_term(Y[i, j]) ^ constant_term(2), order)
                tmp875[i, j] = Taylor1(constant_term(tmp872[i, j]) + constant_term(tmp874[i, j]), order)
                tmp877[i, j] = Taylor1(constant_term(Z[i, j]) ^ constant_term(2), order)
                r_p2[i, j] = Taylor1(constant_term(tmp875[i, j]) + constant_term(tmp877[i, j]), order)
                r_p1d2[i, j] = Taylor1(sqrt(constant_term(r_p2[i, j])), order)
                r_p3d2[i, j] = Taylor1(constant_term(r_p2[i, j]) ^ constant_term(1.5), order)
                r_p7d2[i, j] = Taylor1(constant_term(r_p2[i, j]) ^ constant_term(3.5), order)
                newtonianCoeff[i, j] = Taylor1(constant_term(μ[i]) / constant_term(r_p3d2[i, j]), order)
                tmp885[i, j] = Taylor1(constant_term(pn2x[i, j]) + constant_term(pn2y[i, j]), order)
                tmp886[i, j] = Taylor1(constant_term(tmp885[i, j]) + constant_term(pn2z[i, j]), order)
                pn2[i, j] = Taylor1(constant_term(newtonianCoeff[i, j]) * constant_term(tmp886[i, j]), order)
                newton_acc_X[i, j] = Taylor1(constant_term(X[i, j]) * constant_term(newtonianCoeff[i, j]), order)
                newton_acc_Y[i, j] = Taylor1(constant_term(Y[i, j]) * constant_term(newtonianCoeff[i, j]), order)
                newton_acc_Z[i, j] = Taylor1(constant_term(Z[i, j]) * constant_term(newtonianCoeff[i, j]), order)
                newtonian1b_Potential[i, j] = Taylor1(constant_term(μ[i]) / constant_term(r_p1d2[i, j]), order)
                pn3[i, j] = Taylor1(constant_term(3.5) * constant_term(newtonian1b_Potential[i, j]), order)
                U_t_pn2[i, j] = Taylor1(constant_term(pn2[i, j]) * constant_term(U[i, j]), order)
                V_t_pn2[i, j] = Taylor1(constant_term(pn2[i, j]) * constant_term(V[i, j]), order)
                W_t_pn2[i, j] = Taylor1(constant_term(pn2[i, j]) * constant_term(W[i, j]), order)
                if UJ_interaction[i, j]
                    t31[i, j] = Taylor1(constant_term(X[i, j]) * constant_term(M_[1, 3, j]), order)
                    t32[i, j] = Taylor1(constant_term(Y[i, j]) * constant_term(M_[2, 3, j]), order)
                    t33[i, j] = Taylor1(constant_term(Z[i, j]) * constant_term(M_[3, 3, j]), order)
                    tmp900[i, j] = Taylor1(constant_term(t31[i, j]) + constant_term(t32[i, j]), order)
                    r_sin_ϕ[i, j] = Taylor1(constant_term(tmp900[i, j]) + constant_term(t33[i, j]), order)
                    sin_ϕ[i, j] = Taylor1(constant_term(r_sin_ϕ[i, j]) / constant_term(r_p1d2[i, j]), order)
                    ϕ[i, j] = Taylor1(asin(constant_term(sin_ϕ[i, j])), order)
                    tmp1049[i, j] = Taylor1(sqrt(1 - constant_term(sin_ϕ[i, j]) ^ 2), order)
                    cos_ϕ[i, j] = Taylor1(cos(constant_term(ϕ[i, j])), order)
                    tmp1050[i, j] = Taylor1(sin(constant_term(ϕ[i, j])), order)
                    sin2_ϕ[i, j] = Taylor1(constant_term(sin_ϕ[i, j]) ^ constant_term(2), order)
                    sin3_ϕ[i, j] = Taylor1(constant_term(sin_ϕ[i, j]) ^ constant_term(3), order)
                    tmp910[i, j] = Taylor1(constant_term(1.5) * constant_term(sin2_ϕ[i, j]), order)
                    P_2_sin_ϕ[i, j] = Taylor1(constant_term(tmp910[i, j]) - constant_term(0.5), order)
                    ∂P_2_sin_ϕ[i, j] = Taylor1(constant_term(3) * constant_term(sin_ϕ[i, j]), order)
                    tmp916[i, j] = Taylor1(constant_term(-1.5) * constant_term(sin_ϕ[i, j]), order)
                    tmp918[i, j] = Taylor1(constant_term(2.5) * constant_term(sin3_ϕ[i, j]), order)
                    P_3_sin_ϕ[i, j] = Taylor1(constant_term(tmp916[i, j]) + constant_term(tmp918[i, j]), order)
                    tmp922[i, j] = Taylor1(constant_term(7.5) * constant_term(sin2_ϕ[i, j]), order)
                    ∂P_3_sin_ϕ[i, j] = Taylor1(constant_term(-1.5) + constant_term(tmp922[i, j]), order)
                    tmp924[j] = Taylor1(-(constant_term(Λ2[j])), order)
                    tmp926[i, j] = Taylor1(constant_term(r_p2[i, j]) ^ constant_term(2), order)
                    Λ2j_div_r4[i, j] = Taylor1(constant_term(tmp924[j]) / constant_term(tmp926[i, j]), order)
                    tmp928[j] = Taylor1(-(constant_term(Λ3[j])), order)
                    tmp930[i, j] = Taylor1(constant_term(r_p1d2[i, j]) ^ constant_term(5), order)
                    Λ3j_div_r5[i, j] = Taylor1(constant_term(tmp928[j]) / constant_term(tmp930[i, j]), order)
                    tmp932[i, j] = Taylor1(-(constant_term(cos_ϕ[i, j])), order)
                    m_c_ϕ_∂P_2[i, j] = Taylor1(constant_term(tmp932[i, j]) * constant_term(∂P_2_sin_ϕ[i, j]), order)
                    tmp934[i, j] = Taylor1(-(constant_term(cos_ϕ[i, j])), order)
                    m_c_ϕ_∂P_3[i, j] = Taylor1(constant_term(tmp934[i, j]) * constant_term(∂P_3_sin_ϕ[i, j]), order)
                    tmp937[i, j] = Taylor1(constant_term(Λ2j_div_r4[i, j]) * constant_term(3), order)
                    F_J2_ξ[i, j] = Taylor1(constant_term(tmp937[i, j]) * constant_term(P_2_sin_ϕ[i, j]), order)
                    F_J2_ζ[i, j] = Taylor1(constant_term(Λ2j_div_r4[i, j]) * constant_term(m_c_ϕ_∂P_2[i, j]), order)
                    tmp941[i, j] = Taylor1(constant_term(Λ3j_div_r5[i, j]) * constant_term(4), order)
                    F_J3_ξ[i, j] = Taylor1(constant_term(tmp941[i, j]) * constant_term(P_3_sin_ϕ[i, j]), order)
                    F_J3_ζ[i, j] = Taylor1(constant_term(Λ3j_div_r5[i, j]) * constant_term(m_c_ϕ_∂P_3[i, j]), order)
                    F_J_ξ[i, j] = Taylor1(constant_term(F_J2_ξ[i, j]) + constant_term(F_J3_ξ[i, j]), order)
                    F_J_ζ[i, j] = Taylor1(constant_term(F_J2_ζ[i, j]) + constant_term(F_J3_ζ[i, j]), order)
                    ξx[i, j] = Taylor1(constant_term(X[i, j]) / constant_term(r_p1d2[i, j]), order)
                    ξy[i, j] = Taylor1(constant_term(Y[i, j]) / constant_term(r_p1d2[i, j]), order)
                    ξz[i, j] = Taylor1(constant_term(Z[i, j]) / constant_term(r_p1d2[i, j]), order)
                    ηx1[i, j] = Taylor1(constant_term(M_[2, 3, j]) * constant_term(ξz[i, j]), order)
                    ηy1[i, j] = Taylor1(constant_term(M_[3, 3, j]) * constant_term(ξx[i, j]), order)
                    ηz1[i, j] = Taylor1(constant_term(M_[1, 3, j]) * constant_term(ξy[i, j]), order)
                    ηx2[i, j] = Taylor1(constant_term(M_[3, 3, j]) * constant_term(ξy[i, j]), order)
                    ηy2[i, j] = Taylor1(constant_term(M_[1, 3, j]) * constant_term(ξz[i, j]), order)
                    ηz2[i, j] = Taylor1(constant_term(M_[2, 3, j]) * constant_term(ξx[i, j]), order)
                    ηx[i, j] = Taylor1(constant_term(ηx1[i, j]) - constant_term(ηx2[i, j]), order)
                    ηy[i, j] = Taylor1(constant_term(ηy1[i, j]) - constant_term(ηy2[i, j]), order)
                    ηz[i, j] = Taylor1(constant_term(ηz1[i, j]) - constant_term(ηz2[i, j]), order)
                    ζx1[i, j] = Taylor1(constant_term(ξy[i, j]) * constant_term(ηz[i, j]), order)
                    ζy1[i, j] = Taylor1(constant_term(ξz[i, j]) * constant_term(ηx[i, j]), order)
                    ζz1[i, j] = Taylor1(constant_term(ξx[i, j]) * constant_term(ηy[i, j]), order)
                    ζx2[i, j] = Taylor1(constant_term(ξz[i, j]) * constant_term(ηy[i, j]), order)
                    ζy2[i, j] = Taylor1(constant_term(ξx[i, j]) * constant_term(ηz[i, j]), order)
                    ζz2[i, j] = Taylor1(constant_term(ξy[i, j]) * constant_term(ηx[i, j]), order)
                    ζx[i, j] = Taylor1(constant_term(ζx1[i, j]) - constant_term(ζx2[i, j]), order)
                    ζy[i, j] = Taylor1(constant_term(ζy1[i, j]) - constant_term(ζy2[i, j]), order)
                    ζz[i, j] = Taylor1(constant_term(ζz1[i, j]) - constant_term(ζz2[i, j]), order)
                    F_J2_x1[i, j] = Taylor1(constant_term(F_J_ξ[i, j]) * constant_term(ξx[i, j]), order)
                    F_J2_y1[i, j] = Taylor1(constant_term(F_J_ξ[i, j]) * constant_term(ξy[i, j]), order)
                    F_J2_z1[i, j] = Taylor1(constant_term(F_J_ξ[i, j]) * constant_term(ξz[i, j]), order)
                    F_J2_x2[i, j] = Taylor1(constant_term(F_J_ζ[i, j]) * constant_term(ζx[i, j]), order)
                    F_J2_y2[i, j] = Taylor1(constant_term(F_J_ζ[i, j]) * constant_term(ζy[i, j]), order)
                    F_J2_z2[i, j] = Taylor1(constant_term(F_J_ζ[i, j]) * constant_term(ζz[i, j]), order)
                    F_J2_x[i, j] = Taylor1(constant_term(F_J2_x1[i, j]) + constant_term(F_J2_x2[i, j]), order)
                    F_J2_y[i, j] = Taylor1(constant_term(F_J2_y1[i, j]) + constant_term(F_J2_y2[i, j]), order)
                    F_J2_z[i, j] = Taylor1(constant_term(F_J2_z1[i, j]) + constant_term(F_J2_z2[i, j]), order)
                end
            end
        end
        tmp977[3j - 2] = Taylor1(constant_term(dq[3j - 2]) ^ constant_term(2), order)
        tmp979[3j - 1] = Taylor1(constant_term(dq[3j - 1]) ^ constant_term(2), order)
        tmp980[3j - 2] = Taylor1(constant_term(tmp977[3j - 2]) + constant_term(tmp979[3j - 1]), order)
        tmp982[3j] = Taylor1(constant_term(dq[3j]) ^ constant_term(2), order)
        v2[j] = Taylor1(constant_term(tmp980[3j - 2]) + constant_term(tmp982[3j]), order)
    end
    for j = _1_to_N
        for i = _1_to_N
            if i == j
            else
                newtonianNb_Potential[j] = Taylor1(constant_term(newtonianNb_Potential[j]) + constant_term(newtonian1b_Potential[i, j]), order)
                newtonX[j] = Taylor1(constant_term(newtonX[j]) + constant_term(newton_acc_X[i, j]), order)
                newtonY[j] = Taylor1(constant_term(newtonY[j]) + constant_term(newton_acc_Y[i, j]), order)
                newtonZ[j] = Taylor1(constant_term(newtonZ[j]) + constant_term(newton_acc_Z[i, j]), order)
                if UJ_interaction[i, j]
                    accX[j] = Taylor1(constant_term(accX[j]) + constant_term(μ[i] * F_J2_x[i, j]), order)
                    accY[j] = Taylor1(constant_term(accY[j]) + constant_term(μ[i] * F_J2_y[i, j]), order)
                    accZ[j] = Taylor1(constant_term(accZ[j]) + constant_term(μ[i] * F_J2_z[i, j]), order)
                    accX[i] = Taylor1(constant_term(accX[i]) - constant_term(μ[j] * F_J2_x[i, j]), order)
                    accY[i] = Taylor1(constant_term(accY[i]) - constant_term(μ[j] * F_J2_y[i, j]), order)
                    accZ[i] = Taylor1(constant_term(accZ[i]) - constant_term(μ[j] * F_J2_z[i, j]), order)
                end
            end
        end
    end
    tmp1004 = Array{Taylor1{_S}}(undef, size(v2))
    tmp1004 .= Taylor1(zero(_S), order)
    tmp1005 = Array{Taylor1{_S}}(undef, size(v2))
    tmp1005 .= Taylor1(zero(_S), order)
    tmp1007 = Array{Taylor1{_S}}(undef, size(vi_dot_vj))
    tmp1007 .= Taylor1(zero(_S), order)
    Xij_t_Ui = Array{Taylor1{_S}}(undef, size(X))
    Xij_t_Ui .= Taylor1(zero(_S), order)
    Yij_t_Vi = Array{Taylor1{_S}}(undef, size(Y))
    Yij_t_Vi .= Taylor1(zero(_S), order)
    Zij_t_Wi = Array{Taylor1{_S}}(undef, size(Z))
    Zij_t_Wi .= Taylor1(zero(_S), order)
    tmp1013 = Array{Taylor1{_S}}(undef, size(Xij_t_Ui))
    tmp1013 .= Taylor1(zero(_S), order)
    Rij_dot_Vi = Array{Taylor1{_S}}(undef, size(tmp1013))
    Rij_dot_Vi .= Taylor1(zero(_S), order)
    tmp1016 = Array{Taylor1{_S}}(undef, size(Rij_dot_Vi))
    tmp1016 .= Taylor1(zero(_S), order)
    pn1t7 = Array{Taylor1{_S}}(undef, size(tmp1016))
    pn1t7 .= Taylor1(zero(_S), order)
    tmp1019 = Array{Taylor1{_S}}(undef, size(pn1t7))
    tmp1019 .= Taylor1(zero(_S), order)
    pn1t2_7 = Array{Taylor1{_S}}(undef, size(ϕs_and_vs))
    pn1t2_7 .= Taylor1(zero(_S), order)
    for j = _1_to_N
        for i = _1_to_N
            if i == j
            else
                _4ϕj[i, j] = Taylor1(constant_term(4) * constant_term(newtonianNb_Potential[j]), order)
                ϕi_plus_4ϕj[i, j] = Taylor1(constant_term(newtonianNb_Potential[i]) + constant_term(_4ϕj[i, j]), order)
                tmp1004[i] = Taylor1(constant_term(2) * constant_term(v2[i]), order)
                tmp1005[j] = Taylor1(constant_term(v2[j]) + constant_term(tmp1004[i]), order)
                tmp1007[i, j] = Taylor1(constant_term(4) * constant_term(vi_dot_vj[i, j]), order)
                sj2_plus_2si2_minus_4vivj[i, j] = Taylor1(constant_term(tmp1005[j]) - constant_term(tmp1007[i, j]), order)
                ϕs_and_vs[i, j] = Taylor1(constant_term(sj2_plus_2si2_minus_4vivj[i, j]) - constant_term(ϕi_plus_4ϕj[i, j]), order)
                Xij_t_Ui[i, j] = Taylor1(constant_term(X[i, j]) * constant_term(dq[3i - 2]), order)
                Yij_t_Vi[i, j] = Taylor1(constant_term(Y[i, j]) * constant_term(dq[3i - 1]), order)
                Zij_t_Wi[i, j] = Taylor1(constant_term(Z[i, j]) * constant_term(dq[3i]), order)
                tmp1013[i, j] = Taylor1(constant_term(Xij_t_Ui[i, j]) + constant_term(Yij_t_Vi[i, j]), order)
                Rij_dot_Vi[i, j] = Taylor1(constant_term(tmp1013[i, j]) + constant_term(Zij_t_Wi[i, j]), order)
                tmp1016[i, j] = Taylor1(constant_term(Rij_dot_Vi[i, j]) ^ constant_term(2), order)
                pn1t7[i, j] = Taylor1(constant_term(tmp1016[i, j]) / constant_term(r_p2[i, j]), order)
                tmp1019[i, j] = Taylor1(constant_term(1.5) * constant_term(pn1t7[i, j]), order)
                pn1t2_7[i, j] = Taylor1(constant_term(ϕs_and_vs[i, j]) - constant_term(tmp1019[i, j]), order)
                pn1t1_7[i, j] = Taylor1(constant_term(c_p2) + constant_term(pn1t2_7[i, j]), order)
                for k = Base.OneTo(postnewton_iter)
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
        for k = Base.OneTo(postnewton_iter)
            pntempX[j, k] = Taylor1(identity(constant_term(zero_q_1)), order)
            pntempY[j, k] = Taylor1(identity(constant_term(zero_q_1)), order)
            pntempZ[j, k] = Taylor1(identity(constant_term(zero_q_1)), order)
        end
    end
    tmp1025 = Array{Taylor1{_S}}(undef, size(pNX_t_X))
    tmp1025 .= Taylor1(zero(_S), order)
    tmp1026 = Array{Taylor1{_S}}(undef, size(tmp1025))
    tmp1026 .= Taylor1(zero(_S), order)
    for k = Base.OneTo(postnewton_iter)
        for j = _1_to_N
            for i = _1_to_N
                if i == j
                else
                    pNX_t_X[i, j, k] = Taylor1(constant_term(postNewtonX[i, k]) * constant_term(X[i, j]), order)
                    pNY_t_Y[i, j, k] = Taylor1(constant_term(postNewtonY[i, k]) * constant_term(Y[i, j]), order)
                    pNZ_t_Z[i, j, k] = Taylor1(constant_term(postNewtonZ[i, k]) * constant_term(Z[i, j]), order)
                    tmp1025[i, j, k] = Taylor1(constant_term(pNX_t_X[i, j, k]) + constant_term(pNY_t_Y[i, j, k]), order)
                    tmp1026[i, j, k] = Taylor1(constant_term(tmp1025[i, j, k]) + constant_term(pNZ_t_Z[i, j, k]), order)
                    pn1[i, j, k] = Taylor1(constant_term(pn1t1_7[i, j]) + constant_term(tmp1026[i, j, k]), order)
                    X_t_pn1[i, j, k] = Taylor1(constant_term(newton_acc_X[i, j]) * constant_term(pn1[i, j, k]), order)
                    Y_t_pn1[i, j, k] = Taylor1(constant_term(newton_acc_Y[i, j]) * constant_term(pn1[i, j, k]), order)
                    Z_t_pn1[i, j, k] = Taylor1(constant_term(newton_acc_Z[i, j]) * constant_term(pn1[i, j, k]), order)
                    pNX_t_pn3[i, j, k] = Taylor1(constant_term(postNewtonX[i, k]) * constant_term(pn3[i, j]), order)
                    pNY_t_pn3[i, j, k] = Taylor1(constant_term(postNewtonY[i, k]) * constant_term(pn3[i, j]), order)
                    pNZ_t_pn3[i, j, k] = Taylor1(constant_term(postNewtonZ[i, k]) * constant_term(pn3[i, j]), order)
                end
            end
        end
        for j = _1_to_N
            for i = _1_to_N
                if i == j
                else
                    pntempX[j, k] = Taylor1(constant_term(pntempX[j, k]) + constant_term(X_t_pn1[i, j, k] + (U_t_pn2[i, j] + pNX_t_pn3[i, j, k])), order)
                    pntempY[j, k] = Taylor1(constant_term(pntempY[j, k]) + constant_term(Y_t_pn1[i, j, k] + (V_t_pn2[i, j] + pNY_t_pn3[i, j, k])), order)
                    pntempZ[j, k] = Taylor1(constant_term(pntempZ[j, k]) + constant_term(Z_t_pn1[i, j, k] + (W_t_pn2[i, j] + pNZ_t_pn3[i, j, k])), order)
                end
            end
            postNewtonX[j, k + 1] = Taylor1(constant_term(pntempX[j, k]) * constant_term(c_m2), order)
            postNewtonY[j, k + 1] = Taylor1(constant_term(pntempY[j, k]) * constant_term(c_m2), order)
            postNewtonZ[j, k + 1] = Taylor1(constant_term(pntempZ[j, k]) * constant_term(c_m2), order)
        end
    end
    for i = _1_to_N
        dq[3 * (N + i) - 2] = Taylor1(constant_term(postNewtonX[i, postnewton_iter + 1]) + constant_term(accX[i]), order)
        dq[3 * (N + i) - 1] = Taylor1(constant_term(postNewtonY[i, postnewton_iter + 1]) + constant_term(accY[i]), order)
        dq[3 * (N + i)] = Taylor1(constant_term(postNewtonZ[i, postnewton_iter + 1]) + constant_term(accZ[i]), order)
    end
    for __idx = eachindex(q)
        (q[__idx]).coeffs[2] = (dq[__idx]).coeffs[1]
    end
    for ord = 1:order - 1
        ordnext = ord + 1
        for j = _1_to_N
            for i = _1_to_N
                if i == j
                else
                    TaylorSeries.identity!(pn2[i, j], zero_q_1, ord)
                    TaylorSeries.identity!(pn3[i, j], zero_q_1, ord)
                    TaylorSeries.identity!(_4ϕj[i, j], zero_q_1, ord)
                    TaylorSeries.identity!(ϕi_plus_4ϕj[i, j], zero_q_1, ord)
                    TaylorSeries.identity!(sj2_plus_2si2_minus_4vivj[i, j], zero_q_1, ord)
                    TaylorSeries.identity!(ϕs_and_vs[i, j], zero_q_1, ord)
                    TaylorSeries.identity!(UU[i, j], zero_q_1, ord)
                    TaylorSeries.identity!(VV[i, j], zero_q_1, ord)
                    TaylorSeries.identity!(WW[i, j], zero_q_1, ord)
                    TaylorSeries.identity!(U_t_pn2[i, j], zero_q_1, ord)
                    TaylorSeries.identity!(V_t_pn2[i, j], zero_q_1, ord)
                    TaylorSeries.identity!(W_t_pn2[i, j], zero_q_1, ord)
                    TaylorSeries.identity!(vi_dot_vj[i, j], zero_q_1, ord)
                    TaylorSeries.identity!(newton_acc_X[i, j], zero_q_1, ord)
                    TaylorSeries.identity!(newton_acc_Y[i, j], zero_q_1, ord)
                    TaylorSeries.identity!(newton_acc_Z[i, j], zero_q_1, ord)
                    TaylorSeries.identity!(pn1t1_7[i, j], zero_q_1, ord)
                end
            end
        end
        for j = _1_to_N
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
        for j = j2_body_index
            for i = _1_to_N
                if i == j
                else
                    TaylorSeries.identity!(t31[i, j], zero_q_1, ord)
                    TaylorSeries.identity!(t32[i, j], zero_q_1, ord)
                    TaylorSeries.identity!(t33[i, j], zero_q_1, ord)
                    TaylorSeries.identity!(F_J2_x[i, j], zero_q_1, ord)
                    TaylorSeries.identity!(F_J2_y[i, j], zero_q_1, ord)
                    TaylorSeries.identity!(F_J2_z[i, j], zero_q_1, ord)
                    TaylorSeries.identity!(F_J2_x1[i, j], zero_q_1, ord)
                    TaylorSeries.identity!(F_J2_y1[i, j], zero_q_1, ord)
                    TaylorSeries.identity!(F_J2_z1[i, j], zero_q_1, ord)
                    TaylorSeries.identity!(F_J2_x2[i, j], zero_q_1, ord)
                    TaylorSeries.identity!(F_J2_y2[i, j], zero_q_1, ord)
                    TaylorSeries.identity!(F_J2_z2[i, j], zero_q_1, ord)
                    TaylorSeries.identity!(sin_ϕ[i, j], zero_q_1, ord)
                    TaylorSeries.identity!(sin2_ϕ[i, j], zero_q_1, ord)
                    TaylorSeries.identity!(sin3_ϕ[i, j], zero_q_1, ord)
                    TaylorSeries.identity!(sin4_ϕ[i, j], zero_q_1, ord)
                    TaylorSeries.identity!(ϕ[i, j], zero_q_1, ord)
                    TaylorSeries.identity!(cos_ϕ[i, j], zero_q_1, ord)
                    TaylorSeries.identity!(P_2_sin_ϕ[i, j], zero_q_1, ord)
                    TaylorSeries.identity!(∂P_2_sin_ϕ[i, j], zero_q_1, ord)
                    TaylorSeries.identity!(P_3_sin_ϕ[i, j], zero_q_1, ord)
                    TaylorSeries.identity!(∂P_3_sin_ϕ[i, j], zero_q_1, ord)
                    TaylorSeries.identity!(m_c_ϕ_∂P_2[i, j], zero_q_1, ord)
                    TaylorSeries.identity!(m_c_ϕ_∂P_3[i, j], zero_q_1, ord)
                    TaylorSeries.identity!(Λ2j_div_r4[i, j], zero_q_1, ord)
                    TaylorSeries.identity!(Λ3j_div_r5[i, j], zero_q_1, ord)
                    TaylorSeries.identity!(F_J_ξ[i, j], zero_q_1, ord)
                    TaylorSeries.identity!(F_J_η[i, j], zero_q_1, ord)
                    TaylorSeries.identity!(F_J_ζ[i, j], zero_q_1, ord)
                    TaylorSeries.identity!(F_J2_ξ[i, j], zero_q_1, ord)
                    TaylorSeries.identity!(F_J2_η[i, j], zero_q_1, ord)
                    TaylorSeries.identity!(F_J2_ζ[i, j], zero_q_1, ord)
                    TaylorSeries.identity!(F_J3_ξ[i, j], zero_q_1, ord)
                    TaylorSeries.identity!(F_J3_η[i, j], zero_q_1, ord)
                    TaylorSeries.identity!(F_J3_ζ[i, j], zero_q_1, ord)
                    TaylorSeries.identity!(ξx[i, j], zero_q_1, ord)
                    TaylorSeries.identity!(ξy[i, j], zero_q_1, ord)
                    TaylorSeries.identity!(ξz[i, j], zero_q_1, ord)
                    TaylorSeries.identity!(ηx[i, j], zero_q_1, ord)
                    TaylorSeries.identity!(ηy[i, j], zero_q_1, ord)
                    TaylorSeries.identity!(ηz[i, j], zero_q_1, ord)
                    TaylorSeries.identity!(ηx1[i, j], zero_q_1, ord)
                    TaylorSeries.identity!(ηy1[i, j], zero_q_1, ord)
                    TaylorSeries.identity!(ηz1[i, j], zero_q_1, ord)
                    TaylorSeries.identity!(ηx2[i, j], zero_q_1, ord)
                    TaylorSeries.identity!(ηy2[i, j], zero_q_1, ord)
                    TaylorSeries.identity!(ηz2[i, j], zero_q_1, ord)
                    TaylorSeries.identity!(ζx[i, j], zero_q_1, ord)
                    TaylorSeries.identity!(ζy[i, j], zero_q_1, ord)
                    TaylorSeries.identity!(ζz[i, j], zero_q_1, ord)
                    TaylorSeries.identity!(ζx1[i, j], zero_q_1, ord)
                    TaylorSeries.identity!(ζy1[i, j], zero_q_1, ord)
                    TaylorSeries.identity!(ζz1[i, j], zero_q_1, ord)
                    TaylorSeries.identity!(ζx2[i, j], zero_q_1, ord)
                    TaylorSeries.identity!(ζy2[i, j], zero_q_1, ord)
                    TaylorSeries.identity!(ζz2[i, j], zero_q_1, ord)
                end
            end
        end
        for j = _1_to_N
            for i = _1_to_N
                if i == j
                else
                    TaylorSeries.subst!(X[i, j], q[3i - 2], q[3j - 2], ord)
                    TaylorSeries.subst!(Y[i, j], q[3i - 1], q[3j - 1], ord)
                    TaylorSeries.subst!(Z[i, j], q[3i], q[3j], ord)
                    TaylorSeries.subst!(U[i, j], dq[3i - 2], dq[3j - 2], ord)
                    TaylorSeries.subst!(V[i, j], dq[3i - 1], dq[3j - 1], ord)
                    TaylorSeries.subst!(W[i, j], dq[3i], dq[3j], ord)
                    TaylorSeries.mul!(tmp849[3j - 2], 4, dq[3j - 2], ord)
                    TaylorSeries.mul!(tmp851[3i - 2], 3, dq[3i - 2], ord)
                    TaylorSeries.subst!(_4U_m_3X[i, j], tmp849[3j - 2], tmp851[3i - 2], ord)
                    TaylorSeries.mul!(tmp854[3j - 1], 4, dq[3j - 1], ord)
                    TaylorSeries.mul!(tmp856[3i - 1], 3, dq[3i - 1], ord)
                    TaylorSeries.subst!(_4V_m_3Y[i, j], tmp854[3j - 1], tmp856[3i - 1], ord)
                    TaylorSeries.mul!(tmp859[3j], 4, dq[3j], ord)
                    TaylorSeries.mul!(tmp861[3i], 3, dq[3i], ord)
                    TaylorSeries.subst!(_4W_m_3Z[i, j], tmp859[3j], tmp861[3i], ord)
                    TaylorSeries.mul!(pn2x[i, j], X[i, j], _4U_m_3X[i, j], ord)
                    TaylorSeries.mul!(pn2y[i, j], Y[i, j], _4V_m_3Y[i, j], ord)
                    TaylorSeries.mul!(pn2z[i, j], Z[i, j], _4W_m_3Z[i, j], ord)
                    TaylorSeries.mul!(UU[i, j], dq[3i - 2], dq[3j - 2], ord)
                    TaylorSeries.mul!(VV[i, j], dq[3i - 1], dq[3j - 1], ord)
                    TaylorSeries.mul!(WW[i, j], dq[3i], dq[3j], ord)
                    TaylorSeries.add!(tmp869[i, j], UU[i, j], VV[i, j], ord)
                    TaylorSeries.add!(vi_dot_vj[i, j], tmp869[i, j], WW[i, j], ord)
                    TaylorSeries.pow!(tmp872[i, j], X[i, j], 2, ord)
                    TaylorSeries.pow!(tmp874[i, j], Y[i, j], 2, ord)
                    TaylorSeries.add!(tmp875[i, j], tmp872[i, j], tmp874[i, j], ord)
                    TaylorSeries.pow!(tmp877[i, j], Z[i, j], 2, ord)
                    TaylorSeries.add!(r_p2[i, j], tmp875[i, j], tmp877[i, j], ord)
                    TaylorSeries.sqrt!(r_p1d2[i, j], r_p2[i, j], ord)
                    TaylorSeries.pow!(r_p3d2[i, j], r_p2[i, j], 1.5, ord)
                    TaylorSeries.pow!(r_p7d2[i, j], r_p2[i, j], 3.5, ord)
                    TaylorSeries.div!(newtonianCoeff[i, j], μ[i], r_p3d2[i, j], ord)
                    TaylorSeries.add!(tmp885[i, j], pn2x[i, j], pn2y[i, j], ord)
                    TaylorSeries.add!(tmp886[i, j], tmp885[i, j], pn2z[i, j], ord)
                    TaylorSeries.mul!(pn2[i, j], newtonianCoeff[i, j], tmp886[i, j], ord)
                    TaylorSeries.mul!(newton_acc_X[i, j], X[i, j], newtonianCoeff[i, j], ord)
                    TaylorSeries.mul!(newton_acc_Y[i, j], Y[i, j], newtonianCoeff[i, j], ord)
                    TaylorSeries.mul!(newton_acc_Z[i, j], Z[i, j], newtonianCoeff[i, j], ord)
                    TaylorSeries.div!(newtonian1b_Potential[i, j], μ[i], r_p1d2[i, j], ord)
                    TaylorSeries.mul!(pn3[i, j], 3.5, newtonian1b_Potential[i, j], ord)
                    TaylorSeries.mul!(U_t_pn2[i, j], pn2[i, j], U[i, j], ord)
                    TaylorSeries.mul!(V_t_pn2[i, j], pn2[i, j], V[i, j], ord)
                    TaylorSeries.mul!(W_t_pn2[i, j], pn2[i, j], W[i, j], ord)
                    if UJ_interaction[i, j]
                        TaylorSeries.mul!(t31[i, j], X[i, j], M_[1, 3, j], ord)
                        TaylorSeries.mul!(t32[i, j], Y[i, j], M_[2, 3, j], ord)
                        TaylorSeries.mul!(t33[i, j], Z[i, j], M_[3, 3, j], ord)
                        TaylorSeries.add!(tmp900[i, j], t31[i, j], t32[i, j], ord)
                        TaylorSeries.add!(r_sin_ϕ[i, j], tmp900[i, j], t33[i, j], ord)
                        TaylorSeries.div!(sin_ϕ[i, j], r_sin_ϕ[i, j], r_p1d2[i, j], ord)
                        TaylorSeries.asin!(ϕ[i, j], sin_ϕ[i, j], tmp1049[i, j], ord)
                        TaylorSeries.sincos!(tmp1050[i, j], cos_ϕ[i, j], ϕ[i, j], ord)
                        TaylorSeries.pow!(sin2_ϕ[i, j], sin_ϕ[i, j], 2, ord)
                        TaylorSeries.pow!(sin3_ϕ[i, j], sin_ϕ[i, j], 3, ord)
                        TaylorSeries.mul!(tmp910[i, j], 1.5, sin2_ϕ[i, j], ord)
                        TaylorSeries.subst!(P_2_sin_ϕ[i, j], tmp910[i, j], 0.5, ord)
                        TaylorSeries.mul!(∂P_2_sin_ϕ[i, j], 3, sin_ϕ[i, j], ord)
                        TaylorSeries.mul!(tmp916[i, j], -1.5, sin_ϕ[i, j], ord)
                        TaylorSeries.mul!(tmp918[i, j], 2.5, sin3_ϕ[i, j], ord)
                        TaylorSeries.add!(P_3_sin_ϕ[i, j], tmp916[i, j], tmp918[i, j], ord)
                        TaylorSeries.mul!(tmp922[i, j], 7.5, sin2_ϕ[i, j], ord)
                        TaylorSeries.add!(∂P_3_sin_ϕ[i, j], -1.5, tmp922[i, j], ord)
                        TaylorSeries.subst!(tmp924[j], Λ2[j], ord)
                        TaylorSeries.pow!(tmp926[i, j], r_p2[i, j], 2, ord)
                        TaylorSeries.div!(Λ2j_div_r4[i, j], tmp924[j], tmp926[i, j], ord)
                        TaylorSeries.subst!(tmp928[j], Λ3[j], ord)
                        TaylorSeries.pow!(tmp930[i, j], r_p1d2[i, j], 5, ord)
                        TaylorSeries.div!(Λ3j_div_r5[i, j], tmp928[j], tmp930[i, j], ord)
                        TaylorSeries.subst!(tmp932[i, j], cos_ϕ[i, j], ord)
                        TaylorSeries.mul!(m_c_ϕ_∂P_2[i, j], tmp932[i, j], ∂P_2_sin_ϕ[i, j], ord)
                        TaylorSeries.subst!(tmp934[i, j], cos_ϕ[i, j], ord)
                        TaylorSeries.mul!(m_c_ϕ_∂P_3[i, j], tmp934[i, j], ∂P_3_sin_ϕ[i, j], ord)
                        TaylorSeries.mul!(tmp937[i, j], Λ2j_div_r4[i, j], 3, ord)
                        TaylorSeries.mul!(F_J2_ξ[i, j], tmp937[i, j], P_2_sin_ϕ[i, j], ord)
                        TaylorSeries.mul!(F_J2_ζ[i, j], Λ2j_div_r4[i, j], m_c_ϕ_∂P_2[i, j], ord)
                        TaylorSeries.mul!(tmp941[i, j], Λ3j_div_r5[i, j], 4, ord)
                        TaylorSeries.mul!(F_J3_ξ[i, j], tmp941[i, j], P_3_sin_ϕ[i, j], ord)
                        TaylorSeries.mul!(F_J3_ζ[i, j], Λ3j_div_r5[i, j], m_c_ϕ_∂P_3[i, j], ord)
                        TaylorSeries.add!(F_J_ξ[i, j], F_J2_ξ[i, j], F_J3_ξ[i, j], ord)
                        TaylorSeries.add!(F_J_ζ[i, j], F_J2_ζ[i, j], F_J3_ζ[i, j], ord)
                        TaylorSeries.div!(ξx[i, j], X[i, j], r_p1d2[i, j], ord)
                        TaylorSeries.div!(ξy[i, j], Y[i, j], r_p1d2[i, j], ord)
                        TaylorSeries.div!(ξz[i, j], Z[i, j], r_p1d2[i, j], ord)
                        TaylorSeries.mul!(ηx1[i, j], M_[2, 3, j], ξz[i, j], ord)
                        TaylorSeries.mul!(ηy1[i, j], M_[3, 3, j], ξx[i, j], ord)
                        TaylorSeries.mul!(ηz1[i, j], M_[1, 3, j], ξy[i, j], ord)
                        TaylorSeries.mul!(ηx2[i, j], M_[3, 3, j], ξy[i, j], ord)
                        TaylorSeries.mul!(ηy2[i, j], M_[1, 3, j], ξz[i, j], ord)
                        TaylorSeries.mul!(ηz2[i, j], M_[2, 3, j], ξx[i, j], ord)
                        TaylorSeries.subst!(ηx[i, j], ηx1[i, j], ηx2[i, j], ord)
                        TaylorSeries.subst!(ηy[i, j], ηy1[i, j], ηy2[i, j], ord)
                        TaylorSeries.subst!(ηz[i, j], ηz1[i, j], ηz2[i, j], ord)
                        TaylorSeries.mul!(ζx1[i, j], ξy[i, j], ηz[i, j], ord)
                        TaylorSeries.mul!(ζy1[i, j], ξz[i, j], ηx[i, j], ord)
                        TaylorSeries.mul!(ζz1[i, j], ξx[i, j], ηy[i, j], ord)
                        TaylorSeries.mul!(ζx2[i, j], ξz[i, j], ηy[i, j], ord)
                        TaylorSeries.mul!(ζy2[i, j], ξx[i, j], ηz[i, j], ord)
                        TaylorSeries.mul!(ζz2[i, j], ξy[i, j], ηx[i, j], ord)
                        TaylorSeries.subst!(ζx[i, j], ζx1[i, j], ζx2[i, j], ord)
                        TaylorSeries.subst!(ζy[i, j], ζy1[i, j], ζy2[i, j], ord)
                        TaylorSeries.subst!(ζz[i, j], ζz1[i, j], ζz2[i, j], ord)
                        TaylorSeries.mul!(F_J2_x1[i, j], F_J_ξ[i, j], ξx[i, j], ord)
                        TaylorSeries.mul!(F_J2_y1[i, j], F_J_ξ[i, j], ξy[i, j], ord)
                        TaylorSeries.mul!(F_J2_z1[i, j], F_J_ξ[i, j], ξz[i, j], ord)
                        TaylorSeries.mul!(F_J2_x2[i, j], F_J_ζ[i, j], ζx[i, j], ord)
                        TaylorSeries.mul!(F_J2_y2[i, j], F_J_ζ[i, j], ζy[i, j], ord)
                        TaylorSeries.mul!(F_J2_z2[i, j], F_J_ζ[i, j], ζz[i, j], ord)
                        TaylorSeries.add!(F_J2_x[i, j], F_J2_x1[i, j], F_J2_x2[i, j], ord)
                        TaylorSeries.add!(F_J2_y[i, j], F_J2_y1[i, j], F_J2_y2[i, j], ord)
                        TaylorSeries.add!(F_J2_z[i, j], F_J2_z1[i, j], F_J2_z2[i, j], ord)
                    end
                end
            end
            TaylorSeries.pow!(tmp977[3j - 2], dq[3j - 2], 2, ord)
            TaylorSeries.pow!(tmp979[3j - 1], dq[3j - 1], 2, ord)
            TaylorSeries.add!(tmp980[3j - 2], tmp977[3j - 2], tmp979[3j - 1], ord)
            TaylorSeries.pow!(tmp982[3j], dq[3j], 2, ord)
            TaylorSeries.add!(v2[j], tmp980[3j - 2], tmp982[3j], ord)
        end
        for j = _1_to_N
            for i = _1_to_N
                if i == j
                else
                    TaylorSeries.add!(newtonianNb_Potential[j], newtonianNb_Potential[j], newtonian1b_Potential[i, j], ord)
                    TaylorSeries.add!(newtonX[j], newtonX[j], newton_acc_X[i, j], ord)
                    TaylorSeries.add!(newtonY[j], newtonY[j], newton_acc_Y[i, j], ord)
                    TaylorSeries.add!(newtonZ[j], newtonZ[j], newton_acc_Z[i, j], ord)
                    if UJ_interaction[i, j]
                        TaylorSeries.add!(accX[j], accX[j], μ[i] * F_J2_x[i, j], ord)
                        TaylorSeries.add!(accY[j], accY[j], μ[i] * F_J2_y[i, j], ord)
                        TaylorSeries.add!(accZ[j], accZ[j], μ[i] * F_J2_z[i, j], ord)
                        TaylorSeries.subst!(accX[i], accX[i], μ[j] * F_J2_x[i, j], ord)
                        TaylorSeries.subst!(accY[i], accY[i], μ[j] * F_J2_y[i, j], ord)
                        TaylorSeries.subst!(accZ[i], accZ[i], μ[j] * F_J2_z[i, j], ord)
                    end
                end
            end
        end
        for j = _1_to_N
            for i = _1_to_N
                if i == j
                else
                    TaylorSeries.mul!(_4ϕj[i, j], 4, newtonianNb_Potential[j], ord)
                    TaylorSeries.add!(ϕi_plus_4ϕj[i, j], newtonianNb_Potential[i], _4ϕj[i, j], ord)
                    TaylorSeries.mul!(tmp1004[i], 2, v2[i], ord)
                    TaylorSeries.add!(tmp1005[j], v2[j], tmp1004[i], ord)
                    TaylorSeries.mul!(tmp1007[i, j], 4, vi_dot_vj[i, j], ord)
                    TaylorSeries.subst!(sj2_plus_2si2_minus_4vivj[i, j], tmp1005[j], tmp1007[i, j], ord)
                    TaylorSeries.subst!(ϕs_and_vs[i, j], sj2_plus_2si2_minus_4vivj[i, j], ϕi_plus_4ϕj[i, j], ord)
                    TaylorSeries.mul!(Xij_t_Ui[i, j], X[i, j], dq[3i - 2], ord)
                    TaylorSeries.mul!(Yij_t_Vi[i, j], Y[i, j], dq[3i - 1], ord)
                    TaylorSeries.mul!(Zij_t_Wi[i, j], Z[i, j], dq[3i], ord)
                    TaylorSeries.add!(tmp1013[i, j], Xij_t_Ui[i, j], Yij_t_Vi[i, j], ord)
                    TaylorSeries.add!(Rij_dot_Vi[i, j], tmp1013[i, j], Zij_t_Wi[i, j], ord)
                    TaylorSeries.pow!(tmp1016[i, j], Rij_dot_Vi[i, j], 2, ord)
                    TaylorSeries.div!(pn1t7[i, j], tmp1016[i, j], r_p2[i, j], ord)
                    TaylorSeries.mul!(tmp1019[i, j], 1.5, pn1t7[i, j], ord)
                    TaylorSeries.subst!(pn1t2_7[i, j], ϕs_and_vs[i, j], tmp1019[i, j], ord)
                    TaylorSeries.add!(pn1t1_7[i, j], c_p2, pn1t2_7[i, j], ord)
                    for k = Base.OneTo(postnewton_iter)
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
            for k = Base.OneTo(postnewton_iter)
                TaylorSeries.identity!(pntempX[j, k], zero_q_1, ord)
                TaylorSeries.identity!(pntempY[j, k], zero_q_1, ord)
                TaylorSeries.identity!(pntempZ[j, k], zero_q_1, ord)
            end
        end
        for k = Base.OneTo(postnewton_iter)
            for j = _1_to_N
                for i = _1_to_N
                    if i == j
                    else
                        TaylorSeries.mul!(pNX_t_X[i, j, k], postNewtonX[i, k], X[i, j], ord)
                        TaylorSeries.mul!(pNY_t_Y[i, j, k], postNewtonY[i, k], Y[i, j], ord)
                        TaylorSeries.mul!(pNZ_t_Z[i, j, k], postNewtonZ[i, k], Z[i, j], ord)
                        TaylorSeries.add!(tmp1025[i, j, k], pNX_t_X[i, j, k], pNY_t_Y[i, j, k], ord)
                        TaylorSeries.add!(tmp1026[i, j, k], tmp1025[i, j, k], pNZ_t_Z[i, j, k], ord)
                        TaylorSeries.add!(pn1[i, j, k], pn1t1_7[i, j], tmp1026[i, j, k], ord)
                        TaylorSeries.mul!(X_t_pn1[i, j, k], newton_acc_X[i, j], pn1[i, j, k], ord)
                        TaylorSeries.mul!(Y_t_pn1[i, j, k], newton_acc_Y[i, j], pn1[i, j, k], ord)
                        TaylorSeries.mul!(Z_t_pn1[i, j, k], newton_acc_Z[i, j], pn1[i, j, k], ord)
                        TaylorSeries.mul!(pNX_t_pn3[i, j, k], postNewtonX[i, k], pn3[i, j], ord)
                        TaylorSeries.mul!(pNY_t_pn3[i, j, k], postNewtonY[i, k], pn3[i, j], ord)
                        TaylorSeries.mul!(pNZ_t_pn3[i, j, k], postNewtonZ[i, k], pn3[i, j], ord)
                    end
                end
            end
            for j = _1_to_N
                for i = _1_to_N
                    if i == j
                    else
                        TaylorSeries.add!(pntempX[j, k], pntempX[j, k], X_t_pn1[i, j, k] + (U_t_pn2[i, j] + pNX_t_pn3[i, j, k]), ord)
                        TaylorSeries.add!(pntempY[j, k], pntempY[j, k], Y_t_pn1[i, j, k] + (V_t_pn2[i, j] + pNY_t_pn3[i, j, k]), ord)
                        TaylorSeries.add!(pntempZ[j, k], pntempZ[j, k], Z_t_pn1[i, j, k] + (W_t_pn2[i, j] + pNZ_t_pn3[i, j, k]), ord)
                    end
                end
                TaylorSeries.mul!(postNewtonX[j, k + 1], pntempX[j, k], c_m2, ord)
                TaylorSeries.mul!(postNewtonY[j, k + 1], pntempY[j, k], c_m2, ord)
                TaylorSeries.mul!(postNewtonZ[j, k + 1], pntempZ[j, k], c_m2, ord)
            end
        end
        for i = _1_to_N
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

# Multi-threaded version of NBP_pN_A_J23E_J23M_J2S!
function NBP_pN_A_J23E_J23M_J2S_threads!(dq, q, params, t)
    local S = eltype(q[1])
    local N = Int((length(q))/6) # number of bodies, including NEA
    local _1_to_N = Base.OneTo(N) # iterator over all bodies

    #TODO: handle appropiately @taylorize'd version with postnewton_iter>1
    local postnewton_iter = 1 # number of iterations of post-Newtonian subroutine
    local j2_body_index = [su, ea, mo] # indices of bodies with J2 flattening (note: Earth and Moon also have J3)

    # parameters related to speed of light, c
    local c_p2 = 29979.063823897606 # c^2 = 29979.063823897606 au^2/d^2
    local c_m2 = 3.3356611996764786e-5 # c^-2 = 3.3356611996764786e-5 d^2/au^2

    local zero_q_1 = zero(q[1])

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

    Threads.@threads for j in _1_to_N
        for i in _1_to_N
            # i == j && continue
            if i == j
            else
                pn2[i,j] = zero_q_1
                pn3[i,j] = zero_q_1
                _4ϕj[i,j] = zero_q_1
                ϕi_plus_4ϕj[i,j] = zero_q_1
                sj2_plus_2si2_minus_4vivj[i,j] = zero_q_1
                ϕs_and_vs[i,j] = zero_q_1
                UU[i,j] = zero_q_1
                VV[i,j] = zero_q_1
                WW[i,j] = zero_q_1
                U_t_pn2[i,j] = zero_q_1
                V_t_pn2[i,j] = zero_q_1
                W_t_pn2[i,j] = zero_q_1
                vi_dot_vj[i,j] = zero_q_1
                newton_acc_X[i,j] = zero_q_1
                newton_acc_Y[i,j] = zero_q_1
                newton_acc_Z[i,j] = zero_q_1
                pn1t1_7[i,j] = zero_q_1
            end
        end
    end

    # J2 acceleration auxiliaries
    t31 = Array{Taylor1{S}}(undef, N, N)
    t32 = Array{Taylor1{S}}(undef, N, N)
    t33 = Array{Taylor1{S}}(undef, N, N)
    r_sin_ϕ = Array{Taylor1{S}}(undef, N, N)
    F_J2_x = Array{Taylor1{S}}(undef, N, N)
    F_J2_y = Array{Taylor1{S}}(undef, N, N)
    F_J2_z = Array{Taylor1{S}}(undef, N, N)
    F_J2_x1 = Array{Taylor1{S}}(undef, N, N)
    F_J2_y1 = Array{Taylor1{S}}(undef, N, N)
    F_J2_z1 = Array{Taylor1{S}}(undef, N, N)
    F_J2_x2 = Array{Taylor1{S}}(undef, N, N)
    F_J2_y2 = Array{Taylor1{S}}(undef, N, N)
    F_J2_z2 = Array{Taylor1{S}}(undef, N, N)
    sin_ϕ = Array{Taylor1{S}}(undef, N, N)
    sin2_ϕ = Array{Taylor1{S}}(undef, N, N)
    sin3_ϕ = Array{Taylor1{S}}(undef, N, N)
    sin4_ϕ = Array{Taylor1{S}}(undef, N, N)
    ϕ = Array{Taylor1{S}}(undef, N, N)
    cos_ϕ = Array{Taylor1{S}}(undef, N, N)
    P_2_sin_ϕ = Array{Taylor1{S}}(undef, N, N)
    ∂P_2_sin_ϕ = Array{Taylor1{S}}(undef, N, N)
    P_3_sin_ϕ = Array{Taylor1{S}}(undef, N, N)
    ∂P_3_sin_ϕ = Array{Taylor1{S}}(undef, N, N)
    m_c_ϕ_∂P_2 = Array{Taylor1{S}}(undef, N, N)
    m_c_ϕ_∂P_3 = Array{Taylor1{S}}(undef, N, N)
    Λ2j_div_r4 = Array{Taylor1{S}}(undef, N, N)
    Λ3j_div_r5 = Array{Taylor1{S}}(undef, N, N)
    F_J_ξ = Array{Taylor1{S}}(undef, N, N)
    F_J_η = Array{Taylor1{S}}(undef, N, N)
    F_J_ζ = Array{Taylor1{S}}(undef, N, N)
    F_J2_ξ = Array{Taylor1{S}}(undef, N, N)
    F_J2_η = Array{Taylor1{S}}(undef, N, N)
    F_J2_ζ = Array{Taylor1{S}}(undef, N, N)
    F_J3_ξ = Array{Taylor1{S}}(undef, N, N)
    F_J3_η = Array{Taylor1{S}}(undef, N, N)
    F_J3_ζ = Array{Taylor1{S}}(undef, N, N)
    ξx = Array{Taylor1{S}}(undef, N, N)
    ξy = Array{Taylor1{S}}(undef, N, N)
    ξz = Array{Taylor1{S}}(undef, N, N)
    ηx = Array{Taylor1{S}}(undef, N, N)
    ηy = Array{Taylor1{S}}(undef, N, N)
    ηz = Array{Taylor1{S}}(undef, N, N)
    ηx1 = Array{Taylor1{S}}(undef, N, N)
    ηy1 = Array{Taylor1{S}}(undef, N, N)
    ηz1 = Array{Taylor1{S}}(undef, N, N)
    ηx2 = Array{Taylor1{S}}(undef, N, N)
    ηy2 = Array{Taylor1{S}}(undef, N, N)
    ηz2 = Array{Taylor1{S}}(undef, N, N)
    ζx = Array{Taylor1{S}}(undef, N, N)
    ζy = Array{Taylor1{S}}(undef, N, N)
    ζz = Array{Taylor1{S}}(undef, N, N)
    ζx1 = Array{Taylor1{S}}(undef, N, N)
    ζy1 = Array{Taylor1{S}}(undef, N, N)
    ζz1 = Array{Taylor1{S}}(undef, N, N)
    ζx2 = Array{Taylor1{S}}(undef, N, N)
    ζy2 = Array{Taylor1{S}}(undef, N, N)
    ζz2 = Array{Taylor1{S}}(undef, N, N)

    # extended-body accelerations
    accX = Array{Taylor1{S}}(undef, N)
    accY = Array{Taylor1{S}}(undef, N)
    accZ = Array{Taylor1{S}}(undef, N)

    # rotations to and from Earth, Sun and Moon pole-oriented frames
    local dsj2k = t-2.451545e6 # J2000.0 = 2.451545e6
    local αs = deg2rad(α_p_sun*one(t))
    local δs = deg2rad(δ_p_sun*one(t))
    local αm = moon_pole_ra(dsj2k)
    local δm = moon_pole_dec(dsj2k)
    local M_ = Array{Taylor1{S}}(undef, 3, 3, N)
    local M_[:,:,ea] = t2c_jpl_de430(dsj2k)
    local M_[:,:,su] = pole_rotation( αs, δs )
    local M_[:,:,mo] = pole_rotation( αm, δm )

    for j in _1_to_N
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

    Threads.@threads for j in j2_body_index
        for i in _1_to_N
            if i == j
            else
                t31[i,j] = zero_q_1
                t32[i,j] = zero_q_1
                t33[i,j] = zero_q_1
                F_J2_x[i,j] = zero_q_1
                F_J2_y[i,j] = zero_q_1
                F_J2_z[i,j] = zero_q_1
                F_J2_x1[i,j] = zero_q_1
                F_J2_y1[i,j] = zero_q_1
                F_J2_z1[i,j] = zero_q_1
                F_J2_x2[i,j] = zero_q_1
                F_J2_y2[i,j] = zero_q_1
                F_J2_z2[i,j] = zero_q_1
                sin_ϕ[i,j] = zero_q_1
                sin2_ϕ[i,j] = zero_q_1
                sin3_ϕ[i,j] = zero_q_1
                sin4_ϕ[i,j] = zero_q_1
                ϕ[i,j] = zero_q_1
                cos_ϕ[i,j] = zero_q_1
                P_2_sin_ϕ[i,j] = zero_q_1
                ∂P_2_sin_ϕ[i,j] = zero_q_1
                P_3_sin_ϕ[i,j] = zero_q_1
                ∂P_3_sin_ϕ[i,j] = zero_q_1
                m_c_ϕ_∂P_2[i,j] = zero_q_1
                m_c_ϕ_∂P_3[i,j] = zero_q_1
                Λ2j_div_r4[i,j] = zero_q_1
                Λ3j_div_r5[i,j] = zero_q_1
                F_J_ξ[i,j] = zero_q_1
                F_J_η[i,j] = zero_q_1
                F_J_ζ[i,j] = zero_q_1
                F_J2_ξ[i,j] = zero_q_1
                F_J2_η[i,j] = zero_q_1
                F_J2_ζ[i,j] = zero_q_1
                F_J3_ξ[i,j] = zero_q_1
                F_J3_η[i,j] = zero_q_1
                F_J3_ζ[i,j] = zero_q_1
                ξx[i,j] = zero_q_1
                ξy[i,j] = zero_q_1
                ξz[i,j] = zero_q_1
                ηx[i,j] = zero_q_1
                ηy[i,j] = zero_q_1
                ηz[i,j] = zero_q_1
                ηx1[i,j] = zero_q_1
                ηy1[i,j] = zero_q_1
                ηz1[i,j] = zero_q_1
                ηx2[i,j] = zero_q_1
                ηy2[i,j] = zero_q_1
                ηz2[i,j] = zero_q_1
                ζx[i,j] = zero_q_1
                ζy[i,j] = zero_q_1
                ζz[i,j] = zero_q_1
                ζx1[i,j] = zero_q_1
                ζy1[i,j] = zero_q_1
                ζz1[i,j] = zero_q_1
                ζx2[i,j] = zero_q_1
                ζy2[i,j] = zero_q_1
                ζz2[i,j] = zero_q_1
            end #if i == j
        end #for i in _1_to_N
    end #for j in j2_body_index

    #compute point-mass Newtonian accelerations, all bodies
    Threads.@threads for j in _1_to_N
        for i in _1_to_N
            # i == j && continue
            if i == j
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

                #J2 accelerations, if j-th body is flattened
                if UJ_interaction[i,j]
                    # # rotate from inertial frame to extended-body frame
                    t31[i,j] = X[i,j]*M_[1,3,j]
                    t32[i,j] = Y[i,j]*M_[2,3,j]
                    t33[i,j] = Z[i,j]*M_[3,3,j]
                    r_sin_ϕ[i,j] = (t31[i,j]+t32[i,j])+t33[i,j]

                    # compute cartesian coordinates of acceleration due to body figure in body frame
                    sin_ϕ[i,j] = r_sin_ϕ[i,j]/r_p1d2[i,j] # z/r = latitude of point mass wrt body frame
                    ϕ[i,j] = asin(sin_ϕ[i,j])
                    cos_ϕ[i,j] = cos(ϕ[i,j])
                    sin2_ϕ[i,j] = sin_ϕ[i,j]^2
                    sin3_ϕ[i,j] = sin_ϕ[i,j]^3
                    P_2_sin_ϕ[i,j] = 1.5sin2_ϕ[i,j] - 0.5
                    ∂P_2_sin_ϕ[i,j] = 3sin_ϕ[i,j]
                    P_3_sin_ϕ[i,j] = (-1.5sin_ϕ[i,j]) + (2.5sin3_ϕ[i,j])
                    ∂P_3_sin_ϕ[i,j] = -1.5 + 7.5sin2_ϕ[i,j]
                    Λ2j_div_r4[i,j] = (-Λ2[j])/(r_p2[i,j]^2)
                    Λ3j_div_r5[i,j] = (-Λ3[j])/(r_p1d2[i,j]^5)
                    m_c_ϕ_∂P_2[i,j] = (-cos_ϕ[i,j])*∂P_2_sin_ϕ[i,j]
                    m_c_ϕ_∂P_3[i,j] = (-cos_ϕ[i,j])*∂P_3_sin_ϕ[i,j]
                    F_J2_ξ[i,j] = ( Λ2j_div_r4[i,j]*(3P_2_sin_ϕ[i,j]) )
                    #F_J2_η[i,j] = zero_q_1
                    F_J2_ζ[i,j] = Λ2j_div_r4[i,j]*m_c_ϕ_∂P_2[i,j]
                    F_J3_ξ[i,j] = ( Λ3j_div_r5[i,j]*(4P_3_sin_ϕ[i,j]) )
                    #F_J3_η[i,j] = zero_q_1
                    F_J3_ζ[i,j] = Λ3j_div_r5[i,j]*m_c_ϕ_∂P_3[i,j]
                    F_J_ξ[i,j] = F_J2_ξ[i,j] + F_J3_ξ[i,j]
                    #F_J_η[i,j] = zero_q_1
                    F_J_ζ[i,j] = F_J2_ζ[i,j] + F_J3_ζ[i,j]
                    #Compute unit vectors ξ,η,ζ
                    ξx[i,j] = X[i,j]/r_p1d2[i,j]
                    ξy[i,j] = Y[i,j]/r_p1d2[i,j]
                    ξz[i,j] = Z[i,j]/r_p1d2[i,j]
                    #Compute η = p x ξ
                    ηx1[i,j] = M_[2,3,j]*ξz[i,j]
                    ηy1[i,j] = M_[3,3,j]*ξx[i,j]
                    ηz1[i,j] = M_[1,3,j]*ξy[i,j]
                    ηx2[i,j] = M_[3,3,j]*ξy[i,j]
                    ηy2[i,j] = M_[1,3,j]*ξz[i,j]
                    ηz2[i,j] = M_[2,3,j]*ξx[i,j]
                    ηx[i,j] = ηx1[i,j] - ηx2[i,j]
                    ηy[i,j] = ηy1[i,j] - ηy2[i,j]
                    ηz[i,j] = ηz1[i,j] - ηz2[i,j]
                    #Compute ζ = ξ x η
                    ζx1[i,j] = ξy[i,j]*ηz[i,j]
                    ζy1[i,j] = ξz[i,j]*ηx[i,j]
                    ζz1[i,j] = ξx[i,j]*ηy[i,j]
                    ζx2[i,j] = ξz[i,j]*ηy[i,j]
                    ζy2[i,j] = ξx[i,j]*ηz[i,j]
                    ζz2[i,j] = ξy[i,j]*ηx[i,j]
                    ζx[i,j] = ζx1[i,j] - ζx2[i,j]
                    ζy[i,j] = ζy1[i,j] - ζy2[i,j]
                    ζz[i,j] = ζz1[i,j] - ζz2[i,j]
                    # compute cartesian coordinates of acceleration due to body figure in inertial frame
                    F_J2_x1[i,j] = F_J_ξ[i,j]*ξx[i,j]
                    F_J2_y1[i,j] = F_J_ξ[i,j]*ξy[i,j]
                    F_J2_z1[i,j] = F_J_ξ[i,j]*ξz[i,j]
                    F_J2_x2[i,j] = F_J_ζ[i,j]*ζx[i,j]
                    F_J2_y2[i,j] = F_J_ζ[i,j]*ζy[i,j]
                    F_J2_z2[i,j] = F_J_ζ[i,j]*ζz[i,j]
                    F_J2_x[i,j] = F_J2_x1[i,j] + F_J2_x2[i,j]
                    F_J2_y[i,j] = F_J2_y1[i,j] + F_J2_y2[i,j]
                    F_J2_z[i,j] = F_J2_z1[i,j] + F_J2_z2[i,j]
                end # if UJ_interaction[i,j]
            end #if i != j
        end #for, i
        v2[j] = ( (dq[3j-2]^2)+(dq[3j-1]^2) ) + (dq[3j]^2)
    end #for, j

    for j in _1_to_N
        for i in _1_to_N
            # i == j && continue
            if i == j
            else
                newtonianNb_Potential[j] = newtonianNb_Potential[j] + newtonian1b_Potential[i, j]
                newtonX[j] = newtonX[j] + newton_acc_X[i,j]
                newtonY[j] = newtonY[j] + newton_acc_Y[i,j]
                newtonZ[j] = newtonZ[j] + newton_acc_Z[i,j]
                if UJ_interaction[i,j]
                    # # add result to total acceleration on upon j-th body figure due to i-th point mass
                    # @show "acc",j,"+μ",i,"Λ2",j
                    accX[j] = accX[j] + (μ[i]*F_J2_x[i,j])
                    accY[j] = accY[j] + (μ[i]*F_J2_y[i,j])
                    accZ[j] = accZ[j] + (μ[i]*F_J2_z[i,j])
                    # # reaction force on i-th body
                    # @show "acc",i,"-μ",j,"Λ2",j
                    accX[i] = accX[i] - (μ[j]*F_J2_x[i,j])
                    accY[i] = accY[i] - (μ[j]*F_J2_y[i,j])
                    accZ[i] = accZ[i] - (μ[j]*F_J2_z[i,j])
                end
            end
        end
    end

    #post-Newtonian corrections to gravitational acceleration
    #Moyer, 1971, page 7 eq. 35
    # post-Newtonian iterative procedure setup and initialization
    Threads.@threads for j in _1_to_N
        for i in _1_to_N
            # i == j && continue
            if i == j
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
                for k in Base.OneTo(postnewton_iter)
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
            end
        end
        postNewtonX[j,1] = newtonX[j]
        postNewtonY[j,1] = newtonY[j]
        postNewtonZ[j,1] = newtonZ[j]
        for k in Base.OneTo(postnewton_iter)
            pntempX[j,k] = zero_q_1
            pntempY[j,k] = zero_q_1
            pntempZ[j,k] = zero_q_1
        end
    end

    # post-Newtonian iterations
    for k in Base.OneTo(postnewton_iter)
        Threads.@threads for j in _1_to_N
            for i in _1_to_N
                # i == j && continue
                if i == j
                else
                    pNX_t_X[i,j,k] = postNewtonX[i,k]*X[i,j]
                    pNY_t_Y[i,j,k] = postNewtonY[i,k]*Y[i,j]
                    pNZ_t_Z[i,j,k] = postNewtonZ[i,k]*Z[i,j]
                    pn1[i,j,k] = (  pn1t1_7[i,j]  +  ( (pNX_t_X[i,j,k]+pNY_t_Y[i,j,k]) + pNZ_t_Z[i,j,k] )  )

                    X_t_pn1[i,j,k] = newton_acc_X[i,j]*pn1[i,j,k]
                    Y_t_pn1[i,j,k] = newton_acc_Y[i,j]*pn1[i,j,k]
                    Z_t_pn1[i,j,k] = newton_acc_Z[i,j]*pn1[i,j,k]

                    pNX_t_pn3[i,j,k] = postNewtonX[i,k]*pn3[i,j]
                    pNY_t_pn3[i,j,k] = postNewtonY[i,k]*pn3[i,j]
                    pNZ_t_pn3[i,j,k] = postNewtonZ[i,k]*pn3[i,j]
                end # if i != j
            end #for i
        end #for j
        for j in _1_to_N
            for i in _1_to_N
                # i == j && continue
                if i == j
                else
                    pntempX[j,k] = pntempX[j,k] + ( X_t_pn1[i,j,k] + (U_t_pn2[i,j]+pNX_t_pn3[i,j,k]) )
                    pntempY[j,k] = pntempY[j,k] + ( Y_t_pn1[i,j,k] + (V_t_pn2[i,j]+pNY_t_pn3[i,j,k]) )
                    pntempZ[j,k] = pntempZ[j,k] + ( Z_t_pn1[i,j,k] + (W_t_pn2[i,j]+pNZ_t_pn3[i,j,k]) )
                end
            end
            postNewtonX[j,k+1] = pntempX[j,k]*c_m2
            postNewtonY[j,k+1] = pntempY[j,k]*c_m2
            postNewtonZ[j,k+1] = pntempZ[j,k]*c_m2
        end
    end #for k in Base.OneTo(postnewton_iter) # (post-Newtonian iterations)

    #fill accelerations (post-Newtonian and extended body accelerations)
    Threads.@threads for i in _1_to_N
        dq[3(N+i)-2] = postNewtonX[i,postnewton_iter+1] + accX[i]
        dq[3(N+i)-1] = postNewtonY[i,postnewton_iter+1] + accY[i]
        dq[3(N+i)  ] = postNewtonZ[i,postnewton_iter+1] + accZ[i]
    end

    nothing
end

function TaylorIntegration.jetcoeffs!(::Val{NBP_pN_A_J23E_J23M_J2S_threads!}, t::Taylor1{_T}, q::AbstractVector{Taylor1{_S}}, dq::AbstractVector{Taylor1{_S}}, params) where {_T <: Real, _S <: Number}
    order = t.order
    local S = eltype(q[1])
    local N = Int(length(q) / 6)
    local _1_to_N = Base.OneTo(N)
    local postnewton_iter = 1
    local j2_body_index = [su, ea, mo]
    local c_p2 = 29979.063823897606
    local c_m2 = 3.3356611996764786e-5
    local zero_q_1 = zero(q[1])
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
    #= In[16]:83 =# Threads.@threads for j = _1_to_N
            for i = _1_to_N
                if i == j
                else
                    pn2[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                    pn3[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                    _4ϕj[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                    ϕi_plus_4ϕj[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                    sj2_plus_2si2_minus_4vivj[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                    ϕs_and_vs[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                    UU[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                    VV[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                    WW[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                    U_t_pn2[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                    V_t_pn2[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                    W_t_pn2[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                    vi_dot_vj[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                    newton_acc_X[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                    newton_acc_Y[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                    newton_acc_Z[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                    pn1t1_7[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                end
            end
        end
    t31 = Array{Taylor1{S}}(undef, N, N)
    t32 = Array{Taylor1{S}}(undef, N, N)
    t33 = Array{Taylor1{S}}(undef, N, N)
    r_sin_ϕ = Array{Taylor1{S}}(undef, N, N)
    F_J2_x = Array{Taylor1{S}}(undef, N, N)
    F_J2_y = Array{Taylor1{S}}(undef, N, N)
    F_J2_z = Array{Taylor1{S}}(undef, N, N)
    F_J2_x1 = Array{Taylor1{S}}(undef, N, N)
    F_J2_y1 = Array{Taylor1{S}}(undef, N, N)
    F_J2_z1 = Array{Taylor1{S}}(undef, N, N)
    F_J2_x2 = Array{Taylor1{S}}(undef, N, N)
    F_J2_y2 = Array{Taylor1{S}}(undef, N, N)
    F_J2_z2 = Array{Taylor1{S}}(undef, N, N)
    sin_ϕ = Array{Taylor1{S}}(undef, N, N)
    sin2_ϕ = Array{Taylor1{S}}(undef, N, N)
    sin3_ϕ = Array{Taylor1{S}}(undef, N, N)
    sin4_ϕ = Array{Taylor1{S}}(undef, N, N)
    ϕ = Array{Taylor1{S}}(undef, N, N)
    cos_ϕ = Array{Taylor1{S}}(undef, N, N)
    P_2_sin_ϕ = Array{Taylor1{S}}(undef, N, N)
    ∂P_2_sin_ϕ = Array{Taylor1{S}}(undef, N, N)
    P_3_sin_ϕ = Array{Taylor1{S}}(undef, N, N)
    ∂P_3_sin_ϕ = Array{Taylor1{S}}(undef, N, N)
    m_c_ϕ_∂P_2 = Array{Taylor1{S}}(undef, N, N)
    m_c_ϕ_∂P_3 = Array{Taylor1{S}}(undef, N, N)
    Λ2j_div_r4 = Array{Taylor1{S}}(undef, N, N)
    Λ3j_div_r5 = Array{Taylor1{S}}(undef, N, N)
    F_J_ξ = Array{Taylor1{S}}(undef, N, N)
    F_J_η = Array{Taylor1{S}}(undef, N, N)
    F_J_ζ = Array{Taylor1{S}}(undef, N, N)
    F_J2_ξ = Array{Taylor1{S}}(undef, N, N)
    F_J2_η = Array{Taylor1{S}}(undef, N, N)
    F_J2_ζ = Array{Taylor1{S}}(undef, N, N)
    F_J3_ξ = Array{Taylor1{S}}(undef, N, N)
    F_J3_η = Array{Taylor1{S}}(undef, N, N)
    F_J3_ζ = Array{Taylor1{S}}(undef, N, N)
    ξx = Array{Taylor1{S}}(undef, N, N)
    ξy = Array{Taylor1{S}}(undef, N, N)
    ξz = Array{Taylor1{S}}(undef, N, N)
    ηx = Array{Taylor1{S}}(undef, N, N)
    ηy = Array{Taylor1{S}}(undef, N, N)
    ηz = Array{Taylor1{S}}(undef, N, N)
    ηx1 = Array{Taylor1{S}}(undef, N, N)
    ηy1 = Array{Taylor1{S}}(undef, N, N)
    ηz1 = Array{Taylor1{S}}(undef, N, N)
    ηx2 = Array{Taylor1{S}}(undef, N, N)
    ηy2 = Array{Taylor1{S}}(undef, N, N)
    ηz2 = Array{Taylor1{S}}(undef, N, N)
    ζx = Array{Taylor1{S}}(undef, N, N)
    ζy = Array{Taylor1{S}}(undef, N, N)
    ζz = Array{Taylor1{S}}(undef, N, N)
    ζx1 = Array{Taylor1{S}}(undef, N, N)
    ζy1 = Array{Taylor1{S}}(undef, N, N)
    ζz1 = Array{Taylor1{S}}(undef, N, N)
    ζx2 = Array{Taylor1{S}}(undef, N, N)
    ζy2 = Array{Taylor1{S}}(undef, N, N)
    ζz2 = Array{Taylor1{S}}(undef, N, N)
    accX = Array{Taylor1{S}}(undef, N)
    accY = Array{Taylor1{S}}(undef, N)
    accZ = Array{Taylor1{S}}(undef, N)
    local dsj2k = t - 2.451545e6
    local αs = deg2rad(α_p_sun * one(t))
    local δs = deg2rad(δ_p_sun * one(t))
    local αm = moon_pole_ra(dsj2k)
    local δm = moon_pole_dec(dsj2k)
    local M_ = Array{Taylor1{S}}(undef, 3, 3, N)
    local M_[:, :, ea] = t2c_jpl_de430(dsj2k)
    local M_[:, :, su] = pole_rotation(αs, δs)
    local M_[:, :, mo] = pole_rotation(αm, δm)
    for j = _1_to_N
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
    #= In[16]:200 =# Threads.@threads for j = j2_body_index
            for i = _1_to_N
                if i == j
                else
                    t31[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                    t32[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                    t33[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                    F_J2_x[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                    F_J2_y[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                    F_J2_z[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                    F_J2_x1[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                    F_J2_y1[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                    F_J2_z1[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                    F_J2_x2[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                    F_J2_y2[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                    F_J2_z2[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                    sin_ϕ[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                    sin2_ϕ[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                    sin3_ϕ[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                    sin4_ϕ[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                    ϕ[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                    cos_ϕ[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                    P_2_sin_ϕ[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                    ∂P_2_sin_ϕ[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                    P_3_sin_ϕ[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                    ∂P_3_sin_ϕ[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                    m_c_ϕ_∂P_2[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                    m_c_ϕ_∂P_3[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                    Λ2j_div_r4[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                    Λ3j_div_r5[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                    F_J_ξ[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                    F_J_η[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                    F_J_ζ[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                    F_J2_ξ[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                    F_J2_η[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                    F_J2_ζ[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                    F_J3_ξ[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                    F_J3_η[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                    F_J3_ζ[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                    ξx[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                    ξy[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                    ξz[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                    ηx[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                    ηy[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                    ηz[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                    ηx1[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                    ηy1[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                    ηz1[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                    ηx2[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                    ηy2[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                    ηz2[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                    ζx[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                    ζy[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                    ζz[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                    ζx1[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                    ζy1[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                    ζz1[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                    ζx2[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                    ζy2[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                    ζz2[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                end
            end
        end
    tmp1530 = Array{Taylor1{_S}}(undef, size(dq))
    tmp1530 .= Taylor1(zero(_S), order)
    tmp1532 = Array{Taylor1{_S}}(undef, size(dq))
    tmp1532 .= Taylor1(zero(_S), order)
    tmp1535 = Array{Taylor1{_S}}(undef, size(dq))
    tmp1535 .= Taylor1(zero(_S), order)
    tmp1537 = Array{Taylor1{_S}}(undef, size(dq))
    tmp1537 .= Taylor1(zero(_S), order)
    tmp1540 = Array{Taylor1{_S}}(undef, size(dq))
    tmp1540 .= Taylor1(zero(_S), order)
    tmp1542 = Array{Taylor1{_S}}(undef, size(dq))
    tmp1542 .= Taylor1(zero(_S), order)
    pn2x = Array{Taylor1{_S}}(undef, size(X))
    pn2x .= Taylor1(zero(_S), order)
    pn2y = Array{Taylor1{_S}}(undef, size(Y))
    pn2y .= Taylor1(zero(_S), order)
    pn2z = Array{Taylor1{_S}}(undef, size(Z))
    pn2z .= Taylor1(zero(_S), order)
    tmp1550 = Array{Taylor1{_S}}(undef, size(UU))
    tmp1550 .= Taylor1(zero(_S), order)
    tmp1553 = Array{Taylor1{_S}}(undef, size(X))
    tmp1553 .= Taylor1(zero(_S), order)
    tmp1555 = Array{Taylor1{_S}}(undef, size(Y))
    tmp1555 .= Taylor1(zero(_S), order)
    tmp1556 = Array{Taylor1{_S}}(undef, size(tmp1553))
    tmp1556 .= Taylor1(zero(_S), order)
    tmp1558 = Array{Taylor1{_S}}(undef, size(Z))
    tmp1558 .= Taylor1(zero(_S), order)
    tmp1566 = Array{Taylor1{_S}}(undef, size(pn2x))
    tmp1566 .= Taylor1(zero(_S), order)
    tmp1567 = Array{Taylor1{_S}}(undef, size(tmp1566))
    tmp1567 .= Taylor1(zero(_S), order)
    tmp1581 = Array{Taylor1{_S}}(undef, size(t31))
    tmp1581 .= Taylor1(zero(_S), order)
    tmp1730 = Array{Taylor1{_S}}(undef, size(sin_ϕ))
    tmp1730 .= Taylor1(zero(_S), order)
    tmp1731 = Array{Taylor1{_S}}(undef, size(ϕ))
    tmp1731 .= Taylor1(zero(_S), order)
    tmp1591 = Array{Taylor1{_S}}(undef, size(sin2_ϕ))
    tmp1591 .= Taylor1(zero(_S), order)
    tmp1597 = Array{Taylor1{_S}}(undef, size(sin_ϕ))
    tmp1597 .= Taylor1(zero(_S), order)
    tmp1599 = Array{Taylor1{_S}}(undef, size(sin3_ϕ))
    tmp1599 .= Taylor1(zero(_S), order)
    tmp1603 = Array{Taylor1{_S}}(undef, size(sin2_ϕ))
    tmp1603 .= Taylor1(zero(_S), order)
    tmp1605 = Array{Taylor1{_S}}(undef, size(Λ2))
    tmp1605 .= Taylor1(zero(_S), order)
    tmp1607 = Array{Taylor1{_S}}(undef, size(r_p2))
    tmp1607 .= Taylor1(zero(_S), order)
    tmp1609 = Array{Taylor1{_S}}(undef, size(Λ3))
    tmp1609 .= Taylor1(zero(_S), order)
    tmp1611 = Array{Taylor1{_S}}(undef, size(r_p1d2))
    tmp1611 .= Taylor1(zero(_S), order)
    tmp1613 = Array{Taylor1{_S}}(undef, size(cos_ϕ))
    tmp1613 .= Taylor1(zero(_S), order)
    tmp1615 = Array{Taylor1{_S}}(undef, size(cos_ϕ))
    tmp1615 .= Taylor1(zero(_S), order)
    tmp1618 = Array{Taylor1{_S}}(undef, size(Λ2j_div_r4))
    tmp1618 .= Taylor1(zero(_S), order)
    tmp1622 = Array{Taylor1{_S}}(undef, size(Λ3j_div_r5))
    tmp1622 .= Taylor1(zero(_S), order)
    tmp1658 = Array{Taylor1{_S}}(undef, size(dq))
    tmp1658 .= Taylor1(zero(_S), order)
    tmp1660 = Array{Taylor1{_S}}(undef, size(dq))
    tmp1660 .= Taylor1(zero(_S), order)
    tmp1661 = Array{Taylor1{_S}}(undef, size(tmp1658))
    tmp1661 .= Taylor1(zero(_S), order)
    tmp1663 = Array{Taylor1{_S}}(undef, size(dq))
    tmp1663 .= Taylor1(zero(_S), order)
    #= In[16]:265 =# Threads.@threads for j = _1_to_N
            for i = _1_to_N
                if i == j
                else
                    X[i, j] = Taylor1(constant_term(q[3i - 2]) - constant_term(q[3j - 2]), order)
                    Y[i, j] = Taylor1(constant_term(q[3i - 1]) - constant_term(q[3j - 1]), order)
                    Z[i, j] = Taylor1(constant_term(q[3i]) - constant_term(q[3j]), order)
                    U[i, j] = Taylor1(constant_term(dq[3i - 2]) - constant_term(dq[3j - 2]), order)
                    V[i, j] = Taylor1(constant_term(dq[3i - 1]) - constant_term(dq[3j - 1]), order)
                    W[i, j] = Taylor1(constant_term(dq[3i]) - constant_term(dq[3j]), order)
                    tmp1530[3j - 2] = Taylor1(constant_term(4) * constant_term(dq[3j - 2]), order)
                    tmp1532[3i - 2] = Taylor1(constant_term(3) * constant_term(dq[3i - 2]), order)
                    _4U_m_3X[i, j] = Taylor1(constant_term(tmp1530[3j - 2]) - constant_term(tmp1532[3i - 2]), order)
                    tmp1535[3j - 1] = Taylor1(constant_term(4) * constant_term(dq[3j - 1]), order)
                    tmp1537[3i - 1] = Taylor1(constant_term(3) * constant_term(dq[3i - 1]), order)
                    _4V_m_3Y[i, j] = Taylor1(constant_term(tmp1535[3j - 1]) - constant_term(tmp1537[3i - 1]), order)
                    tmp1540[3j] = Taylor1(constant_term(4) * constant_term(dq[3j]), order)
                    tmp1542[3i] = Taylor1(constant_term(3) * constant_term(dq[3i]), order)
                    _4W_m_3Z[i, j] = Taylor1(constant_term(tmp1540[3j]) - constant_term(tmp1542[3i]), order)
                    pn2x[i, j] = Taylor1(constant_term(X[i, j]) * constant_term(_4U_m_3X[i, j]), order)
                    pn2y[i, j] = Taylor1(constant_term(Y[i, j]) * constant_term(_4V_m_3Y[i, j]), order)
                    pn2z[i, j] = Taylor1(constant_term(Z[i, j]) * constant_term(_4W_m_3Z[i, j]), order)
                    UU[i, j] = Taylor1(constant_term(dq[3i - 2]) * constant_term(dq[3j - 2]), order)
                    VV[i, j] = Taylor1(constant_term(dq[3i - 1]) * constant_term(dq[3j - 1]), order)
                    WW[i, j] = Taylor1(constant_term(dq[3i]) * constant_term(dq[3j]), order)
                    tmp1550[i, j] = Taylor1(constant_term(UU[i, j]) + constant_term(VV[i, j]), order)
                    vi_dot_vj[i, j] = Taylor1(constant_term(tmp1550[i, j]) + constant_term(WW[i, j]), order)
                    tmp1553[i, j] = Taylor1(constant_term(X[i, j]) ^ constant_term(2), order)
                    tmp1555[i, j] = Taylor1(constant_term(Y[i, j]) ^ constant_term(2), order)
                    tmp1556[i, j] = Taylor1(constant_term(tmp1553[i, j]) + constant_term(tmp1555[i, j]), order)
                    tmp1558[i, j] = Taylor1(constant_term(Z[i, j]) ^ constant_term(2), order)
                    r_p2[i, j] = Taylor1(constant_term(tmp1556[i, j]) + constant_term(tmp1558[i, j]), order)
                    r_p1d2[i, j] = Taylor1(sqrt(constant_term(r_p2[i, j])), order)
                    r_p3d2[i, j] = Taylor1(constant_term(r_p2[i, j]) ^ constant_term(1.5), order)
                    r_p7d2[i, j] = Taylor1(constant_term(r_p2[i, j]) ^ constant_term(3.5), order)
                    newtonianCoeff[i, j] = Taylor1(constant_term(μ[i]) / constant_term(r_p3d2[i, j]), order)
                    tmp1566[i, j] = Taylor1(constant_term(pn2x[i, j]) + constant_term(pn2y[i, j]), order)
                    tmp1567[i, j] = Taylor1(constant_term(tmp1566[i, j]) + constant_term(pn2z[i, j]), order)
                    pn2[i, j] = Taylor1(constant_term(newtonianCoeff[i, j]) * constant_term(tmp1567[i, j]), order)
                    newton_acc_X[i, j] = Taylor1(constant_term(X[i, j]) * constant_term(newtonianCoeff[i, j]), order)
                    newton_acc_Y[i, j] = Taylor1(constant_term(Y[i, j]) * constant_term(newtonianCoeff[i, j]), order)
                    newton_acc_Z[i, j] = Taylor1(constant_term(Z[i, j]) * constant_term(newtonianCoeff[i, j]), order)
                    newtonian1b_Potential[i, j] = Taylor1(constant_term(μ[i]) / constant_term(r_p1d2[i, j]), order)
                    pn3[i, j] = Taylor1(constant_term(3.5) * constant_term(newtonian1b_Potential[i, j]), order)
                    U_t_pn2[i, j] = Taylor1(constant_term(pn2[i, j]) * constant_term(U[i, j]), order)
                    V_t_pn2[i, j] = Taylor1(constant_term(pn2[i, j]) * constant_term(V[i, j]), order)
                    W_t_pn2[i, j] = Taylor1(constant_term(pn2[i, j]) * constant_term(W[i, j]), order)
                    if UJ_interaction[i, j]
                        t31[i, j] = Taylor1(constant_term(X[i, j]) * constant_term(M_[1, 3, j]), order)
                        t32[i, j] = Taylor1(constant_term(Y[i, j]) * constant_term(M_[2, 3, j]), order)
                        t33[i, j] = Taylor1(constant_term(Z[i, j]) * constant_term(M_[3, 3, j]), order)
                        tmp1581[i, j] = Taylor1(constant_term(t31[i, j]) + constant_term(t32[i, j]), order)
                        r_sin_ϕ[i, j] = Taylor1(constant_term(tmp1581[i, j]) + constant_term(t33[i, j]), order)
                        sin_ϕ[i, j] = Taylor1(constant_term(r_sin_ϕ[i, j]) / constant_term(r_p1d2[i, j]), order)
                        ϕ[i, j] = Taylor1(asin(constant_term(sin_ϕ[i, j])), order)
                        tmp1730[i, j] = Taylor1(sqrt(1 - constant_term(sin_ϕ[i, j]) ^ 2), order)
                        cos_ϕ[i, j] = Taylor1(cos(constant_term(ϕ[i, j])), order)
                        tmp1731[i, j] = Taylor1(sin(constant_term(ϕ[i, j])), order)
                        sin2_ϕ[i, j] = Taylor1(constant_term(sin_ϕ[i, j]) ^ constant_term(2), order)
                        sin3_ϕ[i, j] = Taylor1(constant_term(sin_ϕ[i, j]) ^ constant_term(3), order)
                        tmp1591[i, j] = Taylor1(constant_term(1.5) * constant_term(sin2_ϕ[i, j]), order)
                        P_2_sin_ϕ[i, j] = Taylor1(constant_term(tmp1591[i, j]) - constant_term(0.5), order)
                        ∂P_2_sin_ϕ[i, j] = Taylor1(constant_term(3) * constant_term(sin_ϕ[i, j]), order)
                        tmp1597[i, j] = Taylor1(constant_term(-1.5) * constant_term(sin_ϕ[i, j]), order)
                        tmp1599[i, j] = Taylor1(constant_term(2.5) * constant_term(sin3_ϕ[i, j]), order)
                        P_3_sin_ϕ[i, j] = Taylor1(constant_term(tmp1597[i, j]) + constant_term(tmp1599[i, j]), order)
                        tmp1603[i, j] = Taylor1(constant_term(7.5) * constant_term(sin2_ϕ[i, j]), order)
                        ∂P_3_sin_ϕ[i, j] = Taylor1(constant_term(-1.5) + constant_term(tmp1603[i, j]), order)
                        tmp1605[j] = Taylor1(-(constant_term(Λ2[j])), order)
                        tmp1607[i, j] = Taylor1(constant_term(r_p2[i, j]) ^ constant_term(2), order)
                        Λ2j_div_r4[i, j] = Taylor1(constant_term(tmp1605[j]) / constant_term(tmp1607[i, j]), order)
                        tmp1609[j] = Taylor1(-(constant_term(Λ3[j])), order)
                        tmp1611[i, j] = Taylor1(constant_term(r_p1d2[i, j]) ^ constant_term(5), order)
                        Λ3j_div_r5[i, j] = Taylor1(constant_term(tmp1609[j]) / constant_term(tmp1611[i, j]), order)
                        tmp1613[i, j] = Taylor1(-(constant_term(cos_ϕ[i, j])), order)
                        m_c_ϕ_∂P_2[i, j] = Taylor1(constant_term(tmp1613[i, j]) * constant_term(∂P_2_sin_ϕ[i, j]), order)
                        tmp1615[i, j] = Taylor1(-(constant_term(cos_ϕ[i, j])), order)
                        m_c_ϕ_∂P_3[i, j] = Taylor1(constant_term(tmp1615[i, j]) * constant_term(∂P_3_sin_ϕ[i, j]), order)
                        tmp1618[i, j] = Taylor1(constant_term(Λ2j_div_r4[i, j]) * constant_term(3), order)
                        F_J2_ξ[i, j] = Taylor1(constant_term(tmp1618[i, j]) * constant_term(P_2_sin_ϕ[i, j]), order)
                        F_J2_ζ[i, j] = Taylor1(constant_term(Λ2j_div_r4[i, j]) * constant_term(m_c_ϕ_∂P_2[i, j]), order)
                        tmp1622[i, j] = Taylor1(constant_term(Λ3j_div_r5[i, j]) * constant_term(4), order)
                        F_J3_ξ[i, j] = Taylor1(constant_term(tmp1622[i, j]) * constant_term(P_3_sin_ϕ[i, j]), order)
                        F_J3_ζ[i, j] = Taylor1(constant_term(Λ3j_div_r5[i, j]) * constant_term(m_c_ϕ_∂P_3[i, j]), order)
                        F_J_ξ[i, j] = Taylor1(constant_term(F_J2_ξ[i, j]) + constant_term(F_J3_ξ[i, j]), order)
                        F_J_ζ[i, j] = Taylor1(constant_term(F_J2_ζ[i, j]) + constant_term(F_J3_ζ[i, j]), order)
                        ξx[i, j] = Taylor1(constant_term(X[i, j]) / constant_term(r_p1d2[i, j]), order)
                        ξy[i, j] = Taylor1(constant_term(Y[i, j]) / constant_term(r_p1d2[i, j]), order)
                        ξz[i, j] = Taylor1(constant_term(Z[i, j]) / constant_term(r_p1d2[i, j]), order)
                        ηx1[i, j] = Taylor1(constant_term(M_[2, 3, j]) * constant_term(ξz[i, j]), order)
                        ηy1[i, j] = Taylor1(constant_term(M_[3, 3, j]) * constant_term(ξx[i, j]), order)
                        ηz1[i, j] = Taylor1(constant_term(M_[1, 3, j]) * constant_term(ξy[i, j]), order)
                        ηx2[i, j] = Taylor1(constant_term(M_[3, 3, j]) * constant_term(ξy[i, j]), order)
                        ηy2[i, j] = Taylor1(constant_term(M_[1, 3, j]) * constant_term(ξz[i, j]), order)
                        ηz2[i, j] = Taylor1(constant_term(M_[2, 3, j]) * constant_term(ξx[i, j]), order)
                        ηx[i, j] = Taylor1(constant_term(ηx1[i, j]) - constant_term(ηx2[i, j]), order)
                        ηy[i, j] = Taylor1(constant_term(ηy1[i, j]) - constant_term(ηy2[i, j]), order)
                        ηz[i, j] = Taylor1(constant_term(ηz1[i, j]) - constant_term(ηz2[i, j]), order)
                        ζx1[i, j] = Taylor1(constant_term(ξy[i, j]) * constant_term(ηz[i, j]), order)
                        ζy1[i, j] = Taylor1(constant_term(ξz[i, j]) * constant_term(ηx[i, j]), order)
                        ζz1[i, j] = Taylor1(constant_term(ξx[i, j]) * constant_term(ηy[i, j]), order)
                        ζx2[i, j] = Taylor1(constant_term(ξz[i, j]) * constant_term(ηy[i, j]), order)
                        ζy2[i, j] = Taylor1(constant_term(ξx[i, j]) * constant_term(ηz[i, j]), order)
                        ζz2[i, j] = Taylor1(constant_term(ξy[i, j]) * constant_term(ηx[i, j]), order)
                        ζx[i, j] = Taylor1(constant_term(ζx1[i, j]) - constant_term(ζx2[i, j]), order)
                        ζy[i, j] = Taylor1(constant_term(ζy1[i, j]) - constant_term(ζy2[i, j]), order)
                        ζz[i, j] = Taylor1(constant_term(ζz1[i, j]) - constant_term(ζz2[i, j]), order)
                        F_J2_x1[i, j] = Taylor1(constant_term(F_J_ξ[i, j]) * constant_term(ξx[i, j]), order)
                        F_J2_y1[i, j] = Taylor1(constant_term(F_J_ξ[i, j]) * constant_term(ξy[i, j]), order)
                        F_J2_z1[i, j] = Taylor1(constant_term(F_J_ξ[i, j]) * constant_term(ξz[i, j]), order)
                        F_J2_x2[i, j] = Taylor1(constant_term(F_J_ζ[i, j]) * constant_term(ζx[i, j]), order)
                        F_J2_y2[i, j] = Taylor1(constant_term(F_J_ζ[i, j]) * constant_term(ζy[i, j]), order)
                        F_J2_z2[i, j] = Taylor1(constant_term(F_J_ζ[i, j]) * constant_term(ζz[i, j]), order)
                        F_J2_x[i, j] = Taylor1(constant_term(F_J2_x1[i, j]) + constant_term(F_J2_x2[i, j]), order)
                        F_J2_y[i, j] = Taylor1(constant_term(F_J2_y1[i, j]) + constant_term(F_J2_y2[i, j]), order)
                        F_J2_z[i, j] = Taylor1(constant_term(F_J2_z1[i, j]) + constant_term(F_J2_z2[i, j]), order)
                    end
                end
            end
            tmp1658[3j - 2] = Taylor1(constant_term(dq[3j - 2]) ^ constant_term(2), order)
            tmp1660[3j - 1] = Taylor1(constant_term(dq[3j - 1]) ^ constant_term(2), order)
            tmp1661[3j - 2] = Taylor1(constant_term(tmp1658[3j - 2]) + constant_term(tmp1660[3j - 1]), order)
            tmp1663[3j] = Taylor1(constant_term(dq[3j]) ^ constant_term(2), order)
            v2[j] = Taylor1(constant_term(tmp1661[3j - 2]) + constant_term(tmp1663[3j]), order)
        end
    for j = _1_to_N
        for i = _1_to_N
            if i == j
            else
                newtonianNb_Potential[j] = Taylor1(constant_term(newtonianNb_Potential[j]) + constant_term(newtonian1b_Potential[i, j]), order)
                newtonX[j] = Taylor1(constant_term(newtonX[j]) + constant_term(newton_acc_X[i, j]), order)
                newtonY[j] = Taylor1(constant_term(newtonY[j]) + constant_term(newton_acc_Y[i, j]), order)
                newtonZ[j] = Taylor1(constant_term(newtonZ[j]) + constant_term(newton_acc_Z[i, j]), order)
                if UJ_interaction[i, j]
                    accX[j] = Taylor1(constant_term(accX[j]) + constant_term(μ[i] * F_J2_x[i, j]), order)
                    accY[j] = Taylor1(constant_term(accY[j]) + constant_term(μ[i] * F_J2_y[i, j]), order)
                    accZ[j] = Taylor1(constant_term(accZ[j]) + constant_term(μ[i] * F_J2_z[i, j]), order)
                    accX[i] = Taylor1(constant_term(accX[i]) - constant_term(μ[j] * F_J2_x[i, j]), order)
                    accY[i] = Taylor1(constant_term(accY[i]) - constant_term(μ[j] * F_J2_y[i, j]), order)
                    accZ[i] = Taylor1(constant_term(accZ[i]) - constant_term(μ[j] * F_J2_z[i, j]), order)
                end
            end
        end
    end
    tmp1685 = Array{Taylor1{_S}}(undef, size(v2))
    tmp1685 .= Taylor1(zero(_S), order)
    tmp1686 = Array{Taylor1{_S}}(undef, size(v2))
    tmp1686 .= Taylor1(zero(_S), order)
    tmp1688 = Array{Taylor1{_S}}(undef, size(vi_dot_vj))
    tmp1688 .= Taylor1(zero(_S), order)
    Xij_t_Ui = Array{Taylor1{_S}}(undef, size(X))
    Xij_t_Ui .= Taylor1(zero(_S), order)
    Yij_t_Vi = Array{Taylor1{_S}}(undef, size(Y))
    Yij_t_Vi .= Taylor1(zero(_S), order)
    Zij_t_Wi = Array{Taylor1{_S}}(undef, size(Z))
    Zij_t_Wi .= Taylor1(zero(_S), order)
    tmp1694 = Array{Taylor1{_S}}(undef, size(Xij_t_Ui))
    tmp1694 .= Taylor1(zero(_S), order)
    Rij_dot_Vi = Array{Taylor1{_S}}(undef, size(tmp1694))
    Rij_dot_Vi .= Taylor1(zero(_S), order)
    tmp1697 = Array{Taylor1{_S}}(undef, size(Rij_dot_Vi))
    tmp1697 .= Taylor1(zero(_S), order)
    pn1t7 = Array{Taylor1{_S}}(undef, size(tmp1697))
    pn1t7 .= Taylor1(zero(_S), order)
    tmp1700 = Array{Taylor1{_S}}(undef, size(pn1t7))
    tmp1700 .= Taylor1(zero(_S), order)
    pn1t2_7 = Array{Taylor1{_S}}(undef, size(ϕs_and_vs))
    pn1t2_7 .= Taylor1(zero(_S), order)
    #= In[16]:411 =# Threads.@threads for j = _1_to_N
            for i = _1_to_N
                if i == j
                else
                    _4ϕj[i, j] = Taylor1(constant_term(4) * constant_term(newtonianNb_Potential[j]), order)
                    ϕi_plus_4ϕj[i, j] = Taylor1(constant_term(newtonianNb_Potential[i]) + constant_term(_4ϕj[i, j]), order)
                    tmp1685[i] = Taylor1(constant_term(2) * constant_term(v2[i]), order)
                    tmp1686[j] = Taylor1(constant_term(v2[j]) + constant_term(tmp1685[i]), order)
                    tmp1688[i, j] = Taylor1(constant_term(4) * constant_term(vi_dot_vj[i, j]), order)
                    sj2_plus_2si2_minus_4vivj[i, j] = Taylor1(constant_term(tmp1686[j]) - constant_term(tmp1688[i, j]), order)
                    ϕs_and_vs[i, j] = Taylor1(constant_term(sj2_plus_2si2_minus_4vivj[i, j]) - constant_term(ϕi_plus_4ϕj[i, j]), order)
                    Xij_t_Ui[i, j] = Taylor1(constant_term(X[i, j]) * constant_term(dq[3i - 2]), order)
                    Yij_t_Vi[i, j] = Taylor1(constant_term(Y[i, j]) * constant_term(dq[3i - 1]), order)
                    Zij_t_Wi[i, j] = Taylor1(constant_term(Z[i, j]) * constant_term(dq[3i]), order)
                    tmp1694[i, j] = Taylor1(constant_term(Xij_t_Ui[i, j]) + constant_term(Yij_t_Vi[i, j]), order)
                    Rij_dot_Vi[i, j] = Taylor1(constant_term(tmp1694[i, j]) + constant_term(Zij_t_Wi[i, j]), order)
                    tmp1697[i, j] = Taylor1(constant_term(Rij_dot_Vi[i, j]) ^ constant_term(2), order)
                    pn1t7[i, j] = Taylor1(constant_term(tmp1697[i, j]) / constant_term(r_p2[i, j]), order)
                    tmp1700[i, j] = Taylor1(constant_term(1.5) * constant_term(pn1t7[i, j]), order)
                    pn1t2_7[i, j] = Taylor1(constant_term(ϕs_and_vs[i, j]) - constant_term(tmp1700[i, j]), order)
                    pn1t1_7[i, j] = Taylor1(constant_term(c_p2) + constant_term(pn1t2_7[i, j]), order)
                    for k = Base.OneTo(postnewton_iter)
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
            for k = Base.OneTo(postnewton_iter)
                pntempX[j, k] = Taylor1(identity(constant_term(zero_q_1)), order)
                pntempY[j, k] = Taylor1(identity(constant_term(zero_q_1)), order)
                pntempZ[j, k] = Taylor1(identity(constant_term(zero_q_1)), order)
            end
        end
    tmp1706 = Array{Taylor1{_S}}(undef, size(pNX_t_X))
    tmp1706 .= Taylor1(zero(_S), order)
    tmp1707 = Array{Taylor1{_S}}(undef, size(tmp1706))
    tmp1707 .= Taylor1(zero(_S), order)
    for k = Base.OneTo(postnewton_iter)
        #= In[16]:455 =# Threads.@threads for j = _1_to_N
                for i = _1_to_N
                    if i == j
                    else
                        pNX_t_X[i, j, k] = Taylor1(constant_term(postNewtonX[i, k]) * constant_term(X[i, j]), order)
                        pNY_t_Y[i, j, k] = Taylor1(constant_term(postNewtonY[i, k]) * constant_term(Y[i, j]), order)
                        pNZ_t_Z[i, j, k] = Taylor1(constant_term(postNewtonZ[i, k]) * constant_term(Z[i, j]), order)
                        tmp1706[i, j, k] = Taylor1(constant_term(pNX_t_X[i, j, k]) + constant_term(pNY_t_Y[i, j, k]), order)
                        tmp1707[i, j, k] = Taylor1(constant_term(tmp1706[i, j, k]) + constant_term(pNZ_t_Z[i, j, k]), order)
                        pn1[i, j, k] = Taylor1(constant_term(pn1t1_7[i, j]) + constant_term(tmp1707[i, j, k]), order)
                        X_t_pn1[i, j, k] = Taylor1(constant_term(newton_acc_X[i, j]) * constant_term(pn1[i, j, k]), order)
                        Y_t_pn1[i, j, k] = Taylor1(constant_term(newton_acc_Y[i, j]) * constant_term(pn1[i, j, k]), order)
                        Z_t_pn1[i, j, k] = Taylor1(constant_term(newton_acc_Z[i, j]) * constant_term(pn1[i, j, k]), order)
                        pNX_t_pn3[i, j, k] = Taylor1(constant_term(postNewtonX[i, k]) * constant_term(pn3[i, j]), order)
                        pNY_t_pn3[i, j, k] = Taylor1(constant_term(postNewtonY[i, k]) * constant_term(pn3[i, j]), order)
                        pNZ_t_pn3[i, j, k] = Taylor1(constant_term(postNewtonZ[i, k]) * constant_term(pn3[i, j]), order)
                    end
                end
            end
        for j = _1_to_N
            for i = _1_to_N
                if i == j
                else
                    pntempX[j, k] = Taylor1(constant_term(pntempX[j, k]) + constant_term(X_t_pn1[i, j, k] + (U_t_pn2[i, j] + pNX_t_pn3[i, j, k])), order)
                    pntempY[j, k] = Taylor1(constant_term(pntempY[j, k]) + constant_term(Y_t_pn1[i, j, k] + (V_t_pn2[i, j] + pNY_t_pn3[i, j, k])), order)
                    pntempZ[j, k] = Taylor1(constant_term(pntempZ[j, k]) + constant_term(Z_t_pn1[i, j, k] + (W_t_pn2[i, j] + pNZ_t_pn3[i, j, k])), order)
                end
            end
            postNewtonX[j, k + 1] = Taylor1(constant_term(pntempX[j, k]) * constant_term(c_m2), order)
            postNewtonY[j, k + 1] = Taylor1(constant_term(pntempY[j, k]) * constant_term(c_m2), order)
            postNewtonZ[j, k + 1] = Taylor1(constant_term(pntempZ[j, k]) * constant_term(c_m2), order)
        end
    end
    #= In[16]:492 =# Threads.@threads for i = _1_to_N
            dq[3 * (N + i) - 2] = Taylor1(constant_term(postNewtonX[i, postnewton_iter + 1]) + constant_term(accX[i]), order)
            dq[3 * (N + i) - 1] = Taylor1(constant_term(postNewtonY[i, postnewton_iter + 1]) + constant_term(accY[i]), order)
            dq[3 * (N + i)] = Taylor1(constant_term(postNewtonZ[i, postnewton_iter + 1]) + constant_term(accZ[i]), order)
        end
    for __idx = eachindex(q)
        (q[__idx]).coeffs[2] = (dq[__idx]).coeffs[1]
    end
    for ord = 1:order - 1
        ordnext = ord + 1
        #= In[16]:83 =# Threads.@threads for j = _1_to_N
                for i = _1_to_N
                    if i == j
                    else
                        TaylorSeries.identity!(pn2[i, j], zero_q_1, ord)
                        TaylorSeries.identity!(pn3[i, j], zero_q_1, ord)
                        TaylorSeries.identity!(_4ϕj[i, j], zero_q_1, ord)
                        TaylorSeries.identity!(ϕi_plus_4ϕj[i, j], zero_q_1, ord)
                        TaylorSeries.identity!(sj2_plus_2si2_minus_4vivj[i, j], zero_q_1, ord)
                        TaylorSeries.identity!(ϕs_and_vs[i, j], zero_q_1, ord)
                        TaylorSeries.identity!(UU[i, j], zero_q_1, ord)
                        TaylorSeries.identity!(VV[i, j], zero_q_1, ord)
                        TaylorSeries.identity!(WW[i, j], zero_q_1, ord)
                        TaylorSeries.identity!(U_t_pn2[i, j], zero_q_1, ord)
                        TaylorSeries.identity!(V_t_pn2[i, j], zero_q_1, ord)
                        TaylorSeries.identity!(W_t_pn2[i, j], zero_q_1, ord)
                        TaylorSeries.identity!(vi_dot_vj[i, j], zero_q_1, ord)
                        TaylorSeries.identity!(newton_acc_X[i, j], zero_q_1, ord)
                        TaylorSeries.identity!(newton_acc_Y[i, j], zero_q_1, ord)
                        TaylorSeries.identity!(newton_acc_Z[i, j], zero_q_1, ord)
                        TaylorSeries.identity!(pn1t1_7[i, j], zero_q_1, ord)
                    end
                end
            end
        for j = _1_to_N
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
        #= In[16]:200 =# Threads.@threads for j = j2_body_index
                for i = _1_to_N
                    if i == j
                    else
                        TaylorSeries.identity!(t31[i, j], zero_q_1, ord)
                        TaylorSeries.identity!(t32[i, j], zero_q_1, ord)
                        TaylorSeries.identity!(t33[i, j], zero_q_1, ord)
                        TaylorSeries.identity!(F_J2_x[i, j], zero_q_1, ord)
                        TaylorSeries.identity!(F_J2_y[i, j], zero_q_1, ord)
                        TaylorSeries.identity!(F_J2_z[i, j], zero_q_1, ord)
                        TaylorSeries.identity!(F_J2_x1[i, j], zero_q_1, ord)
                        TaylorSeries.identity!(F_J2_y1[i, j], zero_q_1, ord)
                        TaylorSeries.identity!(F_J2_z1[i, j], zero_q_1, ord)
                        TaylorSeries.identity!(F_J2_x2[i, j], zero_q_1, ord)
                        TaylorSeries.identity!(F_J2_y2[i, j], zero_q_1, ord)
                        TaylorSeries.identity!(F_J2_z2[i, j], zero_q_1, ord)
                        TaylorSeries.identity!(sin_ϕ[i, j], zero_q_1, ord)
                        TaylorSeries.identity!(sin2_ϕ[i, j], zero_q_1, ord)
                        TaylorSeries.identity!(sin3_ϕ[i, j], zero_q_1, ord)
                        TaylorSeries.identity!(sin4_ϕ[i, j], zero_q_1, ord)
                        TaylorSeries.identity!(ϕ[i, j], zero_q_1, ord)
                        TaylorSeries.identity!(cos_ϕ[i, j], zero_q_1, ord)
                        TaylorSeries.identity!(P_2_sin_ϕ[i, j], zero_q_1, ord)
                        TaylorSeries.identity!(∂P_2_sin_ϕ[i, j], zero_q_1, ord)
                        TaylorSeries.identity!(P_3_sin_ϕ[i, j], zero_q_1, ord)
                        TaylorSeries.identity!(∂P_3_sin_ϕ[i, j], zero_q_1, ord)
                        TaylorSeries.identity!(m_c_ϕ_∂P_2[i, j], zero_q_1, ord)
                        TaylorSeries.identity!(m_c_ϕ_∂P_3[i, j], zero_q_1, ord)
                        TaylorSeries.identity!(Λ2j_div_r4[i, j], zero_q_1, ord)
                        TaylorSeries.identity!(Λ3j_div_r5[i, j], zero_q_1, ord)
                        TaylorSeries.identity!(F_J_ξ[i, j], zero_q_1, ord)
                        TaylorSeries.identity!(F_J_η[i, j], zero_q_1, ord)
                        TaylorSeries.identity!(F_J_ζ[i, j], zero_q_1, ord)
                        TaylorSeries.identity!(F_J2_ξ[i, j], zero_q_1, ord)
                        TaylorSeries.identity!(F_J2_η[i, j], zero_q_1, ord)
                        TaylorSeries.identity!(F_J2_ζ[i, j], zero_q_1, ord)
                        TaylorSeries.identity!(F_J3_ξ[i, j], zero_q_1, ord)
                        TaylorSeries.identity!(F_J3_η[i, j], zero_q_1, ord)
                        TaylorSeries.identity!(F_J3_ζ[i, j], zero_q_1, ord)
                        TaylorSeries.identity!(ξx[i, j], zero_q_1, ord)
                        TaylorSeries.identity!(ξy[i, j], zero_q_1, ord)
                        TaylorSeries.identity!(ξz[i, j], zero_q_1, ord)
                        TaylorSeries.identity!(ηx[i, j], zero_q_1, ord)
                        TaylorSeries.identity!(ηy[i, j], zero_q_1, ord)
                        TaylorSeries.identity!(ηz[i, j], zero_q_1, ord)
                        TaylorSeries.identity!(ηx1[i, j], zero_q_1, ord)
                        TaylorSeries.identity!(ηy1[i, j], zero_q_1, ord)
                        TaylorSeries.identity!(ηz1[i, j], zero_q_1, ord)
                        TaylorSeries.identity!(ηx2[i, j], zero_q_1, ord)
                        TaylorSeries.identity!(ηy2[i, j], zero_q_1, ord)
                        TaylorSeries.identity!(ηz2[i, j], zero_q_1, ord)
                        TaylorSeries.identity!(ζx[i, j], zero_q_1, ord)
                        TaylorSeries.identity!(ζy[i, j], zero_q_1, ord)
                        TaylorSeries.identity!(ζz[i, j], zero_q_1, ord)
                        TaylorSeries.identity!(ζx1[i, j], zero_q_1, ord)
                        TaylorSeries.identity!(ζy1[i, j], zero_q_1, ord)
                        TaylorSeries.identity!(ζz1[i, j], zero_q_1, ord)
                        TaylorSeries.identity!(ζx2[i, j], zero_q_1, ord)
                        TaylorSeries.identity!(ζy2[i, j], zero_q_1, ord)
                        TaylorSeries.identity!(ζz2[i, j], zero_q_1, ord)
                    end
                end
            end
        #= In[16]:265 =# Threads.@threads for j = _1_to_N
                for i = _1_to_N
                    if i == j
                    else
                        TaylorSeries.subst!(X[i, j], q[3i - 2], q[3j - 2], ord)
                        TaylorSeries.subst!(Y[i, j], q[3i - 1], q[3j - 1], ord)
                        TaylorSeries.subst!(Z[i, j], q[3i], q[3j], ord)
                        TaylorSeries.subst!(U[i, j], dq[3i - 2], dq[3j - 2], ord)
                        TaylorSeries.subst!(V[i, j], dq[3i - 1], dq[3j - 1], ord)
                        TaylorSeries.subst!(W[i, j], dq[3i], dq[3j], ord)
                        TaylorSeries.mul!(tmp1530[3j - 2], 4, dq[3j - 2], ord)
                        TaylorSeries.mul!(tmp1532[3i - 2], 3, dq[3i - 2], ord)
                        TaylorSeries.subst!(_4U_m_3X[i, j], tmp1530[3j - 2], tmp1532[3i - 2], ord)
                        TaylorSeries.mul!(tmp1535[3j - 1], 4, dq[3j - 1], ord)
                        TaylorSeries.mul!(tmp1537[3i - 1], 3, dq[3i - 1], ord)
                        TaylorSeries.subst!(_4V_m_3Y[i, j], tmp1535[3j - 1], tmp1537[3i - 1], ord)
                        TaylorSeries.mul!(tmp1540[3j], 4, dq[3j], ord)
                        TaylorSeries.mul!(tmp1542[3i], 3, dq[3i], ord)
                        TaylorSeries.subst!(_4W_m_3Z[i, j], tmp1540[3j], tmp1542[3i], ord)
                        TaylorSeries.mul!(pn2x[i, j], X[i, j], _4U_m_3X[i, j], ord)
                        TaylorSeries.mul!(pn2y[i, j], Y[i, j], _4V_m_3Y[i, j], ord)
                        TaylorSeries.mul!(pn2z[i, j], Z[i, j], _4W_m_3Z[i, j], ord)
                        TaylorSeries.mul!(UU[i, j], dq[3i - 2], dq[3j - 2], ord)
                        TaylorSeries.mul!(VV[i, j], dq[3i - 1], dq[3j - 1], ord)
                        TaylorSeries.mul!(WW[i, j], dq[3i], dq[3j], ord)
                        TaylorSeries.add!(tmp1550[i, j], UU[i, j], VV[i, j], ord)
                        TaylorSeries.add!(vi_dot_vj[i, j], tmp1550[i, j], WW[i, j], ord)
                        TaylorSeries.pow!(tmp1553[i, j], X[i, j], 2, ord)
                        TaylorSeries.pow!(tmp1555[i, j], Y[i, j], 2, ord)
                        TaylorSeries.add!(tmp1556[i, j], tmp1553[i, j], tmp1555[i, j], ord)
                        TaylorSeries.pow!(tmp1558[i, j], Z[i, j], 2, ord)
                        TaylorSeries.add!(r_p2[i, j], tmp1556[i, j], tmp1558[i, j], ord)
                        TaylorSeries.sqrt!(r_p1d2[i, j], r_p2[i, j], ord)
                        TaylorSeries.pow!(r_p3d2[i, j], r_p2[i, j], 1.5, ord)
                        TaylorSeries.pow!(r_p7d2[i, j], r_p2[i, j], 3.5, ord)
                        TaylorSeries.div!(newtonianCoeff[i, j], μ[i], r_p3d2[i, j], ord)
                        TaylorSeries.add!(tmp1566[i, j], pn2x[i, j], pn2y[i, j], ord)
                        TaylorSeries.add!(tmp1567[i, j], tmp1566[i, j], pn2z[i, j], ord)
                        TaylorSeries.mul!(pn2[i, j], newtonianCoeff[i, j], tmp1567[i, j], ord)
                        TaylorSeries.mul!(newton_acc_X[i, j], X[i, j], newtonianCoeff[i, j], ord)
                        TaylorSeries.mul!(newton_acc_Y[i, j], Y[i, j], newtonianCoeff[i, j], ord)
                        TaylorSeries.mul!(newton_acc_Z[i, j], Z[i, j], newtonianCoeff[i, j], ord)
                        TaylorSeries.div!(newtonian1b_Potential[i, j], μ[i], r_p1d2[i, j], ord)
                        TaylorSeries.mul!(pn3[i, j], 3.5, newtonian1b_Potential[i, j], ord)
                        TaylorSeries.mul!(U_t_pn2[i, j], pn2[i, j], U[i, j], ord)
                        TaylorSeries.mul!(V_t_pn2[i, j], pn2[i, j], V[i, j], ord)
                        TaylorSeries.mul!(W_t_pn2[i, j], pn2[i, j], W[i, j], ord)
                        if UJ_interaction[i, j]
                            TaylorSeries.mul!(t31[i, j], X[i, j], M_[1, 3, j], ord)
                            TaylorSeries.mul!(t32[i, j], Y[i, j], M_[2, 3, j], ord)
                            TaylorSeries.mul!(t33[i, j], Z[i, j], M_[3, 3, j], ord)
                            TaylorSeries.add!(tmp1581[i, j], t31[i, j], t32[i, j], ord)
                            TaylorSeries.add!(r_sin_ϕ[i, j], tmp1581[i, j], t33[i, j], ord)
                            TaylorSeries.div!(sin_ϕ[i, j], r_sin_ϕ[i, j], r_p1d2[i, j], ord)
                            TaylorSeries.asin!(ϕ[i, j], sin_ϕ[i, j], tmp1730[i, j], ord)
                            TaylorSeries.sincos!(tmp1731[i, j], cos_ϕ[i, j], ϕ[i, j], ord)
                            TaylorSeries.pow!(sin2_ϕ[i, j], sin_ϕ[i, j], 2, ord)
                            TaylorSeries.pow!(sin3_ϕ[i, j], sin_ϕ[i, j], 3, ord)
                            TaylorSeries.mul!(tmp1591[i, j], 1.5, sin2_ϕ[i, j], ord)
                            TaylorSeries.subst!(P_2_sin_ϕ[i, j], tmp1591[i, j], 0.5, ord)
                            TaylorSeries.mul!(∂P_2_sin_ϕ[i, j], 3, sin_ϕ[i, j], ord)
                            TaylorSeries.mul!(tmp1597[i, j], -1.5, sin_ϕ[i, j], ord)
                            TaylorSeries.mul!(tmp1599[i, j], 2.5, sin3_ϕ[i, j], ord)
                            TaylorSeries.add!(P_3_sin_ϕ[i, j], tmp1597[i, j], tmp1599[i, j], ord)
                            TaylorSeries.mul!(tmp1603[i, j], 7.5, sin2_ϕ[i, j], ord)
                            TaylorSeries.add!(∂P_3_sin_ϕ[i, j], -1.5, tmp1603[i, j], ord)
                            TaylorSeries.subst!(tmp1605[j], Λ2[j], ord)
                            TaylorSeries.pow!(tmp1607[i, j], r_p2[i, j], 2, ord)
                            TaylorSeries.div!(Λ2j_div_r4[i, j], tmp1605[j], tmp1607[i, j], ord)
                            TaylorSeries.subst!(tmp1609[j], Λ3[j], ord)
                            TaylorSeries.pow!(tmp1611[i, j], r_p1d2[i, j], 5, ord)
                            TaylorSeries.div!(Λ3j_div_r5[i, j], tmp1609[j], tmp1611[i, j], ord)
                            TaylorSeries.subst!(tmp1613[i, j], cos_ϕ[i, j], ord)
                            TaylorSeries.mul!(m_c_ϕ_∂P_2[i, j], tmp1613[i, j], ∂P_2_sin_ϕ[i, j], ord)
                            TaylorSeries.subst!(tmp1615[i, j], cos_ϕ[i, j], ord)
                            TaylorSeries.mul!(m_c_ϕ_∂P_3[i, j], tmp1615[i, j], ∂P_3_sin_ϕ[i, j], ord)
                            TaylorSeries.mul!(tmp1618[i, j], Λ2j_div_r4[i, j], 3, ord)
                            TaylorSeries.mul!(F_J2_ξ[i, j], tmp1618[i, j], P_2_sin_ϕ[i, j], ord)
                            TaylorSeries.mul!(F_J2_ζ[i, j], Λ2j_div_r4[i, j], m_c_ϕ_∂P_2[i, j], ord)
                            TaylorSeries.mul!(tmp1622[i, j], Λ3j_div_r5[i, j], 4, ord)
                            TaylorSeries.mul!(F_J3_ξ[i, j], tmp1622[i, j], P_3_sin_ϕ[i, j], ord)
                            TaylorSeries.mul!(F_J3_ζ[i, j], Λ3j_div_r5[i, j], m_c_ϕ_∂P_3[i, j], ord)
                            TaylorSeries.add!(F_J_ξ[i, j], F_J2_ξ[i, j], F_J3_ξ[i, j], ord)
                            TaylorSeries.add!(F_J_ζ[i, j], F_J2_ζ[i, j], F_J3_ζ[i, j], ord)
                            TaylorSeries.div!(ξx[i, j], X[i, j], r_p1d2[i, j], ord)
                            TaylorSeries.div!(ξy[i, j], Y[i, j], r_p1d2[i, j], ord)
                            TaylorSeries.div!(ξz[i, j], Z[i, j], r_p1d2[i, j], ord)
                            TaylorSeries.mul!(ηx1[i, j], M_[2, 3, j], ξz[i, j], ord)
                            TaylorSeries.mul!(ηy1[i, j], M_[3, 3, j], ξx[i, j], ord)
                            TaylorSeries.mul!(ηz1[i, j], M_[1, 3, j], ξy[i, j], ord)
                            TaylorSeries.mul!(ηx2[i, j], M_[3, 3, j], ξy[i, j], ord)
                            TaylorSeries.mul!(ηy2[i, j], M_[1, 3, j], ξz[i, j], ord)
                            TaylorSeries.mul!(ηz2[i, j], M_[2, 3, j], ξx[i, j], ord)
                            TaylorSeries.subst!(ηx[i, j], ηx1[i, j], ηx2[i, j], ord)
                            TaylorSeries.subst!(ηy[i, j], ηy1[i, j], ηy2[i, j], ord)
                            TaylorSeries.subst!(ηz[i, j], ηz1[i, j], ηz2[i, j], ord)
                            TaylorSeries.mul!(ζx1[i, j], ξy[i, j], ηz[i, j], ord)
                            TaylorSeries.mul!(ζy1[i, j], ξz[i, j], ηx[i, j], ord)
                            TaylorSeries.mul!(ζz1[i, j], ξx[i, j], ηy[i, j], ord)
                            TaylorSeries.mul!(ζx2[i, j], ξz[i, j], ηy[i, j], ord)
                            TaylorSeries.mul!(ζy2[i, j], ξx[i, j], ηz[i, j], ord)
                            TaylorSeries.mul!(ζz2[i, j], ξy[i, j], ηx[i, j], ord)
                            TaylorSeries.subst!(ζx[i, j], ζx1[i, j], ζx2[i, j], ord)
                            TaylorSeries.subst!(ζy[i, j], ζy1[i, j], ζy2[i, j], ord)
                            TaylorSeries.subst!(ζz[i, j], ζz1[i, j], ζz2[i, j], ord)
                            TaylorSeries.mul!(F_J2_x1[i, j], F_J_ξ[i, j], ξx[i, j], ord)
                            TaylorSeries.mul!(F_J2_y1[i, j], F_J_ξ[i, j], ξy[i, j], ord)
                            TaylorSeries.mul!(F_J2_z1[i, j], F_J_ξ[i, j], ξz[i, j], ord)
                            TaylorSeries.mul!(F_J2_x2[i, j], F_J_ζ[i, j], ζx[i, j], ord)
                            TaylorSeries.mul!(F_J2_y2[i, j], F_J_ζ[i, j], ζy[i, j], ord)
                            TaylorSeries.mul!(F_J2_z2[i, j], F_J_ζ[i, j], ζz[i, j], ord)
                            TaylorSeries.add!(F_J2_x[i, j], F_J2_x1[i, j], F_J2_x2[i, j], ord)
                            TaylorSeries.add!(F_J2_y[i, j], F_J2_y1[i, j], F_J2_y2[i, j], ord)
                            TaylorSeries.add!(F_J2_z[i, j], F_J2_z1[i, j], F_J2_z2[i, j], ord)
                        end
                    end
                end
                TaylorSeries.pow!(tmp1658[3j - 2], dq[3j - 2], 2, ord)
                TaylorSeries.pow!(tmp1660[3j - 1], dq[3j - 1], 2, ord)
                TaylorSeries.add!(tmp1661[3j - 2], tmp1658[3j - 2], tmp1660[3j - 1], ord)
                TaylorSeries.pow!(tmp1663[3j], dq[3j], 2, ord)
                TaylorSeries.add!(v2[j], tmp1661[3j - 2], tmp1663[3j], ord)
            end
        for j = _1_to_N
            for i = _1_to_N
                if i == j
                else
                    TaylorSeries.add!(newtonianNb_Potential[j], newtonianNb_Potential[j], newtonian1b_Potential[i, j], ord)
                    TaylorSeries.add!(newtonX[j], newtonX[j], newton_acc_X[i, j], ord)
                    TaylorSeries.add!(newtonY[j], newtonY[j], newton_acc_Y[i, j], ord)
                    TaylorSeries.add!(newtonZ[j], newtonZ[j], newton_acc_Z[i, j], ord)
                    if UJ_interaction[i, j]
                        TaylorSeries.add!(accX[j], accX[j], μ[i] * F_J2_x[i, j], ord)
                        TaylorSeries.add!(accY[j], accY[j], μ[i] * F_J2_y[i, j], ord)
                        TaylorSeries.add!(accZ[j], accZ[j], μ[i] * F_J2_z[i, j], ord)
                        TaylorSeries.subst!(accX[i], accX[i], μ[j] * F_J2_x[i, j], ord)
                        TaylorSeries.subst!(accY[i], accY[i], μ[j] * F_J2_y[i, j], ord)
                        TaylorSeries.subst!(accZ[i], accZ[i], μ[j] * F_J2_z[i, j], ord)
                    end
                end
            end
        end
        #= In[16]:411 =# Threads.@threads for j = _1_to_N
                for i = _1_to_N
                    if i == j
                    else
                        TaylorSeries.mul!(_4ϕj[i, j], 4, newtonianNb_Potential[j], ord)
                        TaylorSeries.add!(ϕi_plus_4ϕj[i, j], newtonianNb_Potential[i], _4ϕj[i, j], ord)
                        TaylorSeries.mul!(tmp1685[i], 2, v2[i], ord)
                        TaylorSeries.add!(tmp1686[j], v2[j], tmp1685[i], ord)
                        TaylorSeries.mul!(tmp1688[i, j], 4, vi_dot_vj[i, j], ord)
                        TaylorSeries.subst!(sj2_plus_2si2_minus_4vivj[i, j], tmp1686[j], tmp1688[i, j], ord)
                        TaylorSeries.subst!(ϕs_and_vs[i, j], sj2_plus_2si2_minus_4vivj[i, j], ϕi_plus_4ϕj[i, j], ord)
                        TaylorSeries.mul!(Xij_t_Ui[i, j], X[i, j], dq[3i - 2], ord)
                        TaylorSeries.mul!(Yij_t_Vi[i, j], Y[i, j], dq[3i - 1], ord)
                        TaylorSeries.mul!(Zij_t_Wi[i, j], Z[i, j], dq[3i], ord)
                        TaylorSeries.add!(tmp1694[i, j], Xij_t_Ui[i, j], Yij_t_Vi[i, j], ord)
                        TaylorSeries.add!(Rij_dot_Vi[i, j], tmp1694[i, j], Zij_t_Wi[i, j], ord)
                        TaylorSeries.pow!(tmp1697[i, j], Rij_dot_Vi[i, j], 2, ord)
                        TaylorSeries.div!(pn1t7[i, j], tmp1697[i, j], r_p2[i, j], ord)
                        TaylorSeries.mul!(tmp1700[i, j], 1.5, pn1t7[i, j], ord)
                        TaylorSeries.subst!(pn1t2_7[i, j], ϕs_and_vs[i, j], tmp1700[i, j], ord)
                        TaylorSeries.add!(pn1t1_7[i, j], c_p2, pn1t2_7[i, j], ord)
                        for k = Base.OneTo(postnewton_iter)
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
                for k = Base.OneTo(postnewton_iter)
                    TaylorSeries.identity!(pntempX[j, k], zero_q_1, ord)
                    TaylorSeries.identity!(pntempY[j, k], zero_q_1, ord)
                    TaylorSeries.identity!(pntempZ[j, k], zero_q_1, ord)
                end
            end
        for k = Base.OneTo(postnewton_iter)
            #= In[16]:455 =# Threads.@threads for j = _1_to_N
                    for i = _1_to_N
                        if i == j
                        else
                            TaylorSeries.mul!(pNX_t_X[i, j, k], postNewtonX[i, k], X[i, j], ord)
                            TaylorSeries.mul!(pNY_t_Y[i, j, k], postNewtonY[i, k], Y[i, j], ord)
                            TaylorSeries.mul!(pNZ_t_Z[i, j, k], postNewtonZ[i, k], Z[i, j], ord)
                            TaylorSeries.add!(tmp1706[i, j, k], pNX_t_X[i, j, k], pNY_t_Y[i, j, k], ord)
                            TaylorSeries.add!(tmp1707[i, j, k], tmp1706[i, j, k], pNZ_t_Z[i, j, k], ord)
                            TaylorSeries.add!(pn1[i, j, k], pn1t1_7[i, j], tmp1707[i, j, k], ord)
                            TaylorSeries.mul!(X_t_pn1[i, j, k], newton_acc_X[i, j], pn1[i, j, k], ord)
                            TaylorSeries.mul!(Y_t_pn1[i, j, k], newton_acc_Y[i, j], pn1[i, j, k], ord)
                            TaylorSeries.mul!(Z_t_pn1[i, j, k], newton_acc_Z[i, j], pn1[i, j, k], ord)
                            TaylorSeries.mul!(pNX_t_pn3[i, j, k], postNewtonX[i, k], pn3[i, j], ord)
                            TaylorSeries.mul!(pNY_t_pn3[i, j, k], postNewtonY[i, k], pn3[i, j], ord)
                            TaylorSeries.mul!(pNZ_t_pn3[i, j, k], postNewtonZ[i, k], pn3[i, j], ord)
                        end
                    end
                end
            for j = _1_to_N
                for i = _1_to_N
                    if i == j
                    else
                        TaylorSeries.add!(pntempX[j, k], pntempX[j, k], X_t_pn1[i, j, k] + (U_t_pn2[i, j] + pNX_t_pn3[i, j, k]), ord)
                        TaylorSeries.add!(pntempY[j, k], pntempY[j, k], Y_t_pn1[i, j, k] + (V_t_pn2[i, j] + pNY_t_pn3[i, j, k]), ord)
                        TaylorSeries.add!(pntempZ[j, k], pntempZ[j, k], Z_t_pn1[i, j, k] + (W_t_pn2[i, j] + pNZ_t_pn3[i, j, k]), ord)
                    end
                end
                TaylorSeries.mul!(postNewtonX[j, k + 1], pntempX[j, k], c_m2, ord)
                TaylorSeries.mul!(postNewtonY[j, k + 1], pntempY[j, k], c_m2, ord)
                TaylorSeries.mul!(postNewtonZ[j, k + 1], pntempZ[j, k], c_m2, ord)
            end
        end
        #= In[16]:492 =# Threads.@threads for i = _1_to_N
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
