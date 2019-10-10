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
@taylorize function NBP_pN_A_J23E_J23M_J2S!(dq, q, params, t)
    local S = eltype(q[1])
    local N = Int((length(q))/6) # number of bodies
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
    temp_accX_j = Array{Taylor1{S}}(undef, N, N)
    temp_accY_j = Array{Taylor1{S}}(undef, N, N)
    temp_accZ_j = Array{Taylor1{S}}(undef, N, N)
    temp_accX_i = Array{Taylor1{S}}(undef, N, N)
    temp_accY_i = Array{Taylor1{S}}(undef, N, N)
    temp_accZ_i = Array{Taylor1{S}}(undef, N, N)
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
                temp_001 = newtonX[j] + (X[i,j]*newtonianCoeff[i,j])
                newtonX[j] = temp_001
                temp_002 = newtonY[j] + (Y[i,j]*newtonianCoeff[i,j])
                newtonY[j] = temp_002
                temp_003 = newtonZ[j] + (Z[i,j]*newtonianCoeff[i,j])
                newtonZ[j] = temp_003
                temp_004 = newtonianNb_Potential[j] + newtonian1b_Potential[i, j]
                newtonianNb_Potential[j] = temp_004
                if UJ_interaction[i,j]
                    # # add result to total acceleration on upon j-th body figure due to i-th point mass
                    # @show "acc",j,"+μ",i,"Λ2",j
                    temp_accX_j[i,j] = accX[j] + (μ[i]*F_J2_x[i,j])
                    accX[j] = temp_accX_j[i,j]
                    temp_accY_j[i,j] = accY[j] + (μ[i]*F_J2_y[i,j])
                    accY[j] = temp_accY_j[i,j]
                    temp_accZ_j[i,j] = accZ[j] + (μ[i]*F_J2_z[i,j])
                    accZ[j] = temp_accZ_j[i,j]

                    # # reaction force on i-th body
                    # @show "acc",i,"-μ",j,"Λ2",j
                    temp_accX_i[i,j] = accX[i] - (μ[j]*F_J2_x[i,j])
                    accX[i] = temp_accX_i[i,j]
                    temp_accY_i[i,j] = accY[i] - (μ[j]*F_J2_y[i,j])
                    accY[i] = temp_accY_i[i,j]
                    temp_accZ_i[i,j] = accZ[i] - (μ[j]*F_J2_z[i,j])
                    accZ[i] = temp_accZ_i[i,j]
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
                    pn1[i,j,k] = (  pn1t1_7[i,j]  +  0.5*( (pNX_t_X[i,j,k]+pNY_t_Y[i,j,k]) + pNZ_t_Z[i,j,k] )  )

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
                    termpnx = ( X_t_pn1[i,j,k] + (U_t_pn2[i,j]+pNX_t_pn3[i,j,k]) )
                    sumpnx = pntempX[j,k] + termpnx
                    pntempX[j,k] = sumpnx
                    termpny = ( Y_t_pn1[i,j,k] + (V_t_pn2[i,j]+pNY_t_pn3[i,j,k]) )
                    sumpny = pntempY[j,k] + termpny
                    pntempY[j,k] = sumpny
                    termpnz = ( Z_t_pn1[i,j,k] + (W_t_pn2[i,j]+pNZ_t_pn3[i,j,k]) )
                    sumpnz = pntempZ[j,k] + termpnz
                    pntempZ[j,k] = sumpnz
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

@taylorize function NBP_pN_A_J23E_J23M_J2S_threads!(dq, q, params, t)
    local S = eltype(q[1])
    local N = Int((length(q))/6) # number of bodies
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
    temp_accX_j = Array{Taylor1{S}}(undef, N, N)
    temp_accY_j = Array{Taylor1{S}}(undef, N, N)
    temp_accZ_j = Array{Taylor1{S}}(undef, N, N)
    temp_accX_i = Array{Taylor1{S}}(undef, N, N)
    temp_accY_i = Array{Taylor1{S}}(undef, N, N)
    temp_accZ_i = Array{Taylor1{S}}(undef, N, N)
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

    Threads.@threads for j in _1_to_N
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
                temp_001 = newtonX[j] + (X[i,j]*newtonianCoeff[i,j])
                newtonX[j] = temp_001
                temp_002 = newtonY[j] + (Y[i,j]*newtonianCoeff[i,j])
                newtonY[j] = temp_002
                temp_003 = newtonZ[j] + (Z[i,j]*newtonianCoeff[i,j])
                newtonZ[j] = temp_003
                temp_004 = newtonianNb_Potential[j] + newtonian1b_Potential[i, j]
                newtonianNb_Potential[j] = temp_004
                if UJ_interaction[i,j]
                    # # add result to total acceleration on upon j-th body figure due to i-th point mass
                    # @show "acc",j,"+μ",i,"Λ2",j
                    temp_accX_j[i,j] = accX[j] + (μ[i]*F_J2_x[i,j])
                    accX[j] = temp_accX_j[i,j]
                    temp_accY_j[i,j] = accY[j] + (μ[i]*F_J2_y[i,j])
                    accY[j] = temp_accY_j[i,j]
                    temp_accZ_j[i,j] = accZ[j] + (μ[i]*F_J2_z[i,j])
                    accZ[j] = temp_accZ_j[i,j]

                    # # reaction force on i-th body
                    # @show "acc",i,"-μ",j,"Λ2",j
                    temp_accX_i[i,j] = accX[i] - (μ[j]*F_J2_x[i,j])
                    accX[i] = temp_accX_i[i,j]
                    temp_accY_i[i,j] = accY[i] - (μ[j]*F_J2_y[i,j])
                    accY[i] = temp_accY_i[i,j]
                    temp_accZ_i[i,j] = accZ[i] - (μ[j]*F_J2_z[i,j])
                    accZ[i] = temp_accZ_i[i,j]
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
                    pn1[i,j,k] = (  pn1t1_7[i,j]  +  0.5*( (pNX_t_X[i,j,k]+pNY_t_Y[i,j,k]) + pNZ_t_Z[i,j,k] )  )

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
                    termpnx = ( X_t_pn1[i,j,k] + (U_t_pn2[i,j]+pNX_t_pn3[i,j,k]) )
                    sumpnx = pntempX[j,k] + termpnx
                    pntempX[j,k] = sumpnx
                    termpny = ( Y_t_pn1[i,j,k] + (V_t_pn2[i,j]+pNY_t_pn3[i,j,k]) )
                    sumpny = pntempY[j,k] + termpny
                    pntempY[j,k] = sumpny
                    termpnz = ( Z_t_pn1[i,j,k] + (W_t_pn2[i,j]+pNZ_t_pn3[i,j,k]) )
                    sumpnz = pntempZ[j,k] + termpnz
                    pntempZ[j,k] = sumpnz
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