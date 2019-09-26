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

    local succ_approx_iter = 1 # number of iterations of post-Newtonian subroutine
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

    postNewtonX = Array{Taylor1{S}}(undef, N)
    postNewtonY = Array{Taylor1{S}}(undef, N)
    postNewtonZ = Array{Taylor1{S}}(undef, N)

    newtonianNb_Potential = Array{Taylor1{S}}(undef, N)
    newtonian1b_Potential = Array{Taylor1{S}}(undef, N, N)
    newtonianCoeff = Array{Taylor1{S}}(undef, N, N)

    pntempX = Array{Taylor1{S}}(undef, N)
    pntempY = Array{Taylor1{S}}(undef, N)
    pntempZ = Array{Taylor1{S}}(undef, N)

    pn1 = Array{Taylor1{S}}(undef, N, N)
    v2 = Array{Taylor1{S}}(undef, N)
    vi_dot_vj = Array{Taylor1{S}}(undef, N, N)
    pn2 = Array{Taylor1{S}}(undef, N, N)
    pn3 = Array{Taylor1{S}}(undef, N, N)

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
                temp_accX_j[i,j] = zero_q_1
                temp_accY_j[i,j] = zero_q_1
                temp_accZ_j[i,j] = zero_q_1
                temp_accX_i[i,j] = zero_q_1
                temp_accY_i[i,j] = zero_q_1
                temp_accZ_i[i,j] = zero_q_1
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

                temp_001 = newtonX[j] + (X[i,j]*newtonianCoeff[i,j])
                newtonX[j] = temp_001
                temp_002 = newtonY[j] + (Y[i,j]*newtonianCoeff[i,j])
                newtonY[j] = temp_002
                temp_003 = newtonZ[j] + (Z[i,j]*newtonianCoeff[i,j])
                newtonZ[j] = temp_003

                newtonian1b_Potential[i,j] = μ[i]/r_p1d2[i, j]

                temp_004 = newtonianNb_Potential[j] + newtonian1b_Potential[i, j]
                newtonianNb_Potential[j] = temp_004

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
            end #if i != j
        end #for, i
        v2[j] = ( (dq[3j-2]^2)+(dq[3j-1]^2) ) + (dq[3j]^2)
    end #for, j

    for j in _1_to_N
        postNewtonX[j] = newtonX[j]
        postNewtonY[j] = newtonY[j]
        postNewtonZ[j] = newtonZ[j]
    end

    for k in Base.OneTo(succ_approx_iter)
        for j in _1_to_N
            pntempX[j] = zero_q_1
            pntempY[j] = zero_q_1
            pntempZ[j] = zero_q_1
        end
        for j in _1_to_N
            for i in _1_to_N
                # i == j && continue
                if i == j
                else
                    #post-Newtonian corrections to gravitational acceleration
                    #Moyer, 1971, page 7 eq. 35
                    temp_005a = newtonianNb_Potential[i]+(4newtonianNb_Potential[j])
                    temp_005b = (2v2[i])-(4vi_dot_vj[i,j])
                    temp_005c = v2[j]+temp_005b
                    temp_005 = temp_005c-temp_005a
                    temp_006a = X[i,j]*dq[3i-2]
                    temp_006b = Y[i,j]*dq[3i-1]
                    temp_006c = Z[i,j]*dq[3i]
                    temp_006d = ( temp_006a+temp_006b ) + temp_006c
                    # the expression below inside the (...)^2 should have a minus sign in front of the numerator,
                    # but upon squaring it is eliminated, so at the end of the day, it is irrelevant ;)
                    temp_006e = (temp_006d^2)/r_p2[i,j]
                    temp_006 = temp_005-(1.5temp_006e)
                    temp_007a = X[i,j]*postNewtonX[i]
                    temp_007b = Y[i,j]*postNewtonY[i]
                    temp_007c = Z[i,j]*postNewtonZ[i]
                    temp_007d = ( temp_007a+temp_007b ) + temp_007c
                    temp_007 = temp_006 + (0.5temp_007d)
                    temp_008 = c_p2+temp_007
                    pn1[i,j] = newtonianCoeff[i,j]*temp_008

                    temp_009 = X[i,j]*pn1[i,j]
                    temp_010 = Y[i,j]*pn1[i,j]
                    temp_011 = Z[i,j]*pn1[i,j]

                    pn3[i,j] = 3.5*newtonian1b_Potential[i,j]

                    temp_013a = pn2[i,j]*U[i,j]
                    temp_013b = pn3[i,j]*postNewtonX[i]
                    temp_013 = pntempX[j] + (temp_009 + (temp_013a+temp_013b))
                    pntempX[j] = temp_013

                    temp_014a = pn2[i,j]*V[i,j]
                    temp_014b = pn3[i,j]*postNewtonY[i]
                    temp_014 = pntempY[j] + (temp_010 + (temp_014a+temp_014b))
                    pntempY[j] = temp_014

                    temp_015a = pn2[i,j]*W[i,j]
                    temp_015b = pn3[i,j]*postNewtonZ[i]
                    temp_015 = pntempZ[j] + (temp_011 + (temp_015a+temp_015b))
                    pntempZ[j] = temp_015
                end
            end #for i
        end #for j
        for j in _1_to_N
            postNewtonX[j] = pntempX[j]*c_m2
            postNewtonY[j] = pntempY[j]*c_m2
            postNewtonZ[j] = pntempZ[j]*c_m2
        end
    end #for k in Base.OneTo(succ_approx_iter) # (post-Newtonian iterations)

    #fill accelerations (post-Newtonian and extended body accelerations)
    for i in _1_to_N
        dq[3(N+i)-2] = postNewtonX[i] + accX[i]
        dq[3(N+i)-1] = postNewtonY[i] + accY[i]
        dq[3(N+i)  ] = postNewtonZ[i] + accZ[i]
    end

    nothing
end

function TaylorIntegration.jetcoeffs!(::Val{NBP_pN_A_J23E_J23M_J2S!}, t::Taylor1{_T}, q::AbstractVector{Taylor1{_S}}, dq::AbstractVector{Taylor1{_S}}, params) where {_T <: Real, _S <: Number}
    order = t.order
    local S = eltype(q[1])
    local N = Int(length(q) / 6)
    local _1_to_N = Base.OneTo(N)
    local succ_approx_iter = 1
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
    postNewtonX = Array{Taylor1{S}}(undef, N)
    postNewtonY = Array{Taylor1{S}}(undef, N)
    postNewtonZ = Array{Taylor1{S}}(undef, N)
    newtonianNb_Potential = Array{Taylor1{S}}(undef, N)
    newtonian1b_Potential = Array{Taylor1{S}}(undef, N, N)
    pntempX = Array{Taylor1{S}}(undef, N)
    pntempY = Array{Taylor1{S}}(undef, N)
    pntempZ = Array{Taylor1{S}}(undef, N)
    pn1 = Array{Taylor1{S}}(undef, N, N)
    v2 = Array{Taylor1{S}}(undef, N)
    vi_dot_vj = Array{Taylor1{S}}(undef, N, N)
    pn2 = Array{Taylor1{S}}(undef, N, N)
    pn3 = Array{Taylor1{S}}(undef, N, N)
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
                temp_accX_j[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                temp_accY_j[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                temp_accZ_j[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                temp_accX_i[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                temp_accY_i[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                temp_accZ_i[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
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
    tmp793 = Array{Taylor1{_S}}(undef, size(dq))
    tmp793 .= Taylor1(zero(_S), order)
    tmp795 = Array{Taylor1{_S}}(undef, size(dq))
    tmp795 .= Taylor1(zero(_S), order)
    tmp798 = Array{Taylor1{_S}}(undef, size(dq))
    tmp798 .= Taylor1(zero(_S), order)
    tmp800 = Array{Taylor1{_S}}(undef, size(dq))
    tmp800 .= Taylor1(zero(_S), order)
    tmp803 = Array{Taylor1{_S}}(undef, size(dq))
    tmp803 .= Taylor1(zero(_S), order)
    tmp805 = Array{Taylor1{_S}}(undef, size(dq))
    tmp805 .= Taylor1(zero(_S), order)
    pn2x = Array{Taylor1{_S}}(undef, size(X))
    pn2x .= Taylor1(zero(_S), order)
    pn2y = Array{Taylor1{_S}}(undef, size(Y))
    pn2y .= Taylor1(zero(_S), order)
    pn2z = Array{Taylor1{_S}}(undef, size(Z))
    pn2z .= Taylor1(zero(_S), order)
    tmp813 = Array{Taylor1{_S}}(undef, size(UU))
    tmp813 .= Taylor1(zero(_S), order)
    tmp816 = Array{Taylor1{_S}}(undef, size(X))
    tmp816 .= Taylor1(zero(_S), order)
    tmp818 = Array{Taylor1{_S}}(undef, size(Y))
    tmp818 .= Taylor1(zero(_S), order)
    tmp819 = Array{Taylor1{_S}}(undef, size(tmp816))
    tmp819 .= Taylor1(zero(_S), order)
    tmp821 = Array{Taylor1{_S}}(undef, size(Z))
    tmp821 .= Taylor1(zero(_S), order)
    tmp829 = Array{Taylor1{_S}}(undef, size(pn2x))
    tmp829 .= Taylor1(zero(_S), order)
    tmp830 = Array{Taylor1{_S}}(undef, size(tmp829))
    tmp830 .= Taylor1(zero(_S), order)
    tmp832 = Array{Taylor1{_S}}(undef, size(X))
    tmp832 .= Taylor1(zero(_S), order)
    temp_001 = Array{Taylor1{_S}}(undef, size(tmp832))
    temp_001 .= Taylor1(zero(_S), order)
    tmp834 = Array{Taylor1{_S}}(undef, size(Y))
    tmp834 .= Taylor1(zero(_S), order)
    temp_002 = Array{Taylor1{_S}}(undef, size(tmp834))
    temp_002 .= Taylor1(zero(_S), order)
    tmp836 = Array{Taylor1{_S}}(undef, size(Z))
    tmp836 .= Taylor1(zero(_S), order)
    temp_003 = Array{Taylor1{_S}}(undef, size(tmp836))
    temp_003 .= Taylor1(zero(_S), order)
    temp_004 = Array{Taylor1{_S}}(undef, size(newtonian1b_Potential))
    temp_004 .= Taylor1(zero(_S), order)
    tmp843 = Array{Taylor1{_S}}(undef, size(t31))
    tmp843 .= Taylor1(zero(_S), order)
    tmp996 = Array{Taylor1{_S}}(undef, size(sin_ϕ))
    tmp996 .= Taylor1(zero(_S), order)
    tmp997 = Array{Taylor1{_S}}(undef, size(ϕ))
    tmp997 .= Taylor1(zero(_S), order)
    tmp853 = Array{Taylor1{_S}}(undef, size(sin2_ϕ))
    tmp853 .= Taylor1(zero(_S), order)
    tmp859 = Array{Taylor1{_S}}(undef, size(sin_ϕ))
    tmp859 .= Taylor1(zero(_S), order)
    tmp861 = Array{Taylor1{_S}}(undef, size(sin3_ϕ))
    tmp861 .= Taylor1(zero(_S), order)
    tmp865 = Array{Taylor1{_S}}(undef, size(sin2_ϕ))
    tmp865 .= Taylor1(zero(_S), order)
    tmp867 = Array{Taylor1{_S}}(undef, size(Λ2))
    tmp867 .= Taylor1(zero(_S), order)
    tmp869 = Array{Taylor1{_S}}(undef, size(r_p2))
    tmp869 .= Taylor1(zero(_S), order)
    tmp871 = Array{Taylor1{_S}}(undef, size(Λ3))
    tmp871 .= Taylor1(zero(_S), order)
    tmp873 = Array{Taylor1{_S}}(undef, size(r_p1d2))
    tmp873 .= Taylor1(zero(_S), order)
    tmp875 = Array{Taylor1{_S}}(undef, size(cos_ϕ))
    tmp875 .= Taylor1(zero(_S), order)
    tmp877 = Array{Taylor1{_S}}(undef, size(cos_ϕ))
    tmp877 .= Taylor1(zero(_S), order)
    tmp880 = Array{Taylor1{_S}}(undef, size(Λ2j_div_r4))
    tmp880 .= Taylor1(zero(_S), order)
    tmp884 = Array{Taylor1{_S}}(undef, size(Λ3j_div_r5))
    tmp884 .= Taylor1(zero(_S), order)
    tmp919 = Array{Taylor1{_S}}(undef, size(F_J2_x))
    tmp919 .= Taylor1(zero(_S), order)
    tmp921 = Array{Taylor1{_S}}(undef, size(F_J2_y))
    tmp921 .= Taylor1(zero(_S), order)
    tmp923 = Array{Taylor1{_S}}(undef, size(F_J2_z))
    tmp923 .= Taylor1(zero(_S), order)
    tmp925 = Array{Taylor1{_S}}(undef, size(F_J2_x))
    tmp925 .= Taylor1(zero(_S), order)
    tmp927 = Array{Taylor1{_S}}(undef, size(F_J2_y))
    tmp927 .= Taylor1(zero(_S), order)
    tmp929 = Array{Taylor1{_S}}(undef, size(F_J2_z))
    tmp929 .= Taylor1(zero(_S), order)
    tmp932 = Array{Taylor1{_S}}(undef, size(dq))
    tmp932 .= Taylor1(zero(_S), order)
    tmp934 = Array{Taylor1{_S}}(undef, size(dq))
    tmp934 .= Taylor1(zero(_S), order)
    tmp935 = Array{Taylor1{_S}}(undef, size(tmp932))
    tmp935 .= Taylor1(zero(_S), order)
    tmp937 = Array{Taylor1{_S}}(undef, size(dq))
    tmp937 .= Taylor1(zero(_S), order)
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
                tmp793[3j - 2] = Taylor1(constant_term(4) * constant_term(dq[3j - 2]), order)
                tmp795[3i - 2] = Taylor1(constant_term(3) * constant_term(dq[3i - 2]), order)
                _4U_m_3X[i, j] = Taylor1(constant_term(tmp793[3j - 2]) - constant_term(tmp795[3i - 2]), order)
                tmp798[3j - 1] = Taylor1(constant_term(4) * constant_term(dq[3j - 1]), order)
                tmp800[3i - 1] = Taylor1(constant_term(3) * constant_term(dq[3i - 1]), order)
                _4V_m_3Y[i, j] = Taylor1(constant_term(tmp798[3j - 1]) - constant_term(tmp800[3i - 1]), order)
                tmp803[3j] = Taylor1(constant_term(4) * constant_term(dq[3j]), order)
                tmp805[3i] = Taylor1(constant_term(3) * constant_term(dq[3i]), order)
                _4W_m_3Z[i, j] = Taylor1(constant_term(tmp803[3j]) - constant_term(tmp805[3i]), order)
                pn2x[i, j] = Taylor1(constant_term(X[i, j]) * constant_term(_4U_m_3X[i, j]), order)
                pn2y[i, j] = Taylor1(constant_term(Y[i, j]) * constant_term(_4V_m_3Y[i, j]), order)
                pn2z[i, j] = Taylor1(constant_term(Z[i, j]) * constant_term(_4W_m_3Z[i, j]), order)
                UU[i, j] = Taylor1(constant_term(dq[3i - 2]) * constant_term(dq[3j - 2]), order)
                VV[i, j] = Taylor1(constant_term(dq[3i - 1]) * constant_term(dq[3j - 1]), order)
                WW[i, j] = Taylor1(constant_term(dq[3i]) * constant_term(dq[3j]), order)
                tmp813[i, j] = Taylor1(constant_term(UU[i, j]) + constant_term(VV[i, j]), order)
                vi_dot_vj[i, j] = Taylor1(constant_term(tmp813[i, j]) + constant_term(WW[i, j]), order)
                tmp816[i, j] = Taylor1(constant_term(X[i, j]) ^ constant_term(2), order)
                tmp818[i, j] = Taylor1(constant_term(Y[i, j]) ^ constant_term(2), order)
                tmp819[i, j] = Taylor1(constant_term(tmp816[i, j]) + constant_term(tmp818[i, j]), order)
                tmp821[i, j] = Taylor1(constant_term(Z[i, j]) ^ constant_term(2), order)
                r_p2[i, j] = Taylor1(constant_term(tmp819[i, j]) + constant_term(tmp821[i, j]), order)
                r_p1d2[i, j] = Taylor1(sqrt(constant_term(r_p2[i, j])), order)
                r_p3d2[i, j] = Taylor1(constant_term(r_p2[i, j]) ^ constant_term(1.5), order)
                r_p7d2[i, j] = Taylor1(constant_term(r_p2[i, j]) ^ constant_term(3.5), order)
                newtonianCoeff[i, j] = Taylor1(constant_term(μ[i]) / constant_term(r_p3d2[i, j]), order)
                tmp829[i, j] = Taylor1(constant_term(pn2x[i, j]) + constant_term(pn2y[i, j]), order)
                tmp830[i, j] = Taylor1(constant_term(tmp829[i, j]) + constant_term(pn2z[i, j]), order)
                pn2[i, j] = Taylor1(constant_term(newtonianCoeff[i, j]) * constant_term(tmp830[i, j]), order)
                tmp832[i, j] = Taylor1(constant_term(X[i, j]) * constant_term(newtonianCoeff[i, j]), order)
                temp_001[i, j] = Taylor1(constant_term(newtonX[j]) + constant_term(tmp832[i, j]), order)
                newtonX[j] = Taylor1(identity(constant_term(temp_001[i, j])), order)
                tmp834[i, j] = Taylor1(constant_term(Y[i, j]) * constant_term(newtonianCoeff[i, j]), order)
                temp_002[i, j] = Taylor1(constant_term(newtonY[j]) + constant_term(tmp834[i, j]), order)
                newtonY[j] = Taylor1(identity(constant_term(temp_002[i, j])), order)
                tmp836[i, j] = Taylor1(constant_term(Z[i, j]) * constant_term(newtonianCoeff[i, j]), order)
                temp_003[i, j] = Taylor1(constant_term(newtonZ[j]) + constant_term(tmp836[i, j]), order)
                newtonZ[j] = Taylor1(identity(constant_term(temp_003[i, j])), order)
                newtonian1b_Potential[i, j] = Taylor1(constant_term(μ[i]) / constant_term(r_p1d2[i, j]), order)
                temp_004[i, j] = Taylor1(constant_term(newtonianNb_Potential[j]) + constant_term(newtonian1b_Potential[i, j]), order)
                newtonianNb_Potential[j] = Taylor1(identity(constant_term(temp_004[i, j])), order)
                if UJ_interaction[i, j]
                    t31[i, j] = Taylor1(constant_term(X[i, j]) * constant_term(M_[1, 3, j]), order)
                    t32[i, j] = Taylor1(constant_term(Y[i, j]) * constant_term(M_[2, 3, j]), order)
                    t33[i, j] = Taylor1(constant_term(Z[i, j]) * constant_term(M_[3, 3, j]), order)
                    tmp843[i, j] = Taylor1(constant_term(t31[i, j]) + constant_term(t32[i, j]), order)
                    r_sin_ϕ[i, j] = Taylor1(constant_term(tmp843[i, j]) + constant_term(t33[i, j]), order)
                    sin_ϕ[i, j] = Taylor1(constant_term(r_sin_ϕ[i, j]) / constant_term(r_p1d2[i, j]), order)
                    ϕ[i, j] = Taylor1(asin(constant_term(sin_ϕ[i, j])), order)
                    tmp996[i, j] = Taylor1(sqrt(1 - constant_term(sin_ϕ[i, j]) ^ 2), order)
                    cos_ϕ[i, j] = Taylor1(cos(constant_term(ϕ[i, j])), order)
                    tmp997[i, j] = Taylor1(sin(constant_term(ϕ[i, j])), order)
                    sin2_ϕ[i, j] = Taylor1(constant_term(sin_ϕ[i, j]) ^ constant_term(2), order)
                    sin3_ϕ[i, j] = Taylor1(constant_term(sin_ϕ[i, j]) ^ constant_term(3), order)
                    tmp853[i, j] = Taylor1(constant_term(1.5) * constant_term(sin2_ϕ[i, j]), order)
                    P_2_sin_ϕ[i, j] = Taylor1(constant_term(tmp853[i, j]) - constant_term(0.5), order)
                    ∂P_2_sin_ϕ[i, j] = Taylor1(constant_term(3) * constant_term(sin_ϕ[i, j]), order)
                    tmp859[i, j] = Taylor1(constant_term(-1.5) * constant_term(sin_ϕ[i, j]), order)
                    tmp861[i, j] = Taylor1(constant_term(2.5) * constant_term(sin3_ϕ[i, j]), order)
                    P_3_sin_ϕ[i, j] = Taylor1(constant_term(tmp859[i, j]) + constant_term(tmp861[i, j]), order)
                    tmp865[i, j] = Taylor1(constant_term(7.5) * constant_term(sin2_ϕ[i, j]), order)
                    ∂P_3_sin_ϕ[i, j] = Taylor1(constant_term(-1.5) + constant_term(tmp865[i, j]), order)
                    tmp867[j] = Taylor1(-(constant_term(Λ2[j])), order)
                    tmp869[i, j] = Taylor1(constant_term(r_p2[i, j]) ^ constant_term(2), order)
                    Λ2j_div_r4[i, j] = Taylor1(constant_term(tmp867[j]) / constant_term(tmp869[i, j]), order)
                    tmp871[j] = Taylor1(-(constant_term(Λ3[j])), order)
                    tmp873[i, j] = Taylor1(constant_term(r_p1d2[i, j]) ^ constant_term(5), order)
                    Λ3j_div_r5[i, j] = Taylor1(constant_term(tmp871[j]) / constant_term(tmp873[i, j]), order)
                    tmp875[i, j] = Taylor1(-(constant_term(cos_ϕ[i, j])), order)
                    m_c_ϕ_∂P_2[i, j] = Taylor1(constant_term(tmp875[i, j]) * constant_term(∂P_2_sin_ϕ[i, j]), order)
                    tmp877[i, j] = Taylor1(-(constant_term(cos_ϕ[i, j])), order)
                    m_c_ϕ_∂P_3[i, j] = Taylor1(constant_term(tmp877[i, j]) * constant_term(∂P_3_sin_ϕ[i, j]), order)
                    tmp880[i, j] = Taylor1(constant_term(Λ2j_div_r4[i, j]) * constant_term(3), order)
                    F_J2_ξ[i, j] = Taylor1(constant_term(tmp880[i, j]) * constant_term(P_2_sin_ϕ[i, j]), order)
                    F_J2_ζ[i, j] = Taylor1(constant_term(Λ2j_div_r4[i, j]) * constant_term(m_c_ϕ_∂P_2[i, j]), order)
                    tmp884[i, j] = Taylor1(constant_term(Λ3j_div_r5[i, j]) * constant_term(4), order)
                    F_J3_ξ[i, j] = Taylor1(constant_term(tmp884[i, j]) * constant_term(P_3_sin_ϕ[i, j]), order)
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
                    tmp919[i, j] = Taylor1(constant_term(μ[i]) * constant_term(F_J2_x[i, j]), order)
                    temp_accX_j[i, j] = Taylor1(constant_term(accX[j]) + constant_term(tmp919[i, j]), order)
                    accX[j] = Taylor1(identity(constant_term(temp_accX_j[i, j])), order)
                    tmp921[i, j] = Taylor1(constant_term(μ[i]) * constant_term(F_J2_y[i, j]), order)
                    temp_accY_j[i, j] = Taylor1(constant_term(accY[j]) + constant_term(tmp921[i, j]), order)
                    accY[j] = Taylor1(identity(constant_term(temp_accY_j[i, j])), order)
                    tmp923[i, j] = Taylor1(constant_term(μ[i]) * constant_term(F_J2_z[i, j]), order)
                    temp_accZ_j[i, j] = Taylor1(constant_term(accZ[j]) + constant_term(tmp923[i, j]), order)
                    accZ[j] = Taylor1(identity(constant_term(temp_accZ_j[i, j])), order)
                    tmp925[i, j] = Taylor1(constant_term(μ[j]) * constant_term(F_J2_x[i, j]), order)
                    temp_accX_i[i, j] = Taylor1(constant_term(accX[i]) - constant_term(tmp925[i, j]), order)
                    accX[i] = Taylor1(identity(constant_term(temp_accX_i[i, j])), order)
                    tmp927[i, j] = Taylor1(constant_term(μ[j]) * constant_term(F_J2_y[i, j]), order)
                    temp_accY_i[i, j] = Taylor1(constant_term(accY[i]) - constant_term(tmp927[i, j]), order)
                    accY[i] = Taylor1(identity(constant_term(temp_accY_i[i, j])), order)
                    tmp929[i, j] = Taylor1(constant_term(μ[j]) * constant_term(F_J2_z[i, j]), order)
                    temp_accZ_i[i, j] = Taylor1(constant_term(accZ[i]) - constant_term(tmp929[i, j]), order)
                    accZ[i] = Taylor1(identity(constant_term(temp_accZ_i[i, j])), order)
                end
            end
        end
        tmp932[3j - 2] = Taylor1(constant_term(dq[3j - 2]) ^ constant_term(2), order)
        tmp934[3j - 1] = Taylor1(constant_term(dq[3j - 1]) ^ constant_term(2), order)
        tmp935[3j - 2] = Taylor1(constant_term(tmp932[3j - 2]) + constant_term(tmp934[3j - 1]), order)
        tmp937[3j] = Taylor1(constant_term(dq[3j]) ^ constant_term(2), order)
        v2[j] = Taylor1(constant_term(tmp935[3j - 2]) + constant_term(tmp937[3j]), order)
    end
    for j = _1_to_N
        postNewtonX[j] = Taylor1(identity(constant_term(newtonX[j])), order)
        postNewtonY[j] = Taylor1(identity(constant_term(newtonY[j])), order)
        postNewtonZ[j] = Taylor1(identity(constant_term(newtonZ[j])), order)
    end
    tmp940 = Array{Taylor1{_S}}(undef, size(newtonianNb_Potential))
    tmp940 .= Taylor1(zero(_S), order)
    temp_005a = Array{Taylor1{_S}}(undef, size(newtonianNb_Potential))
    temp_005a .= Taylor1(zero(_S), order)
    tmp943 = Array{Taylor1{_S}}(undef, size(v2))
    tmp943 .= Taylor1(zero(_S), order)
    tmp945 = Array{Taylor1{_S}}(undef, size(vi_dot_vj))
    tmp945 .= Taylor1(zero(_S), order)
    temp_005b = Array{Taylor1{_S}}(undef, size(tmp945))
    temp_005b .= Taylor1(zero(_S), order)
    temp_005c = Array{Taylor1{_S}}(undef, size(temp_005b))
    temp_005c .= Taylor1(zero(_S), order)
    temp_005 = Array{Taylor1{_S}}(undef, size(temp_005c))
    temp_005 .= Taylor1(zero(_S), order)
    temp_006a = Array{Taylor1{_S}}(undef, size(X))
    temp_006a .= Taylor1(zero(_S), order)
    temp_006b = Array{Taylor1{_S}}(undef, size(Y))
    temp_006b .= Taylor1(zero(_S), order)
    temp_006c = Array{Taylor1{_S}}(undef, size(Z))
    temp_006c .= Taylor1(zero(_S), order)
    tmp952 = Array{Taylor1{_S}}(undef, size(temp_006a))
    tmp952 .= Taylor1(zero(_S), order)
    temp_006d = Array{Taylor1{_S}}(undef, size(tmp952))
    temp_006d .= Taylor1(zero(_S), order)
    tmp955 = Array{Taylor1{_S}}(undef, size(temp_006d))
    tmp955 .= Taylor1(zero(_S), order)
    temp_006e = Array{Taylor1{_S}}(undef, size(tmp955))
    temp_006e .= Taylor1(zero(_S), order)
    tmp958 = Array{Taylor1{_S}}(undef, size(temp_006e))
    tmp958 .= Taylor1(zero(_S), order)
    temp_006 = Array{Taylor1{_S}}(undef, size(temp_005))
    temp_006 .= Taylor1(zero(_S), order)
    temp_007a = Array{Taylor1{_S}}(undef, size(X))
    temp_007a .= Taylor1(zero(_S), order)
    temp_007b = Array{Taylor1{_S}}(undef, size(Y))
    temp_007b .= Taylor1(zero(_S), order)
    temp_007c = Array{Taylor1{_S}}(undef, size(Z))
    temp_007c .= Taylor1(zero(_S), order)
    tmp963 = Array{Taylor1{_S}}(undef, size(temp_007a))
    tmp963 .= Taylor1(zero(_S), order)
    temp_007d = Array{Taylor1{_S}}(undef, size(tmp963))
    temp_007d .= Taylor1(zero(_S), order)
    tmp966 = Array{Taylor1{_S}}(undef, size(temp_007d))
    tmp966 .= Taylor1(zero(_S), order)
    temp_007 = Array{Taylor1{_S}}(undef, size(temp_006))
    temp_007 .= Taylor1(zero(_S), order)
    temp_008 = Array{Taylor1{_S}}(undef, size(temp_007))
    temp_008 .= Taylor1(zero(_S), order)
    temp_009 = Array{Taylor1{_S}}(undef, size(X))
    temp_009 .= Taylor1(zero(_S), order)
    temp_010 = Array{Taylor1{_S}}(undef, size(Y))
    temp_010 .= Taylor1(zero(_S), order)
    temp_011 = Array{Taylor1{_S}}(undef, size(Z))
    temp_011 .= Taylor1(zero(_S), order)
    temp_013a = Array{Taylor1{_S}}(undef, size(pn2))
    temp_013a .= Taylor1(zero(_S), order)
    temp_013b = Array{Taylor1{_S}}(undef, size(pn3))
    temp_013b .= Taylor1(zero(_S), order)
    tmp977 = Array{Taylor1{_S}}(undef, size(temp_013a))
    tmp977 .= Taylor1(zero(_S), order)
    tmp978 = Array{Taylor1{_S}}(undef, size(temp_009))
    tmp978 .= Taylor1(zero(_S), order)
    temp_013 = Array{Taylor1{_S}}(undef, size(tmp978))
    temp_013 .= Taylor1(zero(_S), order)
    temp_014a = Array{Taylor1{_S}}(undef, size(pn2))
    temp_014a .= Taylor1(zero(_S), order)
    temp_014b = Array{Taylor1{_S}}(undef, size(pn3))
    temp_014b .= Taylor1(zero(_S), order)
    tmp982 = Array{Taylor1{_S}}(undef, size(temp_014a))
    tmp982 .= Taylor1(zero(_S), order)
    tmp983 = Array{Taylor1{_S}}(undef, size(temp_010))
    tmp983 .= Taylor1(zero(_S), order)
    temp_014 = Array{Taylor1{_S}}(undef, size(tmp983))
    temp_014 .= Taylor1(zero(_S), order)
    temp_015a = Array{Taylor1{_S}}(undef, size(pn2))
    temp_015a .= Taylor1(zero(_S), order)
    temp_015b = Array{Taylor1{_S}}(undef, size(pn3))
    temp_015b .= Taylor1(zero(_S), order)
    tmp987 = Array{Taylor1{_S}}(undef, size(temp_015a))
    tmp987 .= Taylor1(zero(_S), order)
    tmp988 = Array{Taylor1{_S}}(undef, size(temp_011))
    tmp988 .= Taylor1(zero(_S), order)
    temp_015 = Array{Taylor1{_S}}(undef, size(tmp988))
    temp_015 .= Taylor1(zero(_S), order)
    for k = Base.OneTo(succ_approx_iter)
        for j = _1_to_N
            pntempX[j] = Taylor1(identity(constant_term(zero_q_1)), order)
            pntempY[j] = Taylor1(identity(constant_term(zero_q_1)), order)
            pntempZ[j] = Taylor1(identity(constant_term(zero_q_1)), order)
        end
        for j = _1_to_N
            for i = _1_to_N
                if i == j
                else
                    tmp940[j] = Taylor1(constant_term(4) * constant_term(newtonianNb_Potential[j]), order)
                    temp_005a[i] = Taylor1(constant_term(newtonianNb_Potential[i]) + constant_term(tmp940[j]), order)
                    tmp943[i] = Taylor1(constant_term(2) * constant_term(v2[i]), order)
                    tmp945[i, j] = Taylor1(constant_term(4) * constant_term(vi_dot_vj[i, j]), order)
                    temp_005b[i, j] = Taylor1(constant_term(tmp943[i]) - constant_term(tmp945[i, j]), order)
                    temp_005c[i, j] = Taylor1(constant_term(v2[j]) + constant_term(temp_005b[i, j]), order)
                    temp_005[i, j] = Taylor1(constant_term(temp_005c[i, j]) - constant_term(temp_005a[i]), order)
                    temp_006a[i, j] = Taylor1(constant_term(X[i, j]) * constant_term(dq[3i - 2]), order)
                    temp_006b[i, j] = Taylor1(constant_term(Y[i, j]) * constant_term(dq[3i - 1]), order)
                    temp_006c[i, j] = Taylor1(constant_term(Z[i, j]) * constant_term(dq[3i]), order)
                    tmp952[i, j] = Taylor1(constant_term(temp_006a[i, j]) + constant_term(temp_006b[i, j]), order)
                    temp_006d[i, j] = Taylor1(constant_term(tmp952[i, j]) + constant_term(temp_006c[i, j]), order)
                    tmp955[i, j] = Taylor1(constant_term(temp_006d[i, j]) ^ constant_term(2), order)
                    temp_006e[i, j] = Taylor1(constant_term(tmp955[i, j]) / constant_term(r_p2[i, j]), order)
                    tmp958[i, j] = Taylor1(constant_term(1.5) * constant_term(temp_006e[i, j]), order)
                    temp_006[i, j] = Taylor1(constant_term(temp_005[i, j]) - constant_term(tmp958[i, j]), order)
                    temp_007a[i, j] = Taylor1(constant_term(X[i, j]) * constant_term(postNewtonX[i]), order)
                    temp_007b[i, j] = Taylor1(constant_term(Y[i, j]) * constant_term(postNewtonY[i]), order)
                    temp_007c[i, j] = Taylor1(constant_term(Z[i, j]) * constant_term(postNewtonZ[i]), order)
                    tmp963[i, j] = Taylor1(constant_term(temp_007a[i, j]) + constant_term(temp_007b[i, j]), order)
                    temp_007d[i, j] = Taylor1(constant_term(tmp963[i, j]) + constant_term(temp_007c[i, j]), order)
                    tmp966[i, j] = Taylor1(constant_term(0.5) * constant_term(temp_007d[i, j]), order)
                    temp_007[i, j] = Taylor1(constant_term(temp_006[i, j]) + constant_term(tmp966[i, j]), order)
                    temp_008[i, j] = Taylor1(constant_term(c_p2) + constant_term(temp_007[i, j]), order)
                    pn1[i, j] = Taylor1(constant_term(newtonianCoeff[i, j]) * constant_term(temp_008[i, j]), order)
                    temp_009[i, j] = Taylor1(constant_term(X[i, j]) * constant_term(pn1[i, j]), order)
                    temp_010[i, j] = Taylor1(constant_term(Y[i, j]) * constant_term(pn1[i, j]), order)
                    temp_011[i, j] = Taylor1(constant_term(Z[i, j]) * constant_term(pn1[i, j]), order)
                    pn3[i, j] = Taylor1(constant_term(3.5) * constant_term(newtonian1b_Potential[i, j]), order)
                    temp_013a[i, j] = Taylor1(constant_term(pn2[i, j]) * constant_term(U[i, j]), order)
                    temp_013b[i, j] = Taylor1(constant_term(pn3[i, j]) * constant_term(postNewtonX[i]), order)
                    tmp977[i, j] = Taylor1(constant_term(temp_013a[i, j]) + constant_term(temp_013b[i, j]), order)
                    tmp978[i, j] = Taylor1(constant_term(temp_009[i, j]) + constant_term(tmp977[i, j]), order)
                    temp_013[i, j] = Taylor1(constant_term(pntempX[j]) + constant_term(tmp978[i, j]), order)
                    pntempX[j] = Taylor1(identity(constant_term(temp_013[i, j])), order)
                    temp_014a[i, j] = Taylor1(constant_term(pn2[i, j]) * constant_term(V[i, j]), order)
                    temp_014b[i, j] = Taylor1(constant_term(pn3[i, j]) * constant_term(postNewtonY[i]), order)
                    tmp982[i, j] = Taylor1(constant_term(temp_014a[i, j]) + constant_term(temp_014b[i, j]), order)
                    tmp983[i, j] = Taylor1(constant_term(temp_010[i, j]) + constant_term(tmp982[i, j]), order)
                    temp_014[i, j] = Taylor1(constant_term(pntempY[j]) + constant_term(tmp983[i, j]), order)
                    pntempY[j] = Taylor1(identity(constant_term(temp_014[i, j])), order)
                    temp_015a[i, j] = Taylor1(constant_term(pn2[i, j]) * constant_term(W[i, j]), order)
                    temp_015b[i, j] = Taylor1(constant_term(pn3[i, j]) * constant_term(postNewtonZ[i]), order)
                    tmp987[i, j] = Taylor1(constant_term(temp_015a[i, j]) + constant_term(temp_015b[i, j]), order)
                    tmp988[i, j] = Taylor1(constant_term(temp_011[i, j]) + constant_term(tmp987[i, j]), order)
                    temp_015[i, j] = Taylor1(constant_term(pntempZ[j]) + constant_term(tmp988[i, j]), order)
                    pntempZ[j] = Taylor1(identity(constant_term(temp_015[i, j])), order)
                end
            end
        end
        for j = _1_to_N
            postNewtonX[j] = Taylor1(constant_term(pntempX[j]) * constant_term(c_m2), order)
            postNewtonY[j] = Taylor1(constant_term(pntempY[j]) * constant_term(c_m2), order)
            postNewtonZ[j] = Taylor1(constant_term(pntempZ[j]) * constant_term(c_m2), order)
        end
    end
    for i = _1_to_N
        dq[3 * (N + i) - 2] = Taylor1(constant_term(postNewtonX[i]) + constant_term(accX[i]), order)
        dq[3 * (N + i) - 1] = Taylor1(constant_term(postNewtonY[i]) + constant_term(accY[i]), order)
        dq[3 * (N + i)] = Taylor1(constant_term(postNewtonZ[i]) + constant_term(accZ[i]), order)
    end
    for __idx = eachindex(q)
        (q[__idx]).coeffs[2] = (dq[__idx]).coeffs[1]
    end
    for ord = 1:order - 1
        ordnext = ord + 1
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
                    TaylorSeries.identity!(temp_accX_j[i, j], zero_q_1, ord)
                    TaylorSeries.identity!(temp_accY_j[i, j], zero_q_1, ord)
                    TaylorSeries.identity!(temp_accZ_j[i, j], zero_q_1, ord)
                    TaylorSeries.identity!(temp_accX_i[i, j], zero_q_1, ord)
                    TaylorSeries.identity!(temp_accY_i[i, j], zero_q_1, ord)
                    TaylorSeries.identity!(temp_accZ_i[i, j], zero_q_1, ord)
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
                    TaylorSeries.mul!(tmp793[3j - 2], 4, dq[3j - 2], ord)
                    TaylorSeries.mul!(tmp795[3i - 2], 3, dq[3i - 2], ord)
                    TaylorSeries.subst!(_4U_m_3X[i, j], tmp793[3j - 2], tmp795[3i - 2], ord)
                    TaylorSeries.mul!(tmp798[3j - 1], 4, dq[3j - 1], ord)
                    TaylorSeries.mul!(tmp800[3i - 1], 3, dq[3i - 1], ord)
                    TaylorSeries.subst!(_4V_m_3Y[i, j], tmp798[3j - 1], tmp800[3i - 1], ord)
                    TaylorSeries.mul!(tmp803[3j], 4, dq[3j], ord)
                    TaylorSeries.mul!(tmp805[3i], 3, dq[3i], ord)
                    TaylorSeries.subst!(_4W_m_3Z[i, j], tmp803[3j], tmp805[3i], ord)
                    TaylorSeries.mul!(pn2x[i, j], X[i, j], _4U_m_3X[i, j], ord)
                    TaylorSeries.mul!(pn2y[i, j], Y[i, j], _4V_m_3Y[i, j], ord)
                    TaylorSeries.mul!(pn2z[i, j], Z[i, j], _4W_m_3Z[i, j], ord)
                    TaylorSeries.mul!(UU[i, j], dq[3i - 2], dq[3j - 2], ord)
                    TaylorSeries.mul!(VV[i, j], dq[3i - 1], dq[3j - 1], ord)
                    TaylorSeries.mul!(WW[i, j], dq[3i], dq[3j], ord)
                    TaylorSeries.add!(tmp813[i, j], UU[i, j], VV[i, j], ord)
                    TaylorSeries.add!(vi_dot_vj[i, j], tmp813[i, j], WW[i, j], ord)
                    TaylorSeries.pow!(tmp816[i, j], X[i, j], 2, ord)
                    TaylorSeries.pow!(tmp818[i, j], Y[i, j], 2, ord)
                    TaylorSeries.add!(tmp819[i, j], tmp816[i, j], tmp818[i, j], ord)
                    TaylorSeries.pow!(tmp821[i, j], Z[i, j], 2, ord)
                    TaylorSeries.add!(r_p2[i, j], tmp819[i, j], tmp821[i, j], ord)
                    TaylorSeries.sqrt!(r_p1d2[i, j], r_p2[i, j], ord)
                    TaylorSeries.pow!(r_p3d2[i, j], r_p2[i, j], 1.5, ord)
                    TaylorSeries.pow!(r_p7d2[i, j], r_p2[i, j], 3.5, ord)
                    TaylorSeries.div!(newtonianCoeff[i, j], μ[i], r_p3d2[i, j], ord)
                    TaylorSeries.add!(tmp829[i, j], pn2x[i, j], pn2y[i, j], ord)
                    TaylorSeries.add!(tmp830[i, j], tmp829[i, j], pn2z[i, j], ord)
                    TaylorSeries.mul!(pn2[i, j], newtonianCoeff[i, j], tmp830[i, j], ord)
                    TaylorSeries.mul!(tmp832[i, j], X[i, j], newtonianCoeff[i, j], ord)
                    TaylorSeries.add!(temp_001[i, j], newtonX[j], tmp832[i, j], ord)
                    TaylorSeries.identity!(newtonX[j], temp_001[i, j], ord)
                    TaylorSeries.mul!(tmp834[i, j], Y[i, j], newtonianCoeff[i, j], ord)
                    TaylorSeries.add!(temp_002[i, j], newtonY[j], tmp834[i, j], ord)
                    TaylorSeries.identity!(newtonY[j], temp_002[i, j], ord)
                    TaylorSeries.mul!(tmp836[i, j], Z[i, j], newtonianCoeff[i, j], ord)
                    TaylorSeries.add!(temp_003[i, j], newtonZ[j], tmp836[i, j], ord)
                    TaylorSeries.identity!(newtonZ[j], temp_003[i, j], ord)
                    TaylorSeries.div!(newtonian1b_Potential[i, j], μ[i], r_p1d2[i, j], ord)
                    TaylorSeries.add!(temp_004[i, j], newtonianNb_Potential[j], newtonian1b_Potential[i, j], ord)
                    TaylorSeries.identity!(newtonianNb_Potential[j], temp_004[i, j], ord)
                    if UJ_interaction[i, j]
                        TaylorSeries.mul!(t31[i, j], X[i, j], M_[1, 3, j], ord)
                        TaylorSeries.mul!(t32[i, j], Y[i, j], M_[2, 3, j], ord)
                        TaylorSeries.mul!(t33[i, j], Z[i, j], M_[3, 3, j], ord)
                        TaylorSeries.add!(tmp843[i, j], t31[i, j], t32[i, j], ord)
                        TaylorSeries.add!(r_sin_ϕ[i, j], tmp843[i, j], t33[i, j], ord)
                        TaylorSeries.div!(sin_ϕ[i, j], r_sin_ϕ[i, j], r_p1d2[i, j], ord)
                        TaylorSeries.asin!(ϕ[i, j], sin_ϕ[i, j], tmp996[i, j], ord)
                        TaylorSeries.sincos!(tmp997[i, j], cos_ϕ[i, j], ϕ[i, j], ord)
                        TaylorSeries.pow!(sin2_ϕ[i, j], sin_ϕ[i, j], 2, ord)
                        TaylorSeries.pow!(sin3_ϕ[i, j], sin_ϕ[i, j], 3, ord)
                        TaylorSeries.mul!(tmp853[i, j], 1.5, sin2_ϕ[i, j], ord)
                        TaylorSeries.subst!(P_2_sin_ϕ[i, j], tmp853[i, j], 0.5, ord)
                        TaylorSeries.mul!(∂P_2_sin_ϕ[i, j], 3, sin_ϕ[i, j], ord)
                        TaylorSeries.mul!(tmp859[i, j], -1.5, sin_ϕ[i, j], ord)
                        TaylorSeries.mul!(tmp861[i, j], 2.5, sin3_ϕ[i, j], ord)
                        TaylorSeries.add!(P_3_sin_ϕ[i, j], tmp859[i, j], tmp861[i, j], ord)
                        TaylorSeries.mul!(tmp865[i, j], 7.5, sin2_ϕ[i, j], ord)
                        TaylorSeries.add!(∂P_3_sin_ϕ[i, j], -1.5, tmp865[i, j], ord)
                        TaylorSeries.subst!(tmp867[j], Λ2[j], ord)
                        TaylorSeries.pow!(tmp869[i, j], r_p2[i, j], 2, ord)
                        TaylorSeries.div!(Λ2j_div_r4[i, j], tmp867[j], tmp869[i, j], ord)
                        TaylorSeries.subst!(tmp871[j], Λ3[j], ord)
                        TaylorSeries.pow!(tmp873[i, j], r_p1d2[i, j], 5, ord)
                        TaylorSeries.div!(Λ3j_div_r5[i, j], tmp871[j], tmp873[i, j], ord)
                        TaylorSeries.subst!(tmp875[i, j], cos_ϕ[i, j], ord)
                        TaylorSeries.mul!(m_c_ϕ_∂P_2[i, j], tmp875[i, j], ∂P_2_sin_ϕ[i, j], ord)
                        TaylorSeries.subst!(tmp877[i, j], cos_ϕ[i, j], ord)
                        TaylorSeries.mul!(m_c_ϕ_∂P_3[i, j], tmp877[i, j], ∂P_3_sin_ϕ[i, j], ord)
                        TaylorSeries.mul!(tmp880[i, j], Λ2j_div_r4[i, j], 3, ord)
                        TaylorSeries.mul!(F_J2_ξ[i, j], tmp880[i, j], P_2_sin_ϕ[i, j], ord)
                        TaylorSeries.mul!(F_J2_ζ[i, j], Λ2j_div_r4[i, j], m_c_ϕ_∂P_2[i, j], ord)
                        TaylorSeries.mul!(tmp884[i, j], Λ3j_div_r5[i, j], 4, ord)
                        TaylorSeries.mul!(F_J3_ξ[i, j], tmp884[i, j], P_3_sin_ϕ[i, j], ord)
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
                        TaylorSeries.mul!(tmp919[i, j], μ[i], F_J2_x[i, j], ord)
                        TaylorSeries.add!(temp_accX_j[i, j], accX[j], tmp919[i, j], ord)
                        TaylorSeries.identity!(accX[j], temp_accX_j[i, j], ord)
                        TaylorSeries.mul!(tmp921[i, j], μ[i], F_J2_y[i, j], ord)
                        TaylorSeries.add!(temp_accY_j[i, j], accY[j], tmp921[i, j], ord)
                        TaylorSeries.identity!(accY[j], temp_accY_j[i, j], ord)
                        TaylorSeries.mul!(tmp923[i, j], μ[i], F_J2_z[i, j], ord)
                        TaylorSeries.add!(temp_accZ_j[i, j], accZ[j], tmp923[i, j], ord)
                        TaylorSeries.identity!(accZ[j], temp_accZ_j[i, j], ord)
                        TaylorSeries.mul!(tmp925[i, j], μ[j], F_J2_x[i, j], ord)
                        TaylorSeries.subst!(temp_accX_i[i, j], accX[i], tmp925[i, j], ord)
                        TaylorSeries.identity!(accX[i], temp_accX_i[i, j], ord)
                        TaylorSeries.mul!(tmp927[i, j], μ[j], F_J2_y[i, j], ord)
                        TaylorSeries.subst!(temp_accY_i[i, j], accY[i], tmp927[i, j], ord)
                        TaylorSeries.identity!(accY[i], temp_accY_i[i, j], ord)
                        TaylorSeries.mul!(tmp929[i, j], μ[j], F_J2_z[i, j], ord)
                        TaylorSeries.subst!(temp_accZ_i[i, j], accZ[i], tmp929[i, j], ord)
                        TaylorSeries.identity!(accZ[i], temp_accZ_i[i, j], ord)
                    end
                end
            end
            TaylorSeries.pow!(tmp932[3j - 2], dq[3j - 2], 2, ord)
            TaylorSeries.pow!(tmp934[3j - 1], dq[3j - 1], 2, ord)
            TaylorSeries.add!(tmp935[3j - 2], tmp932[3j - 2], tmp934[3j - 1], ord)
            TaylorSeries.pow!(tmp937[3j], dq[3j], 2, ord)
            TaylorSeries.add!(v2[j], tmp935[3j - 2], tmp937[3j], ord)
        end
        for j = _1_to_N
            TaylorSeries.identity!(postNewtonX[j], newtonX[j], ord)
            TaylorSeries.identity!(postNewtonY[j], newtonY[j], ord)
            TaylorSeries.identity!(postNewtonZ[j], newtonZ[j], ord)
        end
        for k = Base.OneTo(succ_approx_iter)
            for j = _1_to_N
                TaylorSeries.identity!(pntempX[j], zero_q_1, ord)
                TaylorSeries.identity!(pntempY[j], zero_q_1, ord)
                TaylorSeries.identity!(pntempZ[j], zero_q_1, ord)
            end
            for j = _1_to_N
                for i = _1_to_N
                    if i == j
                    else
                        TaylorSeries.mul!(tmp940[j], 4, newtonianNb_Potential[j], ord)
                        TaylorSeries.add!(temp_005a[i], newtonianNb_Potential[i], tmp940[j], ord)
                        TaylorSeries.mul!(tmp943[i], 2, v2[i], ord)
                        TaylorSeries.mul!(tmp945[i, j], 4, vi_dot_vj[i, j], ord)
                        TaylorSeries.subst!(temp_005b[i, j], tmp943[i], tmp945[i, j], ord)
                        TaylorSeries.add!(temp_005c[i, j], v2[j], temp_005b[i, j], ord)
                        TaylorSeries.subst!(temp_005[i, j], temp_005c[i, j], temp_005a[i], ord)
                        TaylorSeries.mul!(temp_006a[i, j], X[i, j], dq[3i - 2], ord)
                        TaylorSeries.mul!(temp_006b[i, j], Y[i, j], dq[3i - 1], ord)
                        TaylorSeries.mul!(temp_006c[i, j], Z[i, j], dq[3i], ord)
                        TaylorSeries.add!(tmp952[i, j], temp_006a[i, j], temp_006b[i, j], ord)
                        TaylorSeries.add!(temp_006d[i, j], tmp952[i, j], temp_006c[i, j], ord)
                        TaylorSeries.pow!(tmp955[i, j], temp_006d[i, j], 2, ord)
                        TaylorSeries.div!(temp_006e[i, j], tmp955[i, j], r_p2[i, j], ord)
                        TaylorSeries.mul!(tmp958[i, j], 1.5, temp_006e[i, j], ord)
                        TaylorSeries.subst!(temp_006[i, j], temp_005[i, j], tmp958[i, j], ord)
                        TaylorSeries.mul!(temp_007a[i, j], X[i, j], postNewtonX[i], ord)
                        TaylorSeries.mul!(temp_007b[i, j], Y[i, j], postNewtonY[i], ord)
                        TaylorSeries.mul!(temp_007c[i, j], Z[i, j], postNewtonZ[i], ord)
                        TaylorSeries.add!(tmp963[i, j], temp_007a[i, j], temp_007b[i, j], ord)
                        TaylorSeries.add!(temp_007d[i, j], tmp963[i, j], temp_007c[i, j], ord)
                        TaylorSeries.mul!(tmp966[i, j], 0.5, temp_007d[i, j], ord)
                        TaylorSeries.add!(temp_007[i, j], temp_006[i, j], tmp966[i, j], ord)
                        TaylorSeries.add!(temp_008[i, j], c_p2, temp_007[i, j], ord)
                        TaylorSeries.mul!(pn1[i, j], newtonianCoeff[i, j], temp_008[i, j], ord)
                        TaylorSeries.mul!(temp_009[i, j], X[i, j], pn1[i, j], ord)
                        TaylorSeries.mul!(temp_010[i, j], Y[i, j], pn1[i, j], ord)
                        TaylorSeries.mul!(temp_011[i, j], Z[i, j], pn1[i, j], ord)
                        TaylorSeries.mul!(pn3[i, j], 3.5, newtonian1b_Potential[i, j], ord)
                        TaylorSeries.mul!(temp_013a[i, j], pn2[i, j], U[i, j], ord)
                        TaylorSeries.mul!(temp_013b[i, j], pn3[i, j], postNewtonX[i], ord)
                        TaylorSeries.add!(tmp977[i, j], temp_013a[i, j], temp_013b[i, j], ord)
                        TaylorSeries.add!(tmp978[i, j], temp_009[i, j], tmp977[i, j], ord)
                        TaylorSeries.add!(temp_013[i, j], pntempX[j], tmp978[i, j], ord)
                        TaylorSeries.identity!(pntempX[j], temp_013[i, j], ord)
                        TaylorSeries.mul!(temp_014a[i, j], pn2[i, j], V[i, j], ord)
                        TaylorSeries.mul!(temp_014b[i, j], pn3[i, j], postNewtonY[i], ord)
                        TaylorSeries.add!(tmp982[i, j], temp_014a[i, j], temp_014b[i, j], ord)
                        TaylorSeries.add!(tmp983[i, j], temp_010[i, j], tmp982[i, j], ord)
                        TaylorSeries.add!(temp_014[i, j], pntempY[j], tmp983[i, j], ord)
                        TaylorSeries.identity!(pntempY[j], temp_014[i, j], ord)
                        TaylorSeries.mul!(temp_015a[i, j], pn2[i, j], W[i, j], ord)
                        TaylorSeries.mul!(temp_015b[i, j], pn3[i, j], postNewtonZ[i], ord)
                        TaylorSeries.add!(tmp987[i, j], temp_015a[i, j], temp_015b[i, j], ord)
                        TaylorSeries.add!(tmp988[i, j], temp_011[i, j], tmp987[i, j], ord)
                        TaylorSeries.add!(temp_015[i, j], pntempZ[j], tmp988[i, j], ord)
                        TaylorSeries.identity!(pntempZ[j], temp_015[i, j], ord)
                    end
                end
            end
            for j = _1_to_N
                TaylorSeries.mul!(postNewtonX[j], pntempX[j], c_m2, ord)
                TaylorSeries.mul!(postNewtonY[j], pntempY[j], c_m2, ord)
                TaylorSeries.mul!(postNewtonZ[j], pntempZ[j], c_m2, ord)
            end
        end
        for i = _1_to_N
            TaylorSeries.add!(dq[3 * (N + i) - 2], postNewtonX[i], accX[i], ord)
            TaylorSeries.add!(dq[3 * (N + i) - 1], postNewtonY[i], accY[i], ord)
            TaylorSeries.add!(dq[3 * (N + i)], postNewtonZ[i], accZ[i], ord)
        end
        for __idx = eachindex(q)
            (q[__idx]).coeffs[ordnext + 1] = (dq[__idx]).coeffs[ordnext] / ordnext
        end
    end
    return nothing
end
