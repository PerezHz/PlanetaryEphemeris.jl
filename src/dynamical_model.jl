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

    #TODO: handle appropiately "taylorized" version with postnewton_iter>1
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

    postNewtonX = Array{Taylor1{S}}(undef, N, postnewton_iter+1)
    postNewtonY = Array{Taylor1{S}}(undef, N, postnewton_iter+1)
    postNewtonZ = Array{Taylor1{S}}(undef, N, postnewton_iter+1)

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
    _4ϕj = Array{Taylor1{S}}(undef, N, N)
    ϕi_plus_4ϕj = Array{Taylor1{S}}(undef, N, N)
    sj2_plus_2si2_minus_4vivj = Array{Taylor1{S}}(undef, N, N)
    ϕs_and_vs = Array{Taylor1{S}}(undef, N, N)
    X_t_pn1 = Array{Taylor1{S}}(undef, N, N)
    Y_t_pn1 = Array{Taylor1{S}}(undef, N, N)
    Z_t_pn1 = Array{Taylor1{S}}(undef, N, N)
    U_t_pn2 = Array{Taylor1{S}}(undef, N, N)
    V_t_pn2 = Array{Taylor1{S}}(undef, N, N)
    W_t_pn2 = Array{Taylor1{S}}(undef, N, N)
    pNX_t_pn3 = Array{Taylor1{S}}(undef, N, N)
    pNY_t_pn3 = Array{Taylor1{S}}(undef, N, N)
    pNZ_t_pn3 = Array{Taylor1{S}}(undef, N, N)
    pNX_t_X = Array{Taylor1{S}}(undef, N, N)
    pNY_t_Y = Array{Taylor1{S}}(undef, N, N)
    pNZ_t_Z = Array{Taylor1{S}}(undef, N, N)

    for j in _1_to_N
        for i in _1_to_N
            # i == j && continue
            if i == j
            else
                pn1[i,j] = zero_q_1
                pn2[i,j] = zero_q_1
                pn3[i,j] = zero_q_1
                _4ϕj[i,j] = zero_q_1
                ϕi_plus_4ϕj[i,j] = zero_q_1
                sj2_plus_2si2_minus_4vivj[i,j] = zero_q_1
                ϕs_and_vs[i,j] = zero_q_1
                X_t_pn1[i,j] = zero_q_1
                Y_t_pn1[i,j] = zero_q_1
                Z_t_pn1[i,j] = zero_q_1
                UU[i,j] = zero_q_1
                VV[i,j] = zero_q_1
                WW[i,j] = zero_q_1
                U_t_pn2[i,j] = zero_q_1
                V_t_pn2[i,j] = zero_q_1
                W_t_pn2[i,j] = zero_q_1
                pNX_t_pn3[i,j] = zero_q_1
                pNY_t_pn3[i,j] = zero_q_1
                pNZ_t_pn3[i,j] = zero_q_1
                pNX_t_X[i,j] = zero_q_1
                pNY_t_Y[i,j] = zero_q_1
                pNZ_t_Z[i,j] = zero_q_1
            end
        end
        vi_dot_vj[j] = zero_q_1
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

                newtonX[j] = newtonX[j] + (X[i,j]*newtonianCoeff[i,j])
                newtonY[j] = newtonY[j] + (Y[i,j]*newtonianCoeff[i,j])
                newtonZ[j] = newtonZ[j] + (Z[i,j]*newtonianCoeff[i,j])

                newtonian1b_Potential[i,j] = μ[i]/r_p1d2[i, j]

                newtonianNb_Potential[j] = newtonianNb_Potential[j] + newtonian1b_Potential[i, j]

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
                    accX[j] = accX[j] + (μ[i]*F_J2_x[i,j])
                    accY[j] = accY[j] + (μ[i]*F_J2_y[i,j])
                    accZ[j] = accZ[j] + (μ[i]*F_J2_z[i,j])

                    # # reaction force on i-th body
                    # @show "acc",i,"-μ",j,"Λ2",j
                    accX[i] = accX[i] - (μ[j]*F_J2_x[i,j])
                    accY[i] = accY[i] - (μ[j]*F_J2_y[i,j])
                    accZ[i] = accZ[i] - (μ[j]*F_J2_z[i,j])
                end
            end #if i != j
        end #for, i
        v2[j] = ( (dq[3j-2]^2)+(dq[3j-1]^2) ) + (dq[3j]^2)
    end #for, j

    for j in _1_to_N
        postNewtonX[j,1] = newtonX[j]
        postNewtonY[j,1] = newtonY[j]
        postNewtonZ[j,1] = newtonZ[j]
    end

    for k in Base.OneTo(postnewton_iter)
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
                    pNX_t_X[i,j] = postNewtonX[i,k]*X[i,j]
                    pNY_t_Y[i,j] = postNewtonY[i,k]*Y[i,j]
                    pNZ_t_Z[i,j] = postNewtonZ[i,k]*Z[i,j]
                    pn1[i,j] = newtonianCoeff[i,j]*(  c_p2  +  ( pn1t2_7 + ((pNX_t_X[i,j]+pNY_t_Y[i,j])+pNZ_t_Z[i,j]) )  )

                    X_t_pn1[i,j] = X[i,j]*pn1[i,j]
                    Y_t_pn1[i,j] = Y[i,j]*pn1[i,j]
                    Z_t_pn1[i,j] = Z[i,j]*pn1[i,j]

                    pn3[i,j] = 3.5newtonian1b_Potential[i,j]

                    U_t_pn2[i,j] = pn2[i,j]*U[i,j]
                    pNX_t_pn3[i,j] = postNewtonX[i,k]*pn3[i,j]
                    pntempX[j] = pntempX[j] + ( X_t_pn1[i,j] + (U_t_pn2[i,j]+pNX_t_pn3[i,j]) )

                    V_t_pn2[i,j] = pn2[i,j]*V[i,j]
                    pNY_t_pn3[i,j] = postNewtonY[i,k]*pn3[i,j]
                    pntempY[j] = pntempY[j] + ( Y_t_pn1[i,j] + (V_t_pn2[i,j]+pNY_t_pn3[i,j]) )

                    W_t_pn2[i,j] = pn2[i,j]*W[i,j]
                    pNZ_t_pn3[i,j] = postNewtonZ[i,k]*pn3[i,j]
                    pntempZ[j] = pntempZ[j] + ( Z_t_pn1[i,j] + (W_t_pn2[i,j]+pNZ_t_pn3[i,j]) )
                end # if i != j
            end #for i
        end #for j
        for j in _1_to_N
            postNewtonX[j,k+1] = pntempX[j]*c_m2
            postNewtonY[j,k+1] = pntempY[j]*c_m2
            postNewtonZ[j,k+1] = pntempZ[j]*c_m2
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
    postNewtonX = Array{Taylor1{S}}(undef, N, postnewton_iter + 1)
    postNewtonY = Array{Taylor1{S}}(undef, N, postnewton_iter + 1)
    postNewtonZ = Array{Taylor1{S}}(undef, N, postnewton_iter + 1)
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
    _4ϕj = Array{Taylor1{S}}(undef, N, N)
    ϕi_plus_4ϕj = Array{Taylor1{S}}(undef, N, N)
    sj2_plus_2si2_minus_4vivj = Array{Taylor1{S}}(undef, N, N)
    ϕs_and_vs = Array{Taylor1{S}}(undef, N, N)
    X_t_pn1 = Array{Taylor1{S}}(undef, N, N)
    Y_t_pn1 = Array{Taylor1{S}}(undef, N, N)
    Z_t_pn1 = Array{Taylor1{S}}(undef, N, N)
    U_t_pn2 = Array{Taylor1{S}}(undef, N, N)
    V_t_pn2 = Array{Taylor1{S}}(undef, N, N)
    W_t_pn2 = Array{Taylor1{S}}(undef, N, N)
    pNX_t_pn3 = Array{Taylor1{S}}(undef, N, N)
    pNY_t_pn3 = Array{Taylor1{S}}(undef, N, N)
    pNZ_t_pn3 = Array{Taylor1{S}}(undef, N, N)
    pNX_t_X = Array{Taylor1{S}}(undef, N, N)
    pNY_t_Y = Array{Taylor1{S}}(undef, N, N)
    pNZ_t_Z = Array{Taylor1{S}}(undef, N, N)
    for j = _1_to_N
        for i = _1_to_N
            if i == j
            else
                pn1[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                pn2[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                pn3[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                _4ϕj[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                ϕi_plus_4ϕj[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                sj2_plus_2si2_minus_4vivj[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                ϕs_and_vs[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                X_t_pn1[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                Y_t_pn1[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                Z_t_pn1[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                UU[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                VV[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                WW[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                U_t_pn2[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                V_t_pn2[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                W_t_pn2[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                pNX_t_pn3[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                pNY_t_pn3[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                pNZ_t_pn3[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                pNX_t_X[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                pNY_t_Y[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
                pNZ_t_Z[i, j] = Taylor1(identity(constant_term(zero_q_1)), order)
            end
        end
        vi_dot_vj[j] = Taylor1(identity(constant_term(zero_q_1)), order)
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
    tmp2184 = Array{Taylor1{_S}}(undef, size(dq))
    tmp2184 .= Taylor1(zero(_S), order)
    tmp2186 = Array{Taylor1{_S}}(undef, size(dq))
    tmp2186 .= Taylor1(zero(_S), order)
    tmp2189 = Array{Taylor1{_S}}(undef, size(dq))
    tmp2189 .= Taylor1(zero(_S), order)
    tmp2191 = Array{Taylor1{_S}}(undef, size(dq))
    tmp2191 .= Taylor1(zero(_S), order)
    tmp2194 = Array{Taylor1{_S}}(undef, size(dq))
    tmp2194 .= Taylor1(zero(_S), order)
    tmp2196 = Array{Taylor1{_S}}(undef, size(dq))
    tmp2196 .= Taylor1(zero(_S), order)
    pn2x = Array{Taylor1{_S}}(undef, size(X))
    pn2x .= Taylor1(zero(_S), order)
    pn2y = Array{Taylor1{_S}}(undef, size(Y))
    pn2y .= Taylor1(zero(_S), order)
    pn2z = Array{Taylor1{_S}}(undef, size(Z))
    pn2z .= Taylor1(zero(_S), order)
    tmp2204 = Array{Taylor1{_S}}(undef, size(UU))
    tmp2204 .= Taylor1(zero(_S), order)
    tmp2207 = Array{Taylor1{_S}}(undef, size(X))
    tmp2207 .= Taylor1(zero(_S), order)
    tmp2209 = Array{Taylor1{_S}}(undef, size(Y))
    tmp2209 .= Taylor1(zero(_S), order)
    tmp2210 = Array{Taylor1{_S}}(undef, size(tmp2207))
    tmp2210 .= Taylor1(zero(_S), order)
    tmp2212 = Array{Taylor1{_S}}(undef, size(Z))
    tmp2212 .= Taylor1(zero(_S), order)
    tmp2220 = Array{Taylor1{_S}}(undef, size(pn2x))
    tmp2220 .= Taylor1(zero(_S), order)
    tmp2221 = Array{Taylor1{_S}}(undef, size(tmp2220))
    tmp2221 .= Taylor1(zero(_S), order)
    tmp2234 = Array{Taylor1{_S}}(undef, size(t31))
    tmp2234 .= Taylor1(zero(_S), order)
    tmp2385 = Array{Taylor1{_S}}(undef, size(sin_ϕ))
    tmp2385 .= Taylor1(zero(_S), order)
    tmp2386 = Array{Taylor1{_S}}(undef, size(ϕ))
    tmp2386 .= Taylor1(zero(_S), order)
    tmp2244 = Array{Taylor1{_S}}(undef, size(sin2_ϕ))
    tmp2244 .= Taylor1(zero(_S), order)
    tmp2250 = Array{Taylor1{_S}}(undef, size(sin_ϕ))
    tmp2250 .= Taylor1(zero(_S), order)
    tmp2252 = Array{Taylor1{_S}}(undef, size(sin3_ϕ))
    tmp2252 .= Taylor1(zero(_S), order)
    tmp2256 = Array{Taylor1{_S}}(undef, size(sin2_ϕ))
    tmp2256 .= Taylor1(zero(_S), order)
    tmp2258 = Array{Taylor1{_S}}(undef, size(Λ2))
    tmp2258 .= Taylor1(zero(_S), order)
    tmp2260 = Array{Taylor1{_S}}(undef, size(r_p2))
    tmp2260 .= Taylor1(zero(_S), order)
    tmp2262 = Array{Taylor1{_S}}(undef, size(Λ3))
    tmp2262 .= Taylor1(zero(_S), order)
    tmp2264 = Array{Taylor1{_S}}(undef, size(r_p1d2))
    tmp2264 .= Taylor1(zero(_S), order)
    tmp2266 = Array{Taylor1{_S}}(undef, size(cos_ϕ))
    tmp2266 .= Taylor1(zero(_S), order)
    tmp2268 = Array{Taylor1{_S}}(undef, size(cos_ϕ))
    tmp2268 .= Taylor1(zero(_S), order)
    tmp2271 = Array{Taylor1{_S}}(undef, size(Λ2j_div_r4))
    tmp2271 .= Taylor1(zero(_S), order)
    tmp2275 = Array{Taylor1{_S}}(undef, size(Λ3j_div_r5))
    tmp2275 .= Taylor1(zero(_S), order)
    tmp2323 = Array{Taylor1{_S}}(undef, size(dq))
    tmp2323 .= Taylor1(zero(_S), order)
    tmp2325 = Array{Taylor1{_S}}(undef, size(dq))
    tmp2325 .= Taylor1(zero(_S), order)
    tmp2326 = Array{Taylor1{_S}}(undef, size(tmp2323))
    tmp2326 .= Taylor1(zero(_S), order)
    tmp2328 = Array{Taylor1{_S}}(undef, size(dq))
    tmp2328 .= Taylor1(zero(_S), order)
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
                tmp2184[3j - 2] = Taylor1(constant_term(4) * constant_term(dq[3j - 2]), order)
                tmp2186[3i - 2] = Taylor1(constant_term(3) * constant_term(dq[3i - 2]), order)
                _4U_m_3X[i, j] = Taylor1(constant_term(tmp2184[3j - 2]) - constant_term(tmp2186[3i - 2]), order)
                tmp2189[3j - 1] = Taylor1(constant_term(4) * constant_term(dq[3j - 1]), order)
                tmp2191[3i - 1] = Taylor1(constant_term(3) * constant_term(dq[3i - 1]), order)
                _4V_m_3Y[i, j] = Taylor1(constant_term(tmp2189[3j - 1]) - constant_term(tmp2191[3i - 1]), order)
                tmp2194[3j] = Taylor1(constant_term(4) * constant_term(dq[3j]), order)
                tmp2196[3i] = Taylor1(constant_term(3) * constant_term(dq[3i]), order)
                _4W_m_3Z[i, j] = Taylor1(constant_term(tmp2194[3j]) - constant_term(tmp2196[3i]), order)
                pn2x[i, j] = Taylor1(constant_term(X[i, j]) * constant_term(_4U_m_3X[i, j]), order)
                pn2y[i, j] = Taylor1(constant_term(Y[i, j]) * constant_term(_4V_m_3Y[i, j]), order)
                pn2z[i, j] = Taylor1(constant_term(Z[i, j]) * constant_term(_4W_m_3Z[i, j]), order)
                UU[i, j] = Taylor1(constant_term(dq[3i - 2]) * constant_term(dq[3j - 2]), order)
                VV[i, j] = Taylor1(constant_term(dq[3i - 1]) * constant_term(dq[3j - 1]), order)
                WW[i, j] = Taylor1(constant_term(dq[3i]) * constant_term(dq[3j]), order)
                tmp2204[i, j] = Taylor1(constant_term(UU[i, j]) + constant_term(VV[i, j]), order)
                vi_dot_vj[i, j] = Taylor1(constant_term(tmp2204[i, j]) + constant_term(WW[i, j]), order)
                tmp2207[i, j] = Taylor1(constant_term(X[i, j]) ^ constant_term(2), order)
                tmp2209[i, j] = Taylor1(constant_term(Y[i, j]) ^ constant_term(2), order)
                tmp2210[i, j] = Taylor1(constant_term(tmp2207[i, j]) + constant_term(tmp2209[i, j]), order)
                tmp2212[i, j] = Taylor1(constant_term(Z[i, j]) ^ constant_term(2), order)
                r_p2[i, j] = Taylor1(constant_term(tmp2210[i, j]) + constant_term(tmp2212[i, j]), order)
                r_p1d2[i, j] = Taylor1(sqrt(constant_term(r_p2[i, j])), order)
                r_p3d2[i, j] = Taylor1(constant_term(r_p2[i, j]) ^ constant_term(1.5), order)
                r_p7d2[i, j] = Taylor1(constant_term(r_p2[i, j]) ^ constant_term(3.5), order)
                newtonianCoeff[i, j] = Taylor1(constant_term(μ[i]) / constant_term(r_p3d2[i, j]), order)
                tmp2220[i, j] = Taylor1(constant_term(pn2x[i, j]) + constant_term(pn2y[i, j]), order)
                tmp2221[i, j] = Taylor1(constant_term(tmp2220[i, j]) + constant_term(pn2z[i, j]), order)
                pn2[i, j] = Taylor1(constant_term(newtonianCoeff[i, j]) * constant_term(tmp2221[i, j]), order)
                newtonX[j] = Taylor1(constant_term(newtonX[j]) + constant_term(X[i, j] * newtonianCoeff[i, j]), order)
                newtonY[j] = Taylor1(constant_term(newtonY[j]) + constant_term(Y[i, j] * newtonianCoeff[i, j]), order)
                newtonZ[j] = Taylor1(constant_term(newtonZ[j]) + constant_term(Z[i, j] * newtonianCoeff[i, j]), order)
                newtonian1b_Potential[i, j] = Taylor1(constant_term(μ[i]) / constant_term(r_p1d2[i, j]), order)
                newtonianNb_Potential[j] = Taylor1(constant_term(newtonianNb_Potential[j]) + constant_term(newtonian1b_Potential[i, j]), order)
                if UJ_interaction[i, j]
                    t31[i, j] = Taylor1(constant_term(X[i, j]) * constant_term(M_[1, 3, j]), order)
                    t32[i, j] = Taylor1(constant_term(Y[i, j]) * constant_term(M_[2, 3, j]), order)
                    t33[i, j] = Taylor1(constant_term(Z[i, j]) * constant_term(M_[3, 3, j]), order)
                    tmp2234[i, j] = Taylor1(constant_term(t31[i, j]) + constant_term(t32[i, j]), order)
                    r_sin_ϕ[i, j] = Taylor1(constant_term(tmp2234[i, j]) + constant_term(t33[i, j]), order)
                    sin_ϕ[i, j] = Taylor1(constant_term(r_sin_ϕ[i, j]) / constant_term(r_p1d2[i, j]), order)
                    ϕ[i, j] = Taylor1(asin(constant_term(sin_ϕ[i, j])), order)
                    tmp2385[i, j] = Taylor1(sqrt(1 - constant_term(sin_ϕ[i, j]) ^ 2), order)
                    cos_ϕ[i, j] = Taylor1(cos(constant_term(ϕ[i, j])), order)
                    tmp2386[i, j] = Taylor1(sin(constant_term(ϕ[i, j])), order)
                    sin2_ϕ[i, j] = Taylor1(constant_term(sin_ϕ[i, j]) ^ constant_term(2), order)
                    sin3_ϕ[i, j] = Taylor1(constant_term(sin_ϕ[i, j]) ^ constant_term(3), order)
                    tmp2244[i, j] = Taylor1(constant_term(1.5) * constant_term(sin2_ϕ[i, j]), order)
                    P_2_sin_ϕ[i, j] = Taylor1(constant_term(tmp2244[i, j]) - constant_term(0.5), order)
                    ∂P_2_sin_ϕ[i, j] = Taylor1(constant_term(3) * constant_term(sin_ϕ[i, j]), order)
                    tmp2250[i, j] = Taylor1(constant_term(-1.5) * constant_term(sin_ϕ[i, j]), order)
                    tmp2252[i, j] = Taylor1(constant_term(2.5) * constant_term(sin3_ϕ[i, j]), order)
                    P_3_sin_ϕ[i, j] = Taylor1(constant_term(tmp2250[i, j]) + constant_term(tmp2252[i, j]), order)
                    tmp2256[i, j] = Taylor1(constant_term(7.5) * constant_term(sin2_ϕ[i, j]), order)
                    ∂P_3_sin_ϕ[i, j] = Taylor1(constant_term(-1.5) + constant_term(tmp2256[i, j]), order)
                    tmp2258[j] = Taylor1(-(constant_term(Λ2[j])), order)
                    tmp2260[i, j] = Taylor1(constant_term(r_p2[i, j]) ^ constant_term(2), order)
                    Λ2j_div_r4[i, j] = Taylor1(constant_term(tmp2258[j]) / constant_term(tmp2260[i, j]), order)
                    tmp2262[j] = Taylor1(-(constant_term(Λ3[j])), order)
                    tmp2264[i, j] = Taylor1(constant_term(r_p1d2[i, j]) ^ constant_term(5), order)
                    Λ3j_div_r5[i, j] = Taylor1(constant_term(tmp2262[j]) / constant_term(tmp2264[i, j]), order)
                    tmp2266[i, j] = Taylor1(-(constant_term(cos_ϕ[i, j])), order)
                    m_c_ϕ_∂P_2[i, j] = Taylor1(constant_term(tmp2266[i, j]) * constant_term(∂P_2_sin_ϕ[i, j]), order)
                    tmp2268[i, j] = Taylor1(-(constant_term(cos_ϕ[i, j])), order)
                    m_c_ϕ_∂P_3[i, j] = Taylor1(constant_term(tmp2268[i, j]) * constant_term(∂P_3_sin_ϕ[i, j]), order)
                    tmp2271[i, j] = Taylor1(constant_term(Λ2j_div_r4[i, j]) * constant_term(3), order)
                    F_J2_ξ[i, j] = Taylor1(constant_term(tmp2271[i, j]) * constant_term(P_2_sin_ϕ[i, j]), order)
                    F_J2_ζ[i, j] = Taylor1(constant_term(Λ2j_div_r4[i, j]) * constant_term(m_c_ϕ_∂P_2[i, j]), order)
                    tmp2275[i, j] = Taylor1(constant_term(Λ3j_div_r5[i, j]) * constant_term(4), order)
                    F_J3_ξ[i, j] = Taylor1(constant_term(tmp2275[i, j]) * constant_term(P_3_sin_ϕ[i, j]), order)
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
                    accX[j] = Taylor1(constant_term(accX[j]) + constant_term(μ[i] * F_J2_x[i, j]), order)
                    accY[j] = Taylor1(constant_term(accY[j]) + constant_term(μ[i] * F_J2_y[i, j]), order)
                    accZ[j] = Taylor1(constant_term(accZ[j]) + constant_term(μ[i] * F_J2_z[i, j]), order)
                    accX[i] = Taylor1(constant_term(accX[i]) - constant_term(μ[j] * F_J2_x[i, j]), order)
                    accY[i] = Taylor1(constant_term(accY[i]) - constant_term(μ[j] * F_J2_y[i, j]), order)
                    accZ[i] = Taylor1(constant_term(accZ[i]) - constant_term(μ[j] * F_J2_z[i, j]), order)
                end
            end
        end
        tmp2323[3j - 2] = Taylor1(constant_term(dq[3j - 2]) ^ constant_term(2), order)
        tmp2325[3j - 1] = Taylor1(constant_term(dq[3j - 1]) ^ constant_term(2), order)
        tmp2326[3j - 2] = Taylor1(constant_term(tmp2323[3j - 2]) + constant_term(tmp2325[3j - 1]), order)
        tmp2328[3j] = Taylor1(constant_term(dq[3j]) ^ constant_term(2), order)
        v2[j] = Taylor1(constant_term(tmp2326[3j - 2]) + constant_term(tmp2328[3j]), order)
    end
    for j = _1_to_N
        postNewtonX[j, 1] = Taylor1(identity(constant_term(newtonX[j])), order)
        postNewtonY[j, 1] = Taylor1(identity(constant_term(newtonY[j])), order)
        postNewtonZ[j, 1] = Taylor1(identity(constant_term(newtonZ[j])), order)
    end
    tmp2334 = Array{Taylor1{_S}}(undef, size(v2))
    tmp2334 .= Taylor1(zero(_S), order)
    tmp2335 = Array{Taylor1{_S}}(undef, size(v2))
    tmp2335 .= Taylor1(zero(_S), order)
    tmp2337 = Array{Taylor1{_S}}(undef, size(vi_dot_vj))
    tmp2337 .= Taylor1(zero(_S), order)
    Xij_t_Ui = Array{Taylor1{_S}}(undef, size(X))
    Xij_t_Ui .= Taylor1(zero(_S), order)
    Yij_t_Vi = Array{Taylor1{_S}}(undef, size(Y))
    Yij_t_Vi .= Taylor1(zero(_S), order)
    Zij_t_Wi = Array{Taylor1{_S}}(undef, size(Z))
    Zij_t_Wi .= Taylor1(zero(_S), order)
    tmp2343 = Array{Taylor1{_S}}(undef, size(Xij_t_Ui))
    tmp2343 .= Taylor1(zero(_S), order)
    Rij_dot_Vi = Array{Taylor1{_S}}(undef, size(tmp2343))
    Rij_dot_Vi .= Taylor1(zero(_S), order)
    tmp2346 = Array{Taylor1{_S}}(undef, size(Rij_dot_Vi))
    tmp2346 .= Taylor1(zero(_S), order)
    pn1t7 = Array{Taylor1{_S}}(undef, size(tmp2346))
    pn1t7 .= Taylor1(zero(_S), order)
    tmp2349 = Array{Taylor1{_S}}(undef, size(pn1t7))
    tmp2349 .= Taylor1(zero(_S), order)
    pn1t2_7 = Array{Taylor1{_S}}(undef, size(ϕs_and_vs))
    pn1t2_7 .= Taylor1(zero(_S), order)
    tmp2354 = Array{Taylor1{_S}}(undef, size(pNX_t_X))
    tmp2354 .= Taylor1(zero(_S), order)
    tmp2355 = Array{Taylor1{_S}}(undef, size(tmp2354))
    tmp2355 .= Taylor1(zero(_S), order)
    tmp2356 = Array{Taylor1{_S}}(undef, size(pn1t2_7))
    tmp2356 .= Taylor1(zero(_S), order)
    tmp2357 = Array{Taylor1{_S}}(undef, size(tmp2356))
    tmp2357 .= Taylor1(zero(_S), order)
    for k = Base.OneTo(postnewton_iter)
        for j = _1_to_N
            pntempX[j] = Taylor1(identity(constant_term(zero_q_1)), order)
            pntempY[j] = Taylor1(identity(constant_term(zero_q_1)), order)
            pntempZ[j] = Taylor1(identity(constant_term(zero_q_1)), order)
        end
        for j = _1_to_N
            for i = _1_to_N
                if i == j
                else
                    _4ϕj[i, j] = Taylor1(constant_term(4) * constant_term(newtonianNb_Potential[j]), order)
                    ϕi_plus_4ϕj[i, j] = Taylor1(constant_term(newtonianNb_Potential[i]) + constant_term(_4ϕj[i, j]), order)
                    tmp2334[i] = Taylor1(constant_term(2) * constant_term(v2[i]), order)
                    tmp2335[j] = Taylor1(constant_term(v2[j]) + constant_term(tmp2334[i]), order)
                    tmp2337[i, j] = Taylor1(constant_term(4) * constant_term(vi_dot_vj[i, j]), order)
                    sj2_plus_2si2_minus_4vivj[i, j] = Taylor1(constant_term(tmp2335[j]) - constant_term(tmp2337[i, j]), order)
                    ϕs_and_vs[i, j] = Taylor1(constant_term(sj2_plus_2si2_minus_4vivj[i, j]) - constant_term(ϕi_plus_4ϕj[i, j]), order)
                    Xij_t_Ui[i, j] = Taylor1(constant_term(X[i, j]) * constant_term(dq[3i - 2]), order)
                    Yij_t_Vi[i, j] = Taylor1(constant_term(Y[i, j]) * constant_term(dq[3i - 1]), order)
                    Zij_t_Wi[i, j] = Taylor1(constant_term(Z[i, j]) * constant_term(dq[3i]), order)
                    tmp2343[i, j] = Taylor1(constant_term(Xij_t_Ui[i, j]) + constant_term(Yij_t_Vi[i, j]), order)
                    Rij_dot_Vi[i, j] = Taylor1(constant_term(tmp2343[i, j]) + constant_term(Zij_t_Wi[i, j]), order)
                    tmp2346[i, j] = Taylor1(constant_term(Rij_dot_Vi[i, j]) ^ constant_term(2), order)
                    pn1t7[i, j] = Taylor1(constant_term(tmp2346[i, j]) / constant_term(r_p2[i, j]), order)
                    tmp2349[i, j] = Taylor1(constant_term(1.5) * constant_term(pn1t7[i, j]), order)
                    pn1t2_7[i, j] = Taylor1(constant_term(ϕs_and_vs[i, j]) - constant_term(tmp2349[i, j]), order)
                    pNX_t_X[i, j] = Taylor1(constant_term(postNewtonX[i, k]) * constant_term(X[i, j]), order)
                    pNY_t_Y[i, j] = Taylor1(constant_term(postNewtonY[i, k]) * constant_term(Y[i, j]), order)
                    pNZ_t_Z[i, j] = Taylor1(constant_term(postNewtonZ[i, k]) * constant_term(Z[i, j]), order)
                    tmp2354[i, j] = Taylor1(constant_term(pNX_t_X[i, j]) + constant_term(pNY_t_Y[i, j]), order)
                    tmp2355[i, j] = Taylor1(constant_term(tmp2354[i, j]) + constant_term(pNZ_t_Z[i, j]), order)
                    tmp2356[i, j] = Taylor1(constant_term(pn1t2_7[i, j]) + constant_term(tmp2355[i, j]), order)
                    tmp2357[i, j] = Taylor1(constant_term(c_p2) + constant_term(tmp2356[i, j]), order)
                    pn1[i, j] = Taylor1(constant_term(newtonianCoeff[i, j]) * constant_term(tmp2357[i, j]), order)
                    X_t_pn1[i, j] = Taylor1(constant_term(X[i, j]) * constant_term(pn1[i, j]), order)
                    Y_t_pn1[i, j] = Taylor1(constant_term(Y[i, j]) * constant_term(pn1[i, j]), order)
                    Z_t_pn1[i, j] = Taylor1(constant_term(Z[i, j]) * constant_term(pn1[i, j]), order)
                    pn3[i, j] = Taylor1(constant_term(3.5) * constant_term(newtonian1b_Potential[i, j]), order)
                    U_t_pn2[i, j] = Taylor1(constant_term(pn2[i, j]) * constant_term(U[i, j]), order)
                    pNX_t_pn3[i, j] = Taylor1(constant_term(postNewtonX[i, k]) * constant_term(pn3[i, j]), order)
                    pntempX[j] = Taylor1(constant_term(pntempX[j]) + constant_term(X_t_pn1[i, j] + (U_t_pn2[i, j] + pNX_t_pn3[i, j])), order)
                    V_t_pn2[i, j] = Taylor1(constant_term(pn2[i, j]) * constant_term(V[i, j]), order)
                    pNY_t_pn3[i, j] = Taylor1(constant_term(postNewtonY[i, k]) * constant_term(pn3[i, j]), order)
                    pntempY[j] = Taylor1(constant_term(pntempY[j]) + constant_term(Y_t_pn1[i, j] + (V_t_pn2[i, j] + pNY_t_pn3[i, j])), order)
                    W_t_pn2[i, j] = Taylor1(constant_term(pn2[i, j]) * constant_term(W[i, j]), order)
                    pNZ_t_pn3[i, j] = Taylor1(constant_term(postNewtonZ[i, k]) * constant_term(pn3[i, j]), order)
                    pntempZ[j] = Taylor1(constant_term(pntempZ[j]) + constant_term(Z_t_pn1[i, j] + (W_t_pn2[i, j] + pNZ_t_pn3[i, j])), order)
                end
            end
        end
        for j = _1_to_N
            postNewtonX[j, k + 1] = Taylor1(constant_term(pntempX[j]) * constant_term(c_m2), order)
            postNewtonY[j, k + 1] = Taylor1(constant_term(pntempY[j]) * constant_term(c_m2), order)
            postNewtonZ[j, k + 1] = Taylor1(constant_term(pntempZ[j]) * constant_term(c_m2), order)
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
                    TaylorSeries.identity!(pn1[i, j], zero_q_1, ord)
                    TaylorSeries.identity!(pn2[i, j], zero_q_1, ord)
                    TaylorSeries.identity!(pn3[i, j], zero_q_1, ord)
                    TaylorSeries.identity!(_4ϕj[i, j], zero_q_1, ord)
                    TaylorSeries.identity!(ϕi_plus_4ϕj[i, j], zero_q_1, ord)
                    TaylorSeries.identity!(sj2_plus_2si2_minus_4vivj[i, j], zero_q_1, ord)
                    TaylorSeries.identity!(ϕs_and_vs[i, j], zero_q_1, ord)
                    TaylorSeries.identity!(X_t_pn1[i, j], zero_q_1, ord)
                    TaylorSeries.identity!(Y_t_pn1[i, j], zero_q_1, ord)
                    TaylorSeries.identity!(Z_t_pn1[i, j], zero_q_1, ord)
                    TaylorSeries.identity!(UU[i, j], zero_q_1, ord)
                    TaylorSeries.identity!(VV[i, j], zero_q_1, ord)
                    TaylorSeries.identity!(WW[i, j], zero_q_1, ord)
                    TaylorSeries.identity!(U_t_pn2[i, j], zero_q_1, ord)
                    TaylorSeries.identity!(V_t_pn2[i, j], zero_q_1, ord)
                    TaylorSeries.identity!(W_t_pn2[i, j], zero_q_1, ord)
                    TaylorSeries.identity!(pNX_t_pn3[i, j], zero_q_1, ord)
                    TaylorSeries.identity!(pNY_t_pn3[i, j], zero_q_1, ord)
                    TaylorSeries.identity!(pNZ_t_pn3[i, j], zero_q_1, ord)
                    TaylorSeries.identity!(pNX_t_X[i, j], zero_q_1, ord)
                    TaylorSeries.identity!(pNY_t_Y[i, j], zero_q_1, ord)
                    TaylorSeries.identity!(pNZ_t_Z[i, j], zero_q_1, ord)
                end
            end
            TaylorSeries.identity!(vi_dot_vj[j], zero_q_1, ord)
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
                    TaylorSeries.mul!(tmp2184[3j - 2], 4, dq[3j - 2], ord)
                    TaylorSeries.mul!(tmp2186[3i - 2], 3, dq[3i - 2], ord)
                    TaylorSeries.subst!(_4U_m_3X[i, j], tmp2184[3j - 2], tmp2186[3i - 2], ord)
                    TaylorSeries.mul!(tmp2189[3j - 1], 4, dq[3j - 1], ord)
                    TaylorSeries.mul!(tmp2191[3i - 1], 3, dq[3i - 1], ord)
                    TaylorSeries.subst!(_4V_m_3Y[i, j], tmp2189[3j - 1], tmp2191[3i - 1], ord)
                    TaylorSeries.mul!(tmp2194[3j], 4, dq[3j], ord)
                    TaylorSeries.mul!(tmp2196[3i], 3, dq[3i], ord)
                    TaylorSeries.subst!(_4W_m_3Z[i, j], tmp2194[3j], tmp2196[3i], ord)
                    TaylorSeries.mul!(pn2x[i, j], X[i, j], _4U_m_3X[i, j], ord)
                    TaylorSeries.mul!(pn2y[i, j], Y[i, j], _4V_m_3Y[i, j], ord)
                    TaylorSeries.mul!(pn2z[i, j], Z[i, j], _4W_m_3Z[i, j], ord)
                    TaylorSeries.mul!(UU[i, j], dq[3i - 2], dq[3j - 2], ord)
                    TaylorSeries.mul!(VV[i, j], dq[3i - 1], dq[3j - 1], ord)
                    TaylorSeries.mul!(WW[i, j], dq[3i], dq[3j], ord)
                    TaylorSeries.add!(tmp2204[i, j], UU[i, j], VV[i, j], ord)
                    TaylorSeries.add!(vi_dot_vj[i, j], tmp2204[i, j], WW[i, j], ord)
                    TaylorSeries.pow!(tmp2207[i, j], X[i, j], 2, ord)
                    TaylorSeries.pow!(tmp2209[i, j], Y[i, j], 2, ord)
                    TaylorSeries.add!(tmp2210[i, j], tmp2207[i, j], tmp2209[i, j], ord)
                    TaylorSeries.pow!(tmp2212[i, j], Z[i, j], 2, ord)
                    TaylorSeries.add!(r_p2[i, j], tmp2210[i, j], tmp2212[i, j], ord)
                    TaylorSeries.sqrt!(r_p1d2[i, j], r_p2[i, j], ord)
                    TaylorSeries.pow!(r_p3d2[i, j], r_p2[i, j], 1.5, ord)
                    TaylorSeries.pow!(r_p7d2[i, j], r_p2[i, j], 3.5, ord)
                    TaylorSeries.div!(newtonianCoeff[i, j], μ[i], r_p3d2[i, j], ord)
                    TaylorSeries.add!(tmp2220[i, j], pn2x[i, j], pn2y[i, j], ord)
                    TaylorSeries.add!(tmp2221[i, j], tmp2220[i, j], pn2z[i, j], ord)
                    TaylorSeries.mul!(pn2[i, j], newtonianCoeff[i, j], tmp2221[i, j], ord)
                    TaylorSeries.add!(newtonX[j], newtonX[j], X[i, j] * newtonianCoeff[i, j], ord)
                    TaylorSeries.add!(newtonY[j], newtonY[j], Y[i, j] * newtonianCoeff[i, j], ord)
                    TaylorSeries.add!(newtonZ[j], newtonZ[j], Z[i, j] * newtonianCoeff[i, j], ord)
                    TaylorSeries.div!(newtonian1b_Potential[i, j], μ[i], r_p1d2[i, j], ord)
                    TaylorSeries.add!(newtonianNb_Potential[j], newtonianNb_Potential[j], newtonian1b_Potential[i, j], ord)
                    if UJ_interaction[i, j]
                        TaylorSeries.mul!(t31[i, j], X[i, j], M_[1, 3, j], ord)
                        TaylorSeries.mul!(t32[i, j], Y[i, j], M_[2, 3, j], ord)
                        TaylorSeries.mul!(t33[i, j], Z[i, j], M_[3, 3, j], ord)
                        TaylorSeries.add!(tmp2234[i, j], t31[i, j], t32[i, j], ord)
                        TaylorSeries.add!(r_sin_ϕ[i, j], tmp2234[i, j], t33[i, j], ord)
                        TaylorSeries.div!(sin_ϕ[i, j], r_sin_ϕ[i, j], r_p1d2[i, j], ord)
                        TaylorSeries.asin!(ϕ[i, j], sin_ϕ[i, j], tmp2385[i, j], ord)
                        TaylorSeries.sincos!(tmp2386[i, j], cos_ϕ[i, j], ϕ[i, j], ord)
                        TaylorSeries.pow!(sin2_ϕ[i, j], sin_ϕ[i, j], 2, ord)
                        TaylorSeries.pow!(sin3_ϕ[i, j], sin_ϕ[i, j], 3, ord)
                        TaylorSeries.mul!(tmp2244[i, j], 1.5, sin2_ϕ[i, j], ord)
                        TaylorSeries.subst!(P_2_sin_ϕ[i, j], tmp2244[i, j], 0.5, ord)
                        TaylorSeries.mul!(∂P_2_sin_ϕ[i, j], 3, sin_ϕ[i, j], ord)
                        TaylorSeries.mul!(tmp2250[i, j], -1.5, sin_ϕ[i, j], ord)
                        TaylorSeries.mul!(tmp2252[i, j], 2.5, sin3_ϕ[i, j], ord)
                        TaylorSeries.add!(P_3_sin_ϕ[i, j], tmp2250[i, j], tmp2252[i, j], ord)
                        TaylorSeries.mul!(tmp2256[i, j], 7.5, sin2_ϕ[i, j], ord)
                        TaylorSeries.add!(∂P_3_sin_ϕ[i, j], -1.5, tmp2256[i, j], ord)
                        TaylorSeries.subst!(tmp2258[j], Λ2[j], ord)
                        TaylorSeries.pow!(tmp2260[i, j], r_p2[i, j], 2, ord)
                        TaylorSeries.div!(Λ2j_div_r4[i, j], tmp2258[j], tmp2260[i, j], ord)
                        TaylorSeries.subst!(tmp2262[j], Λ3[j], ord)
                        TaylorSeries.pow!(tmp2264[i, j], r_p1d2[i, j], 5, ord)
                        TaylorSeries.div!(Λ3j_div_r5[i, j], tmp2262[j], tmp2264[i, j], ord)
                        TaylorSeries.subst!(tmp2266[i, j], cos_ϕ[i, j], ord)
                        TaylorSeries.mul!(m_c_ϕ_∂P_2[i, j], tmp2266[i, j], ∂P_2_sin_ϕ[i, j], ord)
                        TaylorSeries.subst!(tmp2268[i, j], cos_ϕ[i, j], ord)
                        TaylorSeries.mul!(m_c_ϕ_∂P_3[i, j], tmp2268[i, j], ∂P_3_sin_ϕ[i, j], ord)
                        TaylorSeries.mul!(tmp2271[i, j], Λ2j_div_r4[i, j], 3, ord)
                        TaylorSeries.mul!(F_J2_ξ[i, j], tmp2271[i, j], P_2_sin_ϕ[i, j], ord)
                        TaylorSeries.mul!(F_J2_ζ[i, j], Λ2j_div_r4[i, j], m_c_ϕ_∂P_2[i, j], ord)
                        TaylorSeries.mul!(tmp2275[i, j], Λ3j_div_r5[i, j], 4, ord)
                        TaylorSeries.mul!(F_J3_ξ[i, j], tmp2275[i, j], P_3_sin_ϕ[i, j], ord)
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
                        TaylorSeries.add!(accX[j], accX[j], μ[i] * F_J2_x[i, j], ord)
                        TaylorSeries.add!(accY[j], accY[j], μ[i] * F_J2_y[i, j], ord)
                        TaylorSeries.add!(accZ[j], accZ[j], μ[i] * F_J2_z[i, j], ord)
                        TaylorSeries.subst!(accX[i], accX[i], μ[j] * F_J2_x[i, j], ord)
                        TaylorSeries.subst!(accY[i], accY[i], μ[j] * F_J2_y[i, j], ord)
                        TaylorSeries.subst!(accZ[i], accZ[i], μ[j] * F_J2_z[i, j], ord)
                    end
                end
            end
            TaylorSeries.pow!(tmp2323[3j - 2], dq[3j - 2], 2, ord)
            TaylorSeries.pow!(tmp2325[3j - 1], dq[3j - 1], 2, ord)
            TaylorSeries.add!(tmp2326[3j - 2], tmp2323[3j - 2], tmp2325[3j - 1], ord)
            TaylorSeries.pow!(tmp2328[3j], dq[3j], 2, ord)
            TaylorSeries.add!(v2[j], tmp2326[3j - 2], tmp2328[3j], ord)
        end
        for j = _1_to_N
            TaylorSeries.identity!(postNewtonX[j, 1], newtonX[j], ord)
            TaylorSeries.identity!(postNewtonY[j, 1], newtonY[j], ord)
            TaylorSeries.identity!(postNewtonZ[j, 1], newtonZ[j], ord)
        end
        for k = Base.OneTo(postnewton_iter)
            for j = _1_to_N
                TaylorSeries.identity!(pntempX[j], zero_q_1, ord)
                TaylorSeries.identity!(pntempY[j], zero_q_1, ord)
                TaylorSeries.identity!(pntempZ[j], zero_q_1, ord)
            end
            for j = _1_to_N
                for i = _1_to_N
                    if i == j
                    else
                        TaylorSeries.mul!(_4ϕj[i, j], 4, newtonianNb_Potential[j], ord)
                        TaylorSeries.add!(ϕi_plus_4ϕj[i, j], newtonianNb_Potential[i], _4ϕj[i, j], ord)
                        TaylorSeries.mul!(tmp2334[i], 2, v2[i], ord)
                        TaylorSeries.add!(tmp2335[j], v2[j], tmp2334[i], ord)
                        TaylorSeries.mul!(tmp2337[i, j], 4, vi_dot_vj[i, j], ord)
                        TaylorSeries.subst!(sj2_plus_2si2_minus_4vivj[i, j], tmp2335[j], tmp2337[i, j], ord)
                        TaylorSeries.subst!(ϕs_and_vs[i, j], sj2_plus_2si2_minus_4vivj[i, j], ϕi_plus_4ϕj[i, j], ord)
                        TaylorSeries.mul!(Xij_t_Ui[i, j], X[i, j], dq[3i - 2], ord)
                        TaylorSeries.mul!(Yij_t_Vi[i, j], Y[i, j], dq[3i - 1], ord)
                        TaylorSeries.mul!(Zij_t_Wi[i, j], Z[i, j], dq[3i], ord)
                        TaylorSeries.add!(tmp2343[i, j], Xij_t_Ui[i, j], Yij_t_Vi[i, j], ord)
                        TaylorSeries.add!(Rij_dot_Vi[i, j], tmp2343[i, j], Zij_t_Wi[i, j], ord)
                        TaylorSeries.pow!(tmp2346[i, j], Rij_dot_Vi[i, j], 2, ord)
                        TaylorSeries.div!(pn1t7[i, j], tmp2346[i, j], r_p2[i, j], ord)
                        TaylorSeries.mul!(tmp2349[i, j], 1.5, pn1t7[i, j], ord)
                        TaylorSeries.subst!(pn1t2_7[i, j], ϕs_and_vs[i, j], tmp2349[i, j], ord)
                        TaylorSeries.mul!(pNX_t_X[i, j], postNewtonX[i, k], X[i, j], ord)
                        TaylorSeries.mul!(pNY_t_Y[i, j], postNewtonY[i, k], Y[i, j], ord)
                        TaylorSeries.mul!(pNZ_t_Z[i, j], postNewtonZ[i, k], Z[i, j], ord)
                        TaylorSeries.add!(tmp2354[i, j], pNX_t_X[i, j], pNY_t_Y[i, j], ord)
                        TaylorSeries.add!(tmp2355[i, j], tmp2354[i, j], pNZ_t_Z[i, j], ord)
                        TaylorSeries.add!(tmp2356[i, j], pn1t2_7[i, j], tmp2355[i, j], ord)
                        TaylorSeries.add!(tmp2357[i, j], c_p2, tmp2356[i, j], ord)
                        TaylorSeries.mul!(pn1[i, j], newtonianCoeff[i, j], tmp2357[i, j], ord)
                        TaylorSeries.mul!(X_t_pn1[i, j], X[i, j], pn1[i, j], ord)
                        TaylorSeries.mul!(Y_t_pn1[i, j], Y[i, j], pn1[i, j], ord)
                        TaylorSeries.mul!(Z_t_pn1[i, j], Z[i, j], pn1[i, j], ord)
                        TaylorSeries.mul!(pn3[i, j], 3.5, newtonian1b_Potential[i, j], ord)
                        TaylorSeries.mul!(U_t_pn2[i, j], pn2[i, j], U[i, j], ord)
                        TaylorSeries.mul!(pNX_t_pn3[i, j], postNewtonX[i, k], pn3[i, j], ord)
                        TaylorSeries.add!(pntempX[j], pntempX[j], X_t_pn1[i, j] + (U_t_pn2[i, j] + pNX_t_pn3[i, j]), ord)
                        TaylorSeries.mul!(V_t_pn2[i, j], pn2[i, j], V[i, j], ord)
                        TaylorSeries.mul!(pNY_t_pn3[i, j], postNewtonY[i, k], pn3[i, j], ord)
                        TaylorSeries.add!(pntempY[j], pntempY[j], Y_t_pn1[i, j] + (V_t_pn2[i, j] + pNY_t_pn3[i, j]), ord)
                        TaylorSeries.mul!(W_t_pn2[i, j], pn2[i, j], W[i, j], ord)
                        TaylorSeries.mul!(pNZ_t_pn3[i, j], postNewtonZ[i, k], pn3[i, j], ord)
                        TaylorSeries.add!(pntempZ[j], pntempZ[j], Z_t_pn1[i, j] + (W_t_pn2[i, j] + pNZ_t_pn3[i, j]), ord)
                    end
                end
            end
            for j = _1_to_N
                TaylorSeries.mul!(postNewtonX[j, k + 1], pntempX[j], c_m2, ord)
                TaylorSeries.mul!(postNewtonY[j, k + 1], pntempY[j], c_m2, ord)
                TaylorSeries.mul!(postNewtonZ[j, k + 1], pntempZ[j], c_m2, ord)
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