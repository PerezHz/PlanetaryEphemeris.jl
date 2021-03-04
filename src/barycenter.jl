#Fienga et al. (A&A, 2006), Eqs. (2) and (4)
function μ_star_fun(μ, q, i)
    N = (length(q)-13) ÷ 6 # total number of bodies (Sun+planets+Moon+Pluto+asts)
    sum1 = zero(q[1])
    sum2 = zero(q[1])
    for j in 1:N
        i ==j && continue
        #@show j
        ith_body_ind = nbodyind(N, i)
        jth_body_ind = nbodyind(N, j)
        rvec_i = q[ith_body_ind[1:3]]
        rvec_j = q[jth_body_ind[1:3]]
        rvec_ij = rvec_i - rvec_j
        r_ij = sqrt(rvec_ij[1]^2 + rvec_ij[2]^2 + rvec_ij[3]^2)
        sum1 += μ[j]/r_ij
        vvec_i = q[ith_body_ind[4:6]]
        vvec_j = q[jth_body_ind[4:6]]
        vvec_ipj = vvec_i + vvec_j
        sum2 -= μ[j]*( dot(rvec_ij, vvec_ipj) )/(r_ij^3)
    end
    #@show sum1, sum2
    #rvec_i = q[3i-2:3i]
    vvec_i = q[3(N+i)-2:3(N+i)]
    v_i2 = dot(vvec_i, vvec_i)
    mustar = μ[i]*( 1+v_i2/(2c_au_per_day^2)-sum1/(2c_au_per_day^2) )
    mustardot = μ[i]*sum2/(2c_au_per_day^2)
    #@show μ[i]
    return mustar, mustardot
end

#Fienga et al. (A&A, 2006), LHS of Eq. (5)
function ssb_posvel_pN(μ, q)
    N = (length(q)-13) ÷ 6 # total number of bodies (Sun+planets+Moon+Pluto+asts)
    rvec_ssb = zeros(3)
    vvec_ssb = zeros(3)
    μ_star_SSB = zero(q[1])
    for i in 1:N
        μ_star_i, μ_star_dot_i = μ_star_fun(μ, q, i)
        ith_body_ind = nbodyind(N, i)
        rvec_i = q[ith_body_ind[1:3]]
        vvec_i = q[ith_body_ind[4:6]]
        rvec_ssb += μ_star_i*rvec_i
        vvec_ssb += μ_star_i*vvec_i + μ_star_dot_i*rvec_i
        μ_star_SSB += μ_star_i
    end
    return rvec_ssb, vvec_ssb, μ_star_SSB
end

#Fienga et al. (A&A, 2006), Eq. (5), solved for Sun position and velocity
function sun_posvel_pN(μ, q)
    N = (length(q)-13) ÷ 6 # total number of bodies (Sun+planets+Moon+Pluto+asts)
    rvec_ss_wout_sun = zeros(3)
    vvec_ss_wout_sun = zeros(3)
    for i in 2:N
        μ_star_i, μ_star_dot_i = μ_star_fun(μ, q, i)
        ith_body_ind = nbodyind(N, i)
        rvec_i = q[ith_body_ind[1:3]]
        vvec_i = q[ith_body_ind[4:6]]
        rvec_ss_wout_sun += μ_star_i*rvec_i
        vvec_ss_wout_sun += μ_star_i*vvec_i + μ_star_dot_i*rvec_i
    end
    μ_star_1, μ_star_dot_1 = μ_star_fun(μ, q, 1)
    rvec_1 = q[1:3]
    vvec_ss_wout_sun += μ_star_dot_1*rvec_1
    return -rvec_ss_wout_sun/μ_star_1, -vvec_ss_wout_sun/μ_star_1
end