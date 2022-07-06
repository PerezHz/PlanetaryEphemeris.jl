@doc raw"""
    μ_star_fun(μ, q, i)

Returns the mass parameter of the i-th body with relativistic corrections up to ``1/c^2``
```math
\mu_i^* = \mu_i\left(1 + \frac{v_i^2}{2c^2} - \frac{1}{2c^2}\sum_{j\neq i}\frac{\mu_j}{r_{ij}}\right),
```
and its time derivative
```math
\dot{\mu}_i^* = \frac{\mu_i}{2c^2}\left(\sum_{j\neq i}\mu_j\frac{\left(\mathbf{r}_j - \mathbf{r}_i\right)\cdot\left(\dot{\mathbf{r}}_j + \dot{\mathbf{r}_i}\right)}{r_{ij}^3}\right).
```
See equations (2) and (4) of https://doi.org/10.1051/0004-6361:20066607 
"""
function μ_star_fun(μ, q, i)
    N = (length(q)-13) ÷ 6 # Total number of bodies (Sun+planets+Moon+Pluto+asts)
    sum1 = zero(q[1])      # Sum in the mass parameter 
    sum2 = zero(q[1])      # Sum in the time derivative of the mass parameter
    for j in 1:N           # Iterate over the N bodies...
        i ==j && continue  # except for the i-th body
        #@show j
        ith_body_ind = nbodyind(N, i)  # Index of the i-th body in q
        jth_body_ind = nbodyind(N, j)  # Index of the j-th body in q
        rvec_i = q[ith_body_ind[1:3]]  # Position vector of the i-th body
        rvec_j = q[jth_body_ind[1:3]]  # Position vector of the j-th body
        rvec_ij = rvec_i - rvec_j      # Vector joining the i-th and the j-th body
        r_ij = sqrt(rvec_ij[1]^2 + rvec_ij[2]^2 + rvec_ij[3]^2) # Distance between the i-th and the j-th body
        sum1 += μ[j]/r_ij              
        vvec_i = q[ith_body_ind[4:6]]  # Velocity vector of the i-th body
        vvec_j = q[jth_body_ind[4:6]]  # Velocity vector of the j-th body
        vvec_ipj = vvec_i + vvec_j     # Sum of the velocity vectors of the i-th and j-th body
        sum2 -= μ[j]*( dot(rvec_ij, vvec_ipj) )/(r_ij^3)
    end
    #@show sum1, sum2
    #rvec_i = q[3i-2:3i]
    vvec_i = q[3(N+i)-2:3(N+i)]  # Velocity vector of the i-th body
    v_i2 = dot(vvec_i, vvec_i)   # Speed of the i-th body
    mustar = μ[i]*( 1+v_i2/(2c_au_per_day^2)-sum1/(2c_au_per_day^2) ) # Mass parameter 
    mustardot = μ[i]*sum2/(2c_au_per_day^2)                           # Time derivative of the mass parameter  
    #@show μ[i]
    return mustar, mustardot
end

@doc raw"""
    ssb_posvel_pN(μ, q)

Returns the two sums needed to determine the Solar System Barycenter (SSB)
```math
\left\lbrace
\begin{array}{l}
    \sum_i \mu_i^*\mathbf{r}_i = \mathbf{0} \\
    \sum_i \mu_i^*\dot{\mathbf{r}}_i + \dot{\mu}_i^*\mathbf{r}_i = \mathbf{0}
\end{array}
\right.
```
and the sum of the mass parameters of all the bodies in `q`.

See LHS of equation (5) of https://doi.org/10.1051/0004-6361:20066607 
"""
function ssb_posvel_pN(μ, q)
    N = (length(q)-13) ÷ 6   # Total number of bodies (Sun+planets+Moon+Pluto+asts)
    rvec_ssb = zeros(3)      # SSB position
    vvec_ssb = zeros(3)      # SSB velocity
    μ_star_SSB = zero(q[1])  # Sum of the mass parameters of all the bodies in q 
    for i in 1:N             # Iterate over the N bodies  
        μ_star_i, μ_star_dot_i = μ_star_fun(μ, q, i)  # Mass parameter and its time derivative of the i-th body
        ith_body_ind = nbodyind(N, i)                 # Index of the i-th body
        rvec_i = q[ith_body_ind[1:3]]                 # Position vector of the i-th body
        vvec_i = q[ith_body_ind[4:6]]                 # Velocity vector of the i-th body  
        rvec_ssb += μ_star_i*rvec_i                   # Add to the SSB position
        vvec_ssb += μ_star_i*vvec_i + μ_star_dot_i*rvec_i # Add to the SSB velocity
        μ_star_SSB += μ_star_i                        # Add to the sum of mass parameters
    end
    return rvec_ssb, vvec_ssb, μ_star_SSB
end

@doc raw"""
    ssb_posvel_pN(μ, q)

Solves
```math
\left\lbrace
\begin{array}{l}
    \sum_i \mu_i^*\mathbf{r}_i = \mathbf{0} \\
    \sum_i \mu_i^*\dot{\mathbf{r}}_i + \dot{\mu}_i^*\mathbf{r}_i = \mathbf{0}
\end{array}
\right.
```
for Sun position and velocity.

See equation (5) of https://doi.org/10.1051/0004-6361:20066607 
"""
function sun_posvel_pN(μ, q)
    N = (length(q)-13) ÷ 6       # Total number of bodies (Sun+planets+Moon+Pluto+asts)
    rvec_ss_wout_sun = zeros(3)  # SSB position (without Sun's position)
    vvec_ss_wout_sun = zeros(3)  # SSB velocity (without Sun's velocity)
    for i in 2:N                 # Iterate over the N bodies except for the Sun (by default the first body)
        μ_star_i, μ_star_dot_i = μ_star_fun(μ, q, i)  # Mass parameter and its time derivative of the i-th body
        ith_body_ind = nbodyind(N, i)                 # Index of the i-th body
        rvec_i = q[ith_body_ind[1:3]]                 # Position vector of the i-th body
        vvec_i = q[ith_body_ind[4:6]]                 # Velocity vector of the i-th body
        rvec_ss_wout_sun += μ_star_i*rvec_i           # Add to SSB position 
        vvec_ss_wout_sun += μ_star_i*vvec_i + μ_star_dot_i*rvec_i  # Add do SSB velocity
    end
    μ_star_1, μ_star_dot_1 = μ_star_fun(μ, q, 1)  # Mass parameter and its time derivative of the Sun
    rvec_1 = q[1:3]                               # Position vector of the Sun 
    vvec_ss_wout_sun += μ_star_dot_1*rvec_1       # Add Sun to the SSB velocity 
    return -rvec_ss_wout_sun/μ_star_1, -vvec_ss_wout_sun/μ_star_1  # Solve for Sun position and velocity
end