@doc raw"""
    semimajoraxis(x, y, z, u, v, w, m1, m2)

Calculate semimajor axis for the two body problem defined by the relative position 
`(x,y,z)` and velocity `(u,v,w)` vectors between two bodies with masses `m1` and `m2`.
"""
function semimajoraxis(x, y, z, u, v, w, m1, m2)
    r = sqrt(x^2+y^2+z^2)       # Position magnitude
    vsq = u^2+v^2+w^2           # Velocity magnitude squared
    M = m1+m2                   # Total mass
    return 1/( (2/r)-(vsq/M) )  # Semimajor axis 
end

@doc raw"""
    eccentricity(x, y, z, u, v, w, m1, m2)

Calculate eccentricity for the two body problem defined by the relative position 
`(x,y,z)` and velocity `(u,v,w)` vectors between two bodies with masses `m1` and `m2`.
"""
function eccentricity(x, y, z, u, v, w, m1, m2)
    # h: Angular momentum per unit mass
    hsq = (y*w-z*v)^2 + (z*u-x*w)^2 + (x*v-y*u)^2   # h magnitude squared
    a = semimajoraxis(x, y, z, u, v, w, m1, m2)     # Semimajor axis
    M = m1+m2                                       # Total mass   
    quotient = hsq/( M*a )                          
    return sqrt(1 - quotient)                       # Eccentricity 
end

@doc raw"""
    inclination(x, y, z, u, v, w)

Calculate inclination for the two body problem defined by the relative position 
`(x,y,z)` and velocity `(u,v,w)` vectors between two bodies."""
function inclination(x, y, z, u, v, w)
    # h: Angular momentum per unit mass
    hz = x*v-y*u                                         # z-coordinate of h
    h = sqrt( (y*w-z*v)^2 + (z*u-x*w)^2 + (x*v-y*u)^2 )  # h magnitude squared
    return acos(hz/h)                                    # Inclination 
end

@doc raw"""
    aei(x, y, z, u, v, w, m1, m2)

Return semimajor axis `a`, eccentricity `e` and inclination `inc` for the two-body
problem defined by the relative position `(x,y,z)` and velocity `(u,v,w)`vectors
between two bodies with masses `m1` and `m2`.

See also [`semimajoraxis`](@ref), [`eccentricity`](@ref) and [`inclination`](@ref).
"""
function aei(x, y, z, u, v, w, m1, m2)
    r = sqrt(x^2+y^2+z^2)                           # Position magnitude
    vsq = u^2+v^2+w^2                               # Velocity magnitude squared
    M = m1+m2                                       # Total mass 
    a = 1/( (2/r)-(vsq/M) )                         # Semimajor axis
    # h: Angular momentum per unit mass
    hsq = (y*w-z*v)^2 + (z*u-x*w)^2 + (x*v-y*u)^2   # h magnitude squared
    h = sqrt(hsq)                                   # h magnitude
    quotient = hsq/( M*a )         
    e = sqrt(1 - quotient)                          # Eccentricity
    hz = x*v-y*u                                    # z-coordinate of h
    inc = acos(hz/h)                                # Inclination 

    return a, e, inc
end

@doc """
    ae(x, y, z, u, v, w)

Return semimajor axis `a` and eccentricity `e` for the two-body problem defined
by the relative position `(x,y,z)` and velocity `(u,v,w)` vectors between two
bodies with masses `m1` and `m2`.

See also [`semimajoraxis`](@ref) and [`eccentricity`](@ref).
"""
function ae(x, y, z, u, v, w, m1, m2)
    r = sqrt(x^2+y^2+z^2)                           # Position magnitude
    vsq = u^2+v^2+w^2                               # Velocity magnitude squared
    M = m1+m2                                       # Total mass
    a = 1/( (2/r)-(vsq/M) )                         # Semimajor axis 
    # h: Angular momentum per unit mass
    hsq = (y*w-z*v)^2 + (z*u-x*w)^2 + (x*v-y*u)^2   # h magnitude squared
    quotient = hsq/( M*a )                  
    e = sqrt(1 - quotient)                          # Eccentricity
    return a, e
end

@doc raw"""
    longascnode(x, y, z, u, v, w)

Calculate longitude of ascending node for the two body problem defined by the relative 
position `(x,y,z)` and velocity `(u,v,w)` vectors between two bodies.
"""
function longascnode(x, y, z, u, v, w)

    # The longitude of ascending node is computed as the angle between x-axis and the 
    # vector n = (-hy, hx, 0) where hx, hy, are resp., the x and y comps. of angular
    # momentum vector per unit mass, h

    res = atan( y*w-z*v, x*w-z*u)

    if constant_term(res) ≥ zero(constant_term(res))
        return res
    else
        return res+2pi
    end

end

@doc raw"""
    argperi(x, y, z, u, v, w, m1, m2)

Calculate the instantaneous (osculating) argument of pericentre for the two
body problem defined by the relative position `(x,y,z)` and velocity `(u,v,w)` vectors 
between two bodies with masses `m1` and `m2`.
"""
function argperi(x, y, z, u, v, w, m1, m2)
    # h: Angular momentum per unit mass
    # n = (z-axis unit vector)×h = (-hy, hx, 0)
    n = [x*w-z*u, y*w-z*v, zero(x)]
    e = lrlvec(x,y,z,u,v,w,m1,m2)       # Cartesian comps. of Laplace-Runge-Lenz vector
    n = n/sqrt(n[1]^2+n[2]^2+n[3]^2)    # n unit vector
    e = e/sqrt(e[1]^2+e[2]^2+e[3]^2)    # Laplace-Runge-Lenz unit vector
    cosω = dot(n, e)                    # Cosine of argument of pericentre ω

    if constant_term(e[3]) >= zero(constant_term(x))
        return acos(cosω)
    else
        return 2pi-acos(cosω)
    end

end

@doc raw"""
    longperi(x, y, z, u, v, w, m1, m2)

Calculate the instantaneous (osculating) longitude of pericentre for the two
body problem defined by the relative position `(x,y,z)` and velocity `(u,v,w)` vectors 
between two bodies with masses `m1` and `m2`.
"""
function longperi(x, y, z, u, v, w, m1, m2)
    # Argument of pericentre + longitude of ascending node
    return argperi(x, y, z, u, v, w, m1, m2)+longascnode(x, y, z, u, v, w)
end

@doc raw"""
    trueanomaly(x, y, z, u, v, w, m1, m2)

Calculate the instantaneous (osculating) true anomaly for the two
body problem defined by the relative position `(x,y,z)` and velocity `(u,v,w)` vectors 
between two bodies with masses `m1` and `m2`.
"""
function trueanomaly(x, y, z, u, v, w, m1, m2)

    R = sqrt(x^2+y^2+z^2)    # Position magnitude
    V2 = u^2+v^2+w^2         # Velocity mangitude squared
    V = sqrt(V2)             # Velocity magnitude 

    a = 1.0/( (2.0/R)-(V2/(m1+m2)) )   # Semimajor axis

    # h: Angular momentum per unit mass
    h2 = (y*w-z*v)^2 + (z*u-x*w)^2 + (x*v-y*u)^2   # h magnitude squared 
    h = sqrt(h2)                                   # h magnitude  

    quotient=h2/( (m1+m2)*a )

    e = sqrt(1.0 - quotient)  # Eccentricity

    Rdot = (x*u+y*v+z*w)/R    # <Position, Velocity> / Position magnitude

    sin_f = Rdot*a*(1.0-e^2)/( e*h )         # sin( true anomaly )
    cos_f = (1.0/e)*( -1.0+a*(1.0-e^2)/R )   # cos( true anomaly )

    res = atan(sin_f, cos_f)  # True anomaly

    if res ≥ zero(res)
        return res
    else
        return res+2pi
    end
end

@doc raw"""
    ecanomaly(x, y, z, u, v, w, m1, m2)

Calculate the instantaneous (osculating) eccentric anomaly for the two
body problem defined by the relative position `(x,y,z)` and velocity `(u,v,w)` vectors 
between two bodies with masses `m1` and `m2`.

See also [`semimajoraxis`](@ref) and [`eccentricity`](@ref).
"""
function ecanomaly(x, y, z, u, v, w, m1, m2)

    a = semimajoraxis(x, y, z, u, v, w, m1, m2)  # Semimajor axis
    e = eccentricity(x, y, z, u, v, w, m1, m2)   # Eccentricity
    R = sqrt(x^2+y^2+z^2)                        # Position magnitude
    Rdot = (x*u+y*v+z*w)/R                       # <Position, Velocity> / Position magnitude
    mu = m1+m2                                   # Total mass

    cosE = (a-R)/(a*e)                           # sin( eccentric anomaly )
    sinE = Rdot*R/( e*sqrt(mu*a) )               # cos( eccentric anomaly )

    res = atan(sinE, cosE)                       # Eccentric anomaly 

    return res

end

@doc """
    meananomaly(x, y, z, u, v, w, m1, m2)

Calculate the instantaneous (osculating) mean anomaly for the two 
body problem defined by the relative position `(x,y,z)` and velocity `(u,v,w)` vectors 
between two bodies with masses `m1` and `m2`.

See also [`ecanomaly`](@ref) and [`eccentricity`](@ref).
"""
function meananomaly(x, y, z, u, v, w, m1, m2)

    E = ecanomaly(x, y, z, u, v, w, m1, m2)    # Eccentric anomaly
    e = eccentricity(x, y, z, u, v, w, m1, m2) # Eccentricity
    # Mean anomaly
    return E-e*sin(E)
end

@doc """
    timeperipass(t, x, y, z, u, v, w, m1, m2)

Calculate the instantaneous (osculating) time of pericentre passage, at time `t`, for the two
body problem defined by the relative position `(x,y,z)` and velocity `(u,v,w)` vectors between
two bodies with masses `m1` and `m2`.

See also [`meananomaly`](@ref), [`semimajoraxis`](@ref) and [`meanmotion`](@ref).
"""
function timeperipass(t, x, y, z, u, v, w, m1, m2)

    me_an = meananomaly(x, y, z, u, v, w, m1, m2) # Mean anomaly
    a = semimajoraxis(x, y, z, u, v, w, m1, m2)   # Semimajor axis
    me_mo = meanmotion(m1+m2, a)                  # Mean motion 
    # Time of pericentre passage
    return t-me_an/me_mo
end

@doc """
    rungelenzx(x, y, z, u, v, w, m1, m2)

Calculate the instantaneous (osculating) x-component of the Laplace-Runge-Lenz vector for the
two body problem defined by the relative position `(x,y,z)` and velocity `(u,v,w)` vectors 
between two bodies with masses `m1` and `m2`.

See also [`rungelenzy`](@ref), [`rungelenzz`](@ref), [`rungelenzmag`](@ref) and [`lrlvec`](@ref).
"""
function rungelenzx(x, y, z, u, v, w, m1, m2)

    # Velocity v = [u, v, w]
    # Angular momentum per unit mass h = [y.*w-z.*v, z.*u-x.*w, x.*v-y.*u]
    # Total mass
    mu = m1+m2                                         
    # Position r = [x, y, z]
    rmag = sqrt(x^2+y^2+z^2)  # Position magnitude
    # x-component of Laplace-Runge-Lenz vector
    lrl_x = ( -(z*u-x*w)*w+(x*v-y*u)*v )/(mu)-x/rmag

    return lrl_x
end

@doc raw"""
    rungelenzy(x, y, z, u, v, w, m1, m2)

Calculate the instantaneous (osculating) y-component of the Laplace-Runge-Lenz vector for the
two body problem defined by the relative position `(x,y,z)` and velocity `(u,v,w)` vectors 
between two bodies with masses `m1` and `m2`.

See also [`rungelenzx`](@ref), [`rungelenzz`](@ref), [`rungelenzmag`](@ref) and [`lrlvec`](@ref).
"""
function rungelenzy(x, y, z, u, v, w, m1, m2)

    # Velocity v = [u, v, w]
    # Angular momentum per unit mass h = [y.*w-z.*v, z.*u-x.*w, x.*v-y.*u]
    # Total mass
    mu = m1+m2
    # Position r = [x, y, z]
    rmag = sqrt(x^2+y^2+z^2)  # Position magnitude
    # y-component of Laplace-Runge-Lenz vector
    lrl_y = ( -(x*v-y*u)*u+(y*w-z*v)*w )/(mu)-y/rmag

    return lrl_y
end

@doc raw"""
    rungelenzz(x, y, z, u, v, w, m1, m2)

Calculate the instantaneous (osculating) z-component of the Laplace-Runge-Lenz vector for the
two body problem defined by the relative position `(x,y,z)` and velocity `(u,v,w)` vectors 
between two bodies with masses `m1` and `m2`.

See also [`rungelenzx`](@ref), [`rungelenzy`](@ref), [`rungelenzmag`](@ref) and [`lrlvec`](@ref).
"""
function rungelenzz(x, y, z, u, v, w, m1, m2)

    # Velocity v = [u, v, w]
    # Angular momentum per unit mass h = [y.*w-z.*v, z.*u-x.*w, x.*v-y.*u]
    # Total mass
    mu = m1+m2
    # Position r = [x, y, z]
    rmag = sqrt(x^2+y^2+z^2)  # Position magnitude
    # z-component of Laplace-Runge-Lenz vector
    lrl_z = ( -(y*w-z*v)*v+(z*u-x*w)*u )/(mu)-z/rmag

    return lrl_z
end

@doc raw"""
    rungelenzmag(x, y, z, u, v, w, m1, m2)

Calculate the instantaneous (osculating) magnitude of the Laplace-Runge-Lenz vector for the
two body problem defined by the relative position `(x,y,z)` and velocity `(u,v,w)` vectors 
between two bodies with masses `m1` and `m2`.

See also [`rungelenzx`](@ref), [`rungelenzy`](@ref), [`rungelenzz`](@ref) and [`lrlvec`](@ref).
"""
function rungelenzmag(x, y, z, u, v, w, m1, m2)

    # Velocity v = [u, v, w]
    # Angular momentum per unit mass h = [y.*w-z.*v, z.*u-x.*w, x.*v-y.*u]
    # Total mass
    mu = m1+m2
    # Position r = [x, y, z]
    rmag = sqrt(x^2+y^2+z^2)  # Position magnitude
    
    # Components of Laplace-Runge-Lenz vector
    lrl_x = ( -(z*u-x*w)*w+(x*v-y*u)*v )/(mu)-x/rmag  # x-coordinate
    lrl_y = ( -(x*v-y*u)*u+(y*w-z*v)*w )/(mu)-y/rmag  # y-coordinate
    lrl_z = ( -(y*w-z*v)*v+(z*u-x*w)*u )/(mu)-z/rmag  # z-coordinate
    # Mangnitude of Laplace-Runge-Lenz vector
    return sqrt(lrl_x^2+lrl_y^2+lrl_z^2)
end

@doc raw"""
    lrlvec(x, y, z, u, v, w, m1, m2)

Calculate the instantaneous (osculating) cartesian components of the Laplace-Runge-Lenz vector
for the two body problem defined by the relative position `(x,y,z)` and velocity `(u,v,w)` 
vectors between two bodies with masses `m1` and `m2`.

See also [`rungelenzx`](@ref), [`rungelenzy`](@ref), [`rungelenzz`](@ref) and [`rungelenzmag`](@ref).
"""
function lrlvec(x, y, z, u, v, w, m1, m2)
    # Total mass
    mu = m1+m2
    # Position magnitude
    rmag = sqrt(x^2+y^2+z^2)
    # Components of Laplace-Runge-Lenz vector
    lrl_x = ( -(z*u-x*w)*w+(x*v-y*u)*v )/(mu)-x/rmag  # x-coordinate
    lrl_y = ( -(x*v-y*u)*u+(y*w-z*v)*w )/(mu)-y/rmag  # y-coordinate
    lrl_z = ( -(y*w-z*v)*v+(z*u-x*w)*u )/(mu)-z/rmag  # z-coordinate
    # Laplace-Runge-Lenz vector
    return [lrl_x, lrl_y, lrl_z]
end

@doc raw"""
    eccentricanomaly(e::T, M::T) where {T <: Number}

Compute eccentric anomaly (`E`) from eccentricity (`e`) and mean anomaly (`M`). 
See pages 32-37 of https://doi.org/10.1017/CBO9781139174817.
"""
function eccentricanomaly(e::T, M::T) where {T <: Number}
    M0 = mod2pi(M)
    E0 = M0 + sign(sin(M0))*0.85*e # Murray-Dermotts' initial estimate
    # Successive approximations via Newtons' method
    for i in 0:4
        # TODO: implement modified Newton's method for Kepler's equation (Murray-Dermott)
        Eans = E0 - (E0-e*sin(E0)-M0)/(1.0-e*cos(E0))
        E0 = Eans
    end
    return E0
end

@doc raw"""
    trueanomaly(e,E)

Compute true anomaly (`f`) from eccentricity (`e`) and eccentric anomaly (`E`).
"""
function trueanomaly(e,E)
    Enew = mod2pi(E)
    return 2.0*atan(  sqrt((1.0+e)/(1.0-e))*tan(Enew/2)  )
end

@doc raw"""
    meanan2truean(e,M)

Compute true anomaly (`f`) from eccentricity (`e`) and mean anomaly (`M`).

See also [`trueanomaly`](@ref) and [`eccentricanomaly`](@ref).
"""
function meanan2truean(e,M)
    return trueanomaly(e, eccentricanomaly(e, M))
end

@doc raw"""
    meanmotion(mu,a)

Compute mean motion from mass parameter (`mu`) and semimajor axis (`a`).
"""
function meanmotion(mu,a)
    return sqrt(mu/(a^3))
end

@doc raw"""
    meananomaly(n, t, taup)

Compute mean anomaly from mean motion (`n`), time (`t`) and time of pericenter passage (`taup`).
"""
function meananomaly(n, t, taup)
    return n*(t-taup)
end

@doc raw"""
    time2truean(a, e, mu, t, taup)

Compute true anomaly from time (`t`), semimajor axis (`a`), eccentricity (`e`), 
mass parameter (`mu`) and time of pericenter passage (`taup`).

See also [`meanan2truean`](@ref), [`meananomaly`](@ref) and [`meanmotion`](@ref).
"""
function time2truean(a, e, mu, t, taup)
    return meanan2truean(e, meananomaly(meanmotion(mu, a), t, taup))
end
