"""`semimajoraxis(x, y, z, u, v, w, m1, m2)`

Calculates semimajor axis for the two
body problem defined by the relative position `(x,y,z)` and velocity `(u,v,w)` vectors between two bodies
with masses `m1` and `m2`."""
function semimajoraxis(x, y, z, u, v, w, m1, m2)
    r = sqrt(x^2+y^2+z^2)
    vsq = u^2+v^2+w^2
    M = m1+m2
    return 1/( (2/r)-(vsq/M) )
end

"""`eccentricity(x, y, z, u, v, w, m1, m2)`

Calculates eccentricity for the two
body problem defined by the relative position `(x,y,z)` and velocity `(u,v,w)` vectors between two bodies
with masses `m1` and `m2`."""
function eccentricity(x, y, z, u, v, w, m1, m2)
    hsq = (y*w-z*v)^2 + (z*u-x*w)^2 + (x*v-y*u)^2
    a = semimajoraxis(x, y, z, u, v, w, m1, m2)
    M = m1+m2
    quotient = hsq/( M*a )
    return sqrt(1 - quotient)
end

"""`inclination(x, y, z, u, v, w)`

Calculates inclination for the two
body problem defined by the relative position `(x,y,z)` and velocity `(u,v,w)` vectors between two bodies."""
function inclination(x, y, z, u, v, w)
    hz = x*v-y*u
    h = sqrt( (y*w-z*v)^2 + (z*u-x*w)^2 + (x*v-y*u)^2 )
    return acos(hz/h)
end

"""`aei(x, y, z, u, v, w, m1, m2)`

Returns semimajor axis `a`, eccentricity `e` and inclination `I` for the two-body
problem defined by the relative position `(x,y,z)` and velocity `(u,v,w)`vectors
between two bodies with masses `m1` and `m2`."""
function aei(x, y, z, u, v, w, m1, m2)
    r = sqrt(x^2+y^2+z^2)
    vsq = u^2+v^2+w^2
    M = m1+m2
    a = 1/( (2/r)-(vsq/M) )
    hsq = (y*w-z*v)^2 + (z*u-x*w)^2 + (x*v-y*u)^2
    h = sqrt(hsq)
    quotient = hsq/( M*a )
    e = sqrt(1 - quotient)
    hz = x*v-y*u
    inc = acos(hz/h)
    return a, e, inc
end

"""`ae(x, y, z, u, v, w)`

Returns semimajor axis `a` and eccentricity `e` for the two-body problem defined
by the relative position `(x,y,z)` and velocity `(u,v,w)` vectors between two
bodies with masses `m1` and `m2`."""
function ae(x, y, z, u, v, w, m1, m2)
    r = sqrt(x^2+y^2+z^2)
    vsq = u^2+v^2+w^2
    M = m1+m2
    a = 1/( (2/r)-(vsq/M) )
    hsq = (y*w-z*v)^2 + (z*u-x*w)^2 + (x*v-y*u)^2
    quotient = hsq/( M*a )
    e = sqrt(1 - quotient)
    return a, e
end

"""`longascnode(x, y, z, u, v, w)`

Calculates longitude of ascending node for the two
body problem defined by the relative position `(x,y,z)` and velocity `(u,v,w)` vectors between two bodies."""
function longascnode(
    x, y, z,
    u, v, w
    )

    #the longitude of ascending node is computed as
    #the angle between x-axis and the vector n = (-hy,hx,0)
    #where hx, hy, are resp., the x and y comps. of ang mom vector per unit mass, h

    res = atan( y*w-z*v, x*w-z*u)

    if constant_term(res) ≥ zero(constant_term(res))
        return res
    else
        return res+2pi
    end

end

"""`argperi(x, y, z, u, v, w, m1, m2)`

Calculates the instantaneous (osculating) argument of pericentre for the two
body problem defined by the relative position `(x,y,z)` and velocity `(u,v,w)` vectors between two bodies
with masses `m1` and `m2`."""
function argperi(
    x, y, z,
    u, v, w,
    m1, m2
    )

    #n = (z-axis unit vector)×h = (-hy, hx, 0)
    n = [x*w-z*u, y*w-z*v, zero(x)]
    e = lrlvec(x,y,z,u,v,w,m1,m2) #cartesian comps. of Laplace-Runge-Lenz vector
    n = n/sqrt(n[1]^2+n[2]^2+n[3]^2)
    e = e/sqrt(e[1]^2+e[2]^2+e[3]^2)
    cosω = dot(n, e)

    if constant_term(e[3]) >= zero(constant_term(x))
        return acos(cosω)
    else
        return 2pi-acos(cosω)
    end

end

"""`longperi(x, y, z, u, v, w, m1, m2)`

Calculates the instantaneous (osculating) longitude of pericentre for the two
body problem defined by the relative position `(x,y,z)` and velocity `(u,v,w)` vectors between two bodies
with masses `m1` and `m2`."""
function longperi(
    x, y, z,
    u, v, w,
    m1, m2
    )

    return argperi(x, y, z, u, v, w, m1, m2)+longascnode(x, y, z, u, v, w)
end

"""`trueanomaly(x, y, z, u, v, w, m1, m2)`

Calculates the instantaneous (osculating) true anomaly for the two
body problem defined by the relative position `(x,y,z)` and velocity `(u,v,w)` vectors between two bodies
with masses `m1` and `m2`."""
function trueanomaly(
    x, y, z,
    u, v, w,
    m1, m2
    )

    R = sqrt(x^2+y^2+z^2)
    V2 = u^2+v^2+w^2
    V = sqrt(V2)

    a = 1.0/( (2.0/R)-(V2/(m1+m2)) )

    h2 = (y*w-z*v)^2 + (z*u-x*w)^2 + (x*v-y*u)^2
    h = sqrt(h2)

    quotient=h2/( (m1+m2)*a )

    e = sqrt(1.0 - quotient)

    Rdot = (x*u+y*v+z*w)/R

    sin_f = Rdot*a*(1.0-e^2)/( e*h )
    cos_f = (1.0/e)*( -1.0+a*(1.0-e^2)/R )

    res = atan(sin_f, cos_f)

    if res ≥ zero(res)
        return res
    else
        return res+2pi
    end
end

"""`ecanomaly(x, y, z, u, v, w, m1, m2)`

Calculates the instantaneous (osculating) eccentric anomaly for the two
body problem defined by the relative position `(x,y,z)` and velocity `(u,v,w)` vectors between two bodies
with masses `m1` and `m2`."""
function ecanomaly(
    x, y, z,
    u, v, w,
    m1, m2
    )

    a = semimajoraxis(x, y, z, u, v, w, m1, m2)
    e = eccentricity(x, y, z, u, v, w, m1, m2)
    R = sqrt(x^2+y^2+z^2)
    Rdot = (x*u+y*v+z*w)/R
    mu = m1+m2

    cosE = (a-R)/(a*e)
    sinE = Rdot*R/( e*sqrt(mu*a) )

    res = atan(sinE, cosE)

    return res

end

"""`meananomaly(x, y, z, u, v, w, m1, m2)`

Calculates the instantaneous (osculating) mean anomaly for the two
body problem defined by the relative position `(x,y,z)` and velocity `(u,v,w)` vectors between two bodies
with masses `m1` and `m2`."""
function meananomaly(
    x, y, z,
    u, v, w,
    m1, m2
    )

    E = ecanomaly(x, y, z, u, v, w, m1, m2)
    e = eccentricity(x, y, z, u, v, w, m1, m2)

    return E-e*sin(E)
end

"""`timeperipass(t, x, y, z, u, v, w, m1, m2)`

Calculates the instantaneous (osculating) time of pericentre passage, at time `t`, for the two
body problem defined by the relative position `(x,y,z)` and velocity `(u,v,w)` vectors between two bodies
with masses `m1` and `m2`."""
function timeperipass(
    t, x, y, z,
    u, v, w,
    m1, m2
    )

    me_an = meananomaly(x, y, z, u, v, w, m1, m2) #mean anomaly
    a = semimajoraxis(x, y, z, u, v, w, m1, m2)
    me_mo = meanmotion(m1+m2, a) #mean motion

    return t-me_an/me_mo
end

"""`rungelenzx(x, y, z, u, v, w, m1, m2)`

Calculates the instantaneous (osculating) x-component of the Laplace-Runge-Lenz vector for the two
body problem defined by the relative position `(x,y,z)` and velocity `(u,v,w)` vectors between two bodies
with masses `m1` and `m2`."""
function rungelenzx(
    x, y, z,
    u, v, w,
    m1, m2
    )

    # v = [u, v, w]
    # h = [y.*w-z.*v, z.*u-x.*w, x.*v-y.*u]
    mu = m1+m2
    # r = [x, y, z]
    rmag = sqrt(x^2+y^2+z^2)

    lrl_x = ( -(z*u-x*w)*w+(x*v-y*u)*v )/(mu)-x/rmag

    return lrl_x
end

"""`rungelenzy(x, y, z, u, v, w, m1, m2)`

Calculates the instantaneous (osculating) y-component of the Laplace-Runge-Lenz vector for the two
body problem defined by the relative position `(x,y,z)` and velocity `(u,v,w)` vectors between two bodies
with masses `m1` and `m2`."""
function rungelenzy(
    x, y, z,
    u, v, w,
    m1, m2
    )

    # v = [u, v, w]
    # h = [y.*w-z.*v, z.*u-x.*w, x.*v-y.*u]
    mu = m1+m2
    # r = [x, y, z]
    rmag = sqrt(x^2+y^2+z^2)

    lrl_y = ( -(x*v-y*u)*u+(y*w-z*v)*w )/(mu)-y/rmag

    return lrl_y
end

"""`rungelenzz(x, y, z, u, v, w, m1, m2)`

Calculates the instantaneous (osculating) z-component of the Laplace-Runge-Lenz vector for the two
body problem defined by the relative position `(x,y,z)` and velocity `(u,v,w)` vectors between two bodies
with masses `m1` and `m2`."""
function rungelenzz(
    x, y, z,
    u, v, w,
    m1, m2
    )

    # v = [u, v, w]
    # h = [y.*w-z.*v, z.*u-x.*w, x.*v-y.*u]
    mu = m1+m2
    # r = [x, y, z]
    rmag = sqrt(x^2+y^2+z^2)

    lrl_z = ( -(y*w-z*v)*v+(z*u-x*w)*u )/(mu)-z/rmag

    return lrl_z
end

"""`rungelenzmag(x, y, z, u, v, w, m1, m2)`

Calculates the instantaneous (osculating) magnitude of the Laplace-Runge-Lenz vector for the two
body problem defined by the relative position `(x,y,z)` and velocity `(u,v,w)` vectors between two bodies
with masses `m1` and `m2`."""
function rungelenzmag(
    x, y, z,
    u, v, w,
    m1, m2
    )

    # v = [u, v, w]
    # h = [y.*w-z.*v, z.*u-x.*w, x.*v-y.*u]
    mu = m1+m2
    # r = [x, y, z]
    rmag = sqrt(x^2+y^2+z^2)

    lrl_x = ( -(z*u-x*w)*w+(x*v-y*u)*v )/(mu)-x/rmag
    lrl_y = ( -(x*v-y*u)*u+(y*w-z*v)*w )/(mu)-y/rmag
    lrl_z = ( -(y*w-z*v)*v+(z*u-x*w)*u )/(mu)-z/rmag

    return sqrt(lrl_x^2+lrl_y^2+lrl_z^2)
end

"""`lrlvec(x, y, z, u, v, w, m1, m2)`

Calculates the instantaneous (osculating) cartesian components of the Laplace-Runge-Lenz vector for the two
body problem defined by the relative position `(x,y,z)` and velocity `(u,v,w)` vectors between two bodies
with masses `m1` and `m2`."""
function lrlvec(
    x, y, z,
    u, v, w,
    m1, m2
    )

    mu = m1+m2
    rmag = sqrt(x^2+y^2+z^2)

    lrl_x = ( -(z*u-x*w)*w+(x*v-y*u)*v )/(mu)-x/rmag
    lrl_y = ( -(x*v-y*u)*u+(y*w-z*v)*w )/(mu)-y/rmag
    lrl_z = ( -(y*w-z*v)*v+(z*u-x*w)*u )/(mu)-z/rmag

    return [lrl_x, lrl_y, lrl_z]
end

# compute eccentric anomaly (E) from eccentricity (e) and mean anomaly (M)
function eccentricanomaly(e::T, M::T) where {T <: Number}
    M0 = mod2pi(M)
    E0 = M0 + sign(sin(M0))*0.85*e #Murray-Dermotts' initial estimate
    # successive approximations via Newtons' method
    for i in 0:4
        #TODO: implement modified Newton's method for Kepler's equation (Murray-Dermott)
        Eans = E0 - (E0-e*sin(E0)-M0)/(1.0-e*cos(E0))
        E0 = Eans
    end
    return E0
end

# compute true anomaly (f) from eccentricity (e) and eccentric anomaly (E)
function trueanomaly(e,E)
    Enew = mod2pi(E)
    return 2.0*atan(  sqrt((1.0+e)/(1.0-e))*tan(Enew/2)  )
end

# compute true anomaly from eccentricity and mean anomaly
function meanan2truean(e,M)
    return trueanomaly(e, eccentricanomaly(e, M))
end

# get mean motion from mass parameter (mu) and semimajor axis (a)
function meanmotion(mu,a)
    return sqrt(mu/(a^3))
end

# get mean anomaly from mean motion (n), time (t) and time of pericenter passage (taup)
function meananomaly(n, t, taup)
    return n*(t-taup)
end

# compute true anomaly from time, a, e, mu and taup
function time2truean(a, e, mu, t, taup)
    return meanan2truean(e, meananomaly(meanmotion(mu, a), t, taup))
end
