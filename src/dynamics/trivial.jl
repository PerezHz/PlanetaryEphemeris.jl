# This file is part of the PlanetaryEphemeris.jl package; MIT licensed

@taylorize function trivialdynamics!(dq, q, params, t)
    Threads.@threads for j in eachindex(dq)
        dq[j] = zero(q[j])
    end
    nothing
end

@taylorize function freeparticle!(dq, q, params, t)
    # N: number of bodies
    # jd0: initial Julian date
    local N, jd0 = params
    Threads.@threads for j in 1:N
        # Fill first 3N elements of dq with velocities
        dq[3j-2] = q[3(N+j)-2]
        dq[3j-1] = q[3(N+j)-1]
        dq[3j  ] = q[3(N+j)  ]
    end
    Threads.@threads for i in 1:N
        dq[3(N+i)-2] = zero(q[3(N+i)-2])
        dq[3(N+i)-1] = zero(q[3(N+i)-1])
        dq[3(N+i)  ] = zero(q[3(N+i)  ])
    end
    nothing
end
