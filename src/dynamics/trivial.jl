@taylorize function trivialdynamics!(dq, q, params, t)
    Threads.@threads for j in eachindex(dq)
        dq[j] = zero(q[j])
    end
    nothing
end
