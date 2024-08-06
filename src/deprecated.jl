Base.@deprecate propagate(maxsteps, jd0, tspan, ::Val{true}; dynamics = NBP_pN_A_J23E_J23M_J2S_threads!,
                   nast = 343, order = order, abstol = abstol, parse_eqs = true) propagate(maxsteps, jd0, tspan; dynamics = NBP_pN_A_J23E_J23M_J2S_threads!,
                   nast = 343, order = order, abstol = abstol, parse_eqs = true)
Base.@deprecate propagate(maxsteps, jd0, tspan, ::Val{false}; dynamics = NBP_pN_A_J23E_J23M_J2S_threads!,
                   nast = 343, order = order, abstol = abstol, parse_eqs = true) error("PlanetaryEphemeris.propagate may only be called with `Val(true)` in the fourth argument.")
