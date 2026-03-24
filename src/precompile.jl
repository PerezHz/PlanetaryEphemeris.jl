using PrecompileTools

@setup_workload begin
    # Time span
    const tspan = (J2000, J2000 + yr)
    # Initial conditions
    const filename = joinpath(pkgdir(PlanetaryEphemeris), "data", "de430ic_2000Jan1.txt")
    const q0 = read_initial_conditions(filename)
    # Number of bodies
    const N = (length(q0) - 13) ÷ 6
    # Parameters
    const params = (N, J2000)
    # Planetary ephemeris problem
    const PE = PlanetaryEphemerisProblem(DE430!, tspan, q0, params)
    @compile_workload begin
        # All calls in this block will be precompiled, regardless of whether
        # they belong to your package or not (on Julia 1.8 and higher)
        propagate(PE; maxsteps = 1, order = 15, abstol = 1E-12, parse_eqs = true)
    end
end
