using PrecompileTools

@setup_workload begin
    # Initial conditions and their epoch [JD TDB]
    const filename = joinpath(pkgdir(PlanetaryEphemeris), "data", "de430ic.txt")
    const jd0 = 2440400.5
    # Time span
    const tspan = (jd0, jd0 + yr)
    const q0 = read_initial_conditions(filename)
    # Parameters
    const maxsteps = 1
    const precompile_order = 15
    const abstol = 1E-12
    const parse_eqs = true
    const params = DE430Params(jd0, q0, precompile_order; parse_eqs)
    # Planetary ephemeris problem
    const PE = PlanetaryEphemerisProblem(DE430!, tspan, q0, params)
    @compile_workload begin
        # All calls in this block will be precompiled, regardless of whether
        # they belong to your package or not (on Julia 1.8 and higher)
        propagate(PE; maxsteps, order = precompile_order, abstol, parse_eqs)
    end
end
