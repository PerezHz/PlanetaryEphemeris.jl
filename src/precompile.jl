using PrecompileTools

@setup_workload begin
    # Time span
    const tspan = (J2000, J2000 + yr)
    # Initial conditions
    const filename = joinpath(pkgdir(PlanetaryEphemeris), "data", "de430ic_2000Jan1.txt")
    const q0 = read_initial_conditions(filename)
    # Parameters
    const maxsteps = 1
    const order = 15
    const abstol = 1E-12
    const parse_eqs = true
    const params = DE430Params(J2000, q0, order; parse_eqs)
    # Planetary ephemeris problem
    const PE = PlanetaryEphemerisProblem(DE430!, tspan, q0, params)
    @compile_workload begin
        # All calls in this block will be precompiled, regardless of whether
        # they belong to your package or not (on Julia 1.8 and higher)
        propagate(PE; maxsteps, order, abstol, parse_eqs)
    end
end
