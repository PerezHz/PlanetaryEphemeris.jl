using PrecompileTools

@setup_workload begin
    # Initial conditions
    const filename = joinpath(pkgdir(PlanetaryEphemeris), "data", "de430ic_2000Jan1.txt")
    const initcond = read_initial_conditions(filename)
    # Parameters
    const params = (354, J2000)
    # Planetary ephemeris problem
    const PE = PlanetaryEphemerisProblem(DE430!, J2000, initcond, params)
    # Number of years
    const nyears = 1.0
    @compile_workload begin
        # All calls in this block will be precompiled, regardless of whether
        # they belong to your package or not (on Julia 1.8 and higher)
        propagate(PE, nyears; maxsteps = 1, order = 15, abstol = 1E-12, parse_eqs = true)
    end
end
