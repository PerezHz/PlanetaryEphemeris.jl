using PrecompileTools

@setup_workload begin
    # Maximum number of steps
    maxsteps = 2
    # Starting time of integration
    jd0 = datetime2julian(DateTime(2000,1,1,12)) 
    # Number of years 
    nyears = 2031.0 - year(julian2datetime(jd0))
    # Dense output
    dense = Val(true)
    # Dynamical function 
    dynamics = DE430!
    @compile_workload begin
        # All calls in this block will be precompiled, regardless of whether
        # they belong to your package or not (on Julia 1.8 and higher)
        propagate(maxsteps, jd0, nyears, dense; dynamics = dynamics, order = order, abstol = abstol)
    end
end