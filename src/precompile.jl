using PrecompileTools

@setup_workload begin
    # Starting time of integration
    jd0 = datetime2julian(DateTime(2000,1,1,12)) 
    # Number of years 
    nyears = 2031.0 - year(julian2datetime(jd0))
    @compile_workload begin
        # All calls in this block will be precompiled, regardless of whether
        # they belong to your package or not (on Julia 1.8 and higher)
        propagate(1, jd0, nyears, Val(true); dynamics = DE430!, order = 25, abstol = 1.0E-20)
    end
end