using SnoopPrecompile

@info "SnoopPrecompile is analyzing PlanetaryEphemeris.jl code..."

@precompile_setup begin
    # Starting time of integration
    jd0 = datetime2julian(DateTime(2000,1,1,12)) 
    # Number of years 
    nyears = 2031.0 - year(julian2datetime(jd0))

    @precompile_all_calls begin
        propagate(1, jd0, nyears, Val(false); dynamics = DE430!, order = 25, abstol = 1.0E-20)
        propagate(1, jd0, nyears, Val(true); dynamics = DE430!, order = 25, abstol = 1.0E-20)
    end 
end 