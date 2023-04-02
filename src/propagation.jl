@doc raw"""
    Taylor1Serialization{T}

Specialized struct to save `Taylor1{T}` objects to `.jld2` files.
"""
struct Taylor1Serialization{T}
    x::Vector{T}
end

# Tell JLD2 to save Taylor1{T} as Taylor1Serialization{T}
writeas(::Type{Taylor1{T}}) where {T} = Taylor1Serialization{T}
# Convert method to write .jld2 files
convert(::Type{Taylor1Serialization{T}}, a::Taylor1{T}) where {T} = Taylor1Serialization{T}(a.coeffs)
# Convert method to read .jld2 files
convert(::Type{Taylor1{T}}, a::Taylor1Serialization{T}) where {T} = Taylor1{T}(a.x, length(a.x) - 1)

@doc raw"""
    selecteph2jld2(sseph::TaylorInterpolant, bodyind::T, tspan::S, N::Int) where {T <: AbstractVector{Int}, S <: Number}

Save the ephemeris, contained in `sseph`, of the bodies with indices `bodyind`, in a `.jld2` file named as follows

    "sseph" * number of asteroids in sseph * "ast" * number of asteroids to be saved in file * "_"
    * "p" / "m" (forward / backward integration) * number of years in sseph *  "y_et.jld2"

# Arguments

- `sseph::TaylorInterpolant`: ephemeris of all the bodies.
- `bodyind::T`: indices of the bodies to be saved.
- `tspan::S`: time span of the integration (positive -> forward integration / negative -> backward integration).
- `N::Int`: total number of bodies.
"""
function selecteph2jld2(sseph::TaylorInterpolant, bodyind::T, tspan::S, N::Int) where {T <: AbstractVector{Int}, S <: Number}

    # Number of asteroids in sseph
    nast = N - 11
    # indices of the positions and velocities of the bodies to be saved
    indvec = nbodyind(N, bodyind)
    # Number of asteroids to be saved
    nastout = length(bodyind) - 11
    # Check nastout <= nast
    @assert nastout <= nast "Cannot save $nastout asteroids from ephemeris with $nast asteroids"
    # Prefix to distinguish between forward (p) / backward (m) integration
    sgn_yrs = signbit(tspan) ? "m" : "p"
    # Number of years
    nyrs_int = Int(abs(tspan))

    # Write output to .jld2 file

    # Name of the file
    ss16ast_fname = "sseph$(lpad(nast,3,'0'))ast$(lpad(nastout,3,'0'))_" * sgn_yrs * "$(nyrs_int)y_et.jld2"

    # TaylorInterpolant with only the information of the bodies to be saved + Lunar orientation + TT-TDB
    ss16ast_eph = TaylorInterpolant(sseph.t0, sseph.t, sseph.x[:, union(indvec, 6N+1:6N+13)])

    println("Saving solution to file: ", ss16ast_fname)

    # Open file
    JLD2.jldopen(ss16ast_fname, "w") do file
        # Write the ephemeris to file
        write(file, "ss16ast_eph", ss16ast_eph)
    end
    # Check that written output is equal to original variable ss16ast_eph
    recovered_sol_i = JLD2.load(ss16ast_fname, "ss16ast_eph")
    if recovered_sol_i == ss16ast_eph
        println("Solution saved correctly")
    else
        println("Saved and recovered solution are not equal")
    end

    return ss16ast_fname
end

@doc raw"""
    save2jld2andcheck(outfilename::String, sol)

Save `sol` in `outfilename` (.jld2) and check that recovered solution equals `sol`.
"""
function save2jld2andcheck(outfilename::String, sol)

    println("Saving solution to file: ", outfilename)

    # Open file
    JLD2.jldopen(outfilename, "w") do file
        # Loop over solution variables
        for ind in eachindex(sol)
            # Name of the variable
            varname = string(ind)
            println("Saving variable: ", varname)
            # Write the variable
            write(file, varname, sol[ind])
        end
    end

    # Check that saved solution is equal to the original
    println("Checking that all variables were saved correctly...")

    # Loop over solution variables
    for ind in eachindex(sol)
        # Name of the variable
        varname = string(ind)
        # Read varname from files and assign recovered variable to recovered_sol_i
        recovered_sol_i = JLD2.load(outfilename, varname)
        # Check that varname was recovered succesfully
        if recovered_sol_i == sol[ind]
            println("Variable ", varname, " saved correctly" )
        else
            println("Recovered variable ", varname, " is not equal to the original" )
        end
    end

    println("Saved solution")

    return nothing
end

@doc raw"""
    day2sec(x::Matrix{Taylor1{U}}) where {U <: Number}

Convert `x` from days to seconds.
"""
function day2sec(x::Matrix{Taylor1{U}}) where {U <: Number}

    # Order of Taylor polynomials
    order = x[1, 1].order
    # Matrix dimensions
    m, n = size(x)
    # Taylor conversion variable
    t = Taylor1(order) / daysec
    # Allocate memory
    res = Matrix{Taylor1{U}}(undef, m, n)
    # Iterate over the matrix
    for j in 1:n
        for i in 1:m
            @inbounds res[i, j] = x[i, j](t)
        end
    end

    return res
end

@doc raw"""
    propagate(maxsteps::Int, jd0::T, tspan::T, ::Val{false/true}; dynamics::Function = NBP_pN_A_J23E_J23M_J2S!,
              nast::Int = 343, order::Int = order, abstol::T = abstol, parse_eqs::Bool = true) where {T <: Real}

Integrate the Solar System via the Taylor method.

# Arguments

- `maxsteps::Int`: maximum number of steps for the integration.
- `jd0::T`: initial Julian date.
- `tspan::T`: time span of the integration (in Julian days).
- `::Val{false/true}`: whether to save the Taylor polynomials at each step (`true`) or not (`false`).
- `output::Bool`: whether to write the output to a file (`true`) or not.
- `ephfile::String`: name of the file where to save the solution if `ss16ast` is `false`.
- `dynamics::Function`: dynamical model function.
- `nast::Int`: number of asteroids to be considered in the integration.
- `ss16ast::Bool`: whether to save the solution using `selecteph2jld2` (`true`) or not.
- `bodyind::AbstractVector{Int}`: indices of the bodies to be saved.
- `order::Int=order`: order of the Taylor expansions to be used in the integration.
- `abstol::T`: absolute tolerance.
- `parse_eqs::Bool`: whether to use the specialized method of `jetcoeffs!` (`true`) created with `@taylorize` or not.
""" propagate

for V_dense in (:(Val{true}), :(Val{false}))
    @eval begin

        function propagate(maxsteps::Int, jd0::T, tspan::T, ::$V_dense; dynamics::Function = NBP_pN_A_J23E_J23M_J2S!,
                           nast::Int = 343, order::Int = order, abstol::T = abstol, parse_eqs::Bool = true) where {T <: Real}

            # Total number of bodies (Sun + 8 planets + Moon + Pluto + Asteroid)
            N = 11 + nast

            # Get 6N + 13 initial conditions (3N positions + 3N velocities + 6 lunar mantle angles + 6 lunar core angles + TT-TDB)
            q0 = initialcond(N, jd0)

            # Set initial time equal to zero (improves accuracy in data reductions)
            t0 = zero(T)

            # Parameters for dynamical function
            params = (N, jd0)


            # Final time of integration (days)
            tmax = t0 + tspan*yr

            # Integration
            sol_ = @time taylorinteg_threads(dynamics, q0, t0, tmax, order, abstol, $V_dense(), params, maxsteps = maxsteps,
                                             parse_eqs = parse_eqs)

            if $V_dense == Val{true}

                # Parameters for TaylorInterpolant

                # Initial time [ days -> seconds ]
                et0 = ( (jd0 - J2000) * daysec ) :: T

                # Vector of times [ days -> seconds ]
                etv = (sol_.t * daysec) :: Vector{T}

                # Vector of Taylor polynomials [ days -> seconds ]
                sseph_x_et = day2sec(sol_.x) :: Matrix{Taylor1{T}}

                # Save ephemeris in TaylorInterpolant object
                sseph = TaylorInterpolant{T, T, 2}(et0, etv, sseph_x_et)

                return sseph

            else

                sol = (t = sol_[1][:], x = sol_[2][:, :])

                return sol

            end

        end

        function propagate(maxsteps::Int, jd0::T1, tspan::T2, ::$V_dense; dynamics::Function = NBP_pN_A_J23E_J23M_J2S!,
                           nast::Int = 343, order::Int = order, abstol::T3 = abstol, parse_eqs::Bool = true) where {T1, T2, T3 <: Real}

            _jd0, _tspan, _abstol = promote(jd0, tspan, abstol)

            return propagate(maxsteps, _jd0, _tspan, $V_dense(); dynamics = dynamics, nast = nast, order = order,
                             abstol = abstol, parse_eqs = parse_eqs)

        end

    end
end
