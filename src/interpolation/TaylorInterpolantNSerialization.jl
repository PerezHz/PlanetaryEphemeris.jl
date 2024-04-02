@doc raw"""
    TaylorInterpolantNSerialization{T}

Custom serialization struct to save a `TaylorInterpolant{T, TaylorN{T}, 2}` to a `.jld2` file.

# Fields
- `vars::Vector{String}`: jet transport variables.
- `order::Int`: order of Taylor polynomials w.r.t time.
- `varorder::Int`: order of jet transport perturbations.
- `dims::Tuple{Int, Int}`: matrix dimensions.
- `t0::T`: initial time.
- `t::Vector{T}`: vector of times.
- `x::Vector{T}`: vector of coefficients.
"""
struct TaylorInterpolantNSerialization{T}
    vars::Vector{String}
    order::Int
    varorder::Int
    dims::Tuple{Int, Int}
    t0::T
    t::Vector{T}
    x::Vector{T}
end

# Tell JLD2 to save <:TaylorInterpolant{T, TaylorN{T}, 2} as TaylorInterpolantNSerialization{T}
function writeas(::Type{<:TaylorInterpolant{T, TaylorN{T}, 2, Vector{T}, Matrix{Taylor1{TaylorN{T}}}}}) where {T<:Real}
    return TaylorInterpolantNSerialization{T}
end

# Convert method to write .jld2 files
function convert(::Type{TaylorInterpolantNSerialization{T}}, eph::TaylorInterpolant{T, TaylorN{T}, 2, Vector{T}, Matrix{Taylor1{TaylorN{T}}}}) where {T}
    # Variables
    vars = TS.get_variable_names()
    # Number of variables
    n = length(vars)
    # Matrix dimensions
    dims = size(eph.x)
    # Number of elements in matrix
    N = length(eph.x)
    # Taylor1 order
    order = eph.x[1, 1].order
    # Number of coefficients in each Taylor1
    k = order + 1
    # TaylorN order
    varorder = eph.x[1, 1].coeffs[1].order
    # Number of coefficients in each TaylorN
    L = varorder + 1
    # Number of coefficients in each HomogeneousPolynomial
    M = binomial(n + varorder, varorder)

    # Vector of coefficients
    x = Vector{T}(undef, N * k * M)

    # Save coefficients
    i = 1
    # Iterate over matrix elements
    for i_1 in 1:N
        # Iterate over Taylor1 coefficients
        for i_2 in 1:k
            # Iterate over TaylorN coefficients
            for i_3 in 0:varorder
                # Iterate over i_3 order HomogeneousPolynomial
                for i_4 in 1:binomial(n + i_3 - 1, i_3)
                    x[i] = eph.x[i_1].coeffs[i_2].coeffs[i_3+1].coeffs[i_4]
                    i += 1
                end
            end
        end
    end

    return TaylorInterpolantNSerialization{T}(vars, order, varorder, dims, eph.t0, eph.t, x)
end

# Convert method to read .jld2 files
function convert(::Type{TaylorInterpolant{T, TaylorN{T}, 2, Vector{T}, Matrix{Taylor1{TaylorN{T}}}}}, eph::TaylorInterpolantNSerialization{T}) where {T}
    # Variables
    vars = eph.vars
    # Number of variables
    n = length(vars)
    # Matrix dimensions
    dims = eph.dims
    # Number of elements in matrix
    N = dims[1] * dims[2]
    # Taylor1 order
    order = eph.order
    # Number of coefficients in each Taylor1
    k = order + 1
    # TaylorN order
    varorder = eph.varorder
    # Number of coefficients in each TaylorN
    L = varorder + 1
    # Number of coefficients in each HomogeneousPolynomial
    M = binomial(n + varorder, varorder)
    # M = sum(binomial(n + i_3 - 1, i_3) for i_3 in 0:varorder)

    # Set variables
    if TS.get_variable_names() != vars
        TS.set_variables(T, vars, order = varorder)
    end

    # Matrix of Taylor polynomials
    x = Matrix{Taylor1{TaylorN{T}}}(undef, dims[1], dims[2])

    # Reconstruct Taylor polynomials
    i = 1
    # Iterate over matrix elements
    for i_1 in 1:N
        # Reconstruct Taylor1s
        Taylor1_coeffs = Vector{TaylorN{T}}(undef, k)
        for i_2 in 1:k
            # Reconstruct TaylorNs
            TaylorN_coeffs = Vector{HomogeneousPolynomial{T}}(undef, L)
            # Reconstruct HomogeneousPolynomials
            for i_3 in 0:varorder
                TaylorN_coeffs[i_3 + 1] = HomogeneousPolynomial(eph.x[i : i + binomial(n + i_3 - 1, i_3)-1], i_3)
                i += binomial(n + i_3 - 1, i_3)
            end
            Taylor1_coeffs[i_2] = TaylorN(TaylorN_coeffs, varorder)
        end
        x[i_1] = Taylor1{TaylorN{T}}(Taylor1_coeffs, order)
    end

    return TaylorInterpolant{T, TaylorN{T}, 2}(eph.t0, eph.t, x)
end

function convert(::Type{TaylorInterpolant{T, TaylorN{T}, 2}}, eph::TaylorInterpolantSerialization{T}) where {T<:Real}
    return convert(TaylorInterpolant{T, TaylorN{T}, 2, Vector{T}, Matrix{Taylor1{TaylorN{T}}}}, eph)
end