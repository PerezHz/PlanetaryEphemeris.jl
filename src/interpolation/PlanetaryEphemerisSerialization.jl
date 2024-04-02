# Custom serialization

@doc raw"""
    PlanetaryEphemerisSerialization{T}

Custom serialization struct to save a `TaylorInterpolant{T, T, 2}` to a `.jld2` file.

# Fields
- `order::Int`: order of Taylor polynomials.
- `dims::Tuple{Int, Int}`: matrix dimensions.
- `t0::T`: initial time.
- `t::Vector{T}`: vector of times.
- `x::Vector{T}`: vector of coefficients.
"""
struct PlanetaryEphemerisSerialization{T}
    order::Int
    dims::Tuple{Int, Int}
    t0::T
    t::Vector{T}
    x::Vector{T}
end

# Tell JLD2 to save TaylorInterpolant{T, T, 2} as PlanetaryEphemerisSerialization{T}
function writeas(::Type{<:TaylorInterpolant{T, T, 2, Vector{T}, Matrix{Taylor1{T}}}}) where {T<:Real}
    return PlanetaryEphemerisSerialization{T}
end

# Convert method to write .jld2 files
function convert(::Type{PlanetaryEphemerisSerialization{T}}, eph::TaylorInterpolant{T, T, 2, Vector{T}, Matrix{Taylor1{T}}}) where {T <: Real}
    # Taylor polynomials order
    order = eph.x[1, 1].order
    # Number of coefficients in each polynomial
    k = order + 1
    # Matrix dimensions
    dims = size(eph.x)
    # Number of elements in matrix
    N = dims[1] * dims[2]
    # Vector of coefficients
    x = Vector{T}(undef, k * N)
    # Save coefficients
    for i in 1:N
        x[(i-1)*k+1 : i*k] = eph.x[i].coeffs
    end

    return PlanetaryEphemerisSerialization{T}(order, dims, eph.t0, eph.t, x)
end

# Convert method to read .jld2 files
function convert(::Type{TaylorInterpolant{T, T, 2, Vector{T}, Matrix{Taylor1{T}}}}, eph::PlanetaryEphemerisSerialization{T}) where {T<:Real}
    # Taylor polynomials order
    order = eph.order
    # Number of coefficients in each polynomial
    k = order + 1
    # Matrix dimensions
    dims = eph.dims
    # Number of elements in matrix
    N = dims[1] * dims[2]
    # Matrix of Taylor polynomials
    x = Matrix{Taylor1{T}}(undef, dims[1], dims[2])
    # Reconstruct Taylor polynomials
    for i in 1:N
        x[i] = Taylor1{T}(eph.x[(i-1)*k+1 : i*k], order)
    end

    return TaylorInterpolant{T, T, 2}(eph.t0, eph.t, x)
end

function convert(::Type{TaylorInterpolant{T, T, 2}}, eph::PlanetaryEphemerisSerialization{T}) where {T<:Real}
    return convert(TaylorInterpolant{T, T, 2, Vector{T}, Matrix{Taylor1{T}}}, eph)
end