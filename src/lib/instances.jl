@doc raw"""
    instance(:unit, p::Integer, ::Type{K})::Matrix{K}
    instance(:cube, N::Integer, ::Type{K})::Matrix{K}
""" function instance end

function instance(code::Symbol, args...)
    return instance(Val(code), args...)
end

# Unit Circle on p points
function instance(::Val{:unit}, p::Integer, ::Type{K}=Float64) where {K}
    points = Matrix{K}(undef, 2, p)

    for j = 1:p
        z = im * exp(2π * (j - 1) * im / p)

        points[:, j] .= [real(z), imag(z)]
    end

    return points
end

# N-dimensional unit cube
function instance(::Val{:cube}, N::Integer, ::Type{K}=Float64) where {K}
    points = Matrix{K}(undef, N, 2^N)

    for j = 1:2^N
        points[:, j] .= digits(j - 1; base = 2, pad = N)
    end

    return points
end

# function SteinerTree{N,K}(code::Symbol, args...; kw...) where {N,K}
#     return SteinerTree{N,K}(Val(code), args...; kw...)
# end

# function SteinerTree{2,K}(::Val{:unit}, p::Integer) where {K}
#     X = Vector{Vector{K}}(undef, p)

#     for j = 1:p
#         z = im * exp(2π * (j - 1) * im / p)

#         X[j] = [real(z), imag(z)]
#     end

#     return SteinerTree{2,K}(X)
# end