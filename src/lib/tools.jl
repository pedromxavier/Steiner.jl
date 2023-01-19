const POINTS_VIEW{K} = SubArray{K,2,Matrix{K},Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}}, true}

struct DistanceMatrix{K} <: AbstractMatrix{K}
    points::POINTS_VIEW{K}

    function DistanceMatrix{K}(points::Matrix{K}) where {K}
        return new{K}(view(points, :, :))
    end
end

function getindex(m::DistanceMatrix{K}, i::Integer, j::Integer) where {K}
    return norm(m.points[:, i] - m.points[:, j])
end

function angle(x::Vector{K}, y::Vector{K}) where {K}
    u = x ./ norm(x)
    v = y ./ norm(y)
    θ = 2 * atan(norm2(u + v), norm2(u - v))
    
    return ifelse(
        signbit(θ) || signbit(K(π) - θ),
        ifelse(
            signbit(θ),
            zero(K),
            K(π)
        ),
        θ
    )
end

@doc raw"""
""" function solve_fermat(u::Vector{K}, v::Vector{K}, w::Vector{K}) where {K}
    θ = K(2π/3)
    α = angle(v - u, w - u); α >= θ && return u
    β = angle(u - v, w - v); β >= θ && return v
    γ = angle(u - w, v - w); γ >= θ && return w

    x = u - v
    y = w - v
    z = normalize!(cross(x, y))
    LinearAlgebra.Rotation()
end