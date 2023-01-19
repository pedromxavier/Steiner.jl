@doc raw"""
""" abstract type AbstractMethod end

@doc raw"""
    solve(::AbstractMethod, points::Matrix{T}) where {T}
    solve(::AbstractMethod, points::Vector{Vector{T}}) where {T}
""" function solve end