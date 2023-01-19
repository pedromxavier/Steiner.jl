@doc raw"""
    DelaunayTriangulation

Delaunay triangulation method
""" struct DelaunayTriangulation <: ExactMethod end

function solve(::DelaunayTriangulation, points::Matrix{K}) where {K}
    return nothing
end