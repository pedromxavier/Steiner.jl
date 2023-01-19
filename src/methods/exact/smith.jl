@doc raw"""
    SmithMethod

Smith (1992) enumeration method
""" struct SmithMethod <: ExactMethod end

function solve(::SmithMethod, points::Matrix{K}) where {K}
    N = size(points, 2) 
    t = SteinerTree{N,K}(points)

    # 1 Pick Triangle

    # 1.1 Pick random vertex `i`
    i = rand(t.terminal_points)

    # 1.2 Pick the farthest vertex `j` from `i`
    j = argmax(
        j -> t.matrix[i,j],
        filter(j -> j != i, t.terminal_points)
    )

    # 1.3 Pick vertex `k`, whose minimum distance between `i` and `j` is the greatest
    k = argmax(
        k -> min(t.matrix[i,k], t.matrix[j,k]),
        filter(k -> k != i && k != j, t.terminal_points),
    )

    # 2 Solve Fermat's problem for △ = (i, j, k).

    # Embed points in ℝⁿ
    u = t.points[i]
    v = t.points[j]
    w = t.points[k]

    



end