struct SteinerTree{N,K<:AbstractFloat}
    points::Vector{Vector{K}}
    # matrix::DistanceMatrix{K}
    graph::SimpleGraph{Int}
    steiner_points::BitSet
    terminal_points::BitSet

    # Problem Bounds
    lower_bound::K
    upper_bound::K

    function SteinerTree{N,K}(points::Vector{Vector{K}}) where {N,K}
        # Assert that all points have dimension N
        @assert all(length.(points) .== N)

        p = length(points)
        s = p - 2
        n = p + s

        steiner_tree    = SimpleGraph{Int}(n)
        steiner_points  = BitSet()
        terminal_points = BitSet(1:p)

        return new{N,K}(
            [points; Vector{Vector{K}}(undef, s)],
            steiner_tree,
            spanning_tree,
            steiner_points,
            terminal_points,
        )   
    end
end

function SteinerTree{N}(args...; kw...) where {N}
    return SteinerTree{N,Float64}(args...; kw...)
end