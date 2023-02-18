module Steiner

using Graphs:
    AbstractGraph,
    Graphs,
    SimpleGraph,
    add_vertices!,
    edges,
    vertices
using JuMP
using LinearAlgebra
using RecipesBase:
    @recipe,
    @series

# export SteinerTree

# function steiner_tree_length(st::SteinerTree{N,K}) where {N,K}
#     l = zero(K)

#     for e in edges(st.steiner_tree)
#         l += path_length(st, e.src, e.dst)
#     end

#     return l
# end

include("lib/tree.jl")
include("lib/view.jl")
include("lib/tools.jl")
include("lib/instances.jl")

include("interface.jl")

include("methods/exact/exact.jl")
include("methods/heuristic/heuristic.jl")

end # module Steiner
