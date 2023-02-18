@doc raw"""
    NonLinear

NonLinear enumeration method
""" struct NonLinear <: ExactMethod end

function build(::NonLinear, points::Matrix{K}; optimizer) where {K}
    N = size(points, 1) # space dimension
    p = size(points, 2) # terminal points
    s = p - 2

    model = Model(optimizer)

    # Coordinates
    @variable(model, x[1:N,1:s])

end

function solve(method::NonLinear, points::Matrix{K}; optimizer) where {K}
    model = build(method, points; optimizer)

    optimize!(model)

    return model
end