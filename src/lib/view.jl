function superscript(n::Integer)
    if n < 0
        return "⁻$(superscript(-n))"
    elseif n >= 10
        return "$(superscript(n ÷ 10))$(superscript(n % 10))"
    elseif n == 1
        return "¹"
    elseif n == 2
        return "²"
    else
        return "$('⁰' + n)"
    end
end

function Base.show(io::IO, st::SteinerTree{N,K}) where {N,K}
    print(
        io,
        """
        Steiner Tree ∈ ℝ$(superscript(N))

        ▷ Terminal Points…$(length(st.terminal_points))
        ▷ Steiner Points……$(length(st.steiner_points))
        """
    )

    return nothing
end

@recipe function plot(st::SteinerTree{2,K}) where {K}
    title      --> "Steiner Tree (ℓ = $(steiner_tree_length(st)))"
    xlabel     --> "x"
    ylabel     --> "y"

    for e in edges(st.steiner_tree)
        @series begin
            label   := nothing
            color   := :red
            
            i = e.src
            j = e.dst

            x = [st.points[i][1], st.points[j][1]]
            y = [st.points[i][2], st.points[j][2]]

            (x, y)
        end
    end

    # Steiner Points
    @series begin
        label      --> "Steiner"
        seriestype :=  :scatter
        color      :=  :red
        
        x = Vector{K}(undef, length(st.steiner_points))
        y = Vector{K}(undef, length(st.steiner_points))
        
        for (i, s) in enumerate(st.steiner_points)
            x[i], y[i] = st.points[s]
        end

        (x, y)
    end

    # Terminal Points
    @series begin
        label      --> "Terminal"
        seriestype :=  :scatter
        color      :=  :blue

        x = Vector{K}(undef, length(st.terminal_points))
        y = Vector{K}(undef, length(st.terminal_points))
        
        for (i, t) in enumerate(st.terminal_points)
            x[i], y[i] = st.points[t]
        end

        return (x, y)
    end

    return ()
end