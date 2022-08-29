"""
    vector_field_colors(Xs, Ys, f::Function)
    
Map the magnitude of a vector to a color.
I don't know why this works. The answer is taken from
https://discourse.julialang.org/t/quiver-plot-plots-pyplot-with-different-colors-depending-on-the-length-of-the-arrow/59577/5
"""
function vector_field_colors(Xs, Ys, f::Function)
    c = norm.(f.(Xs, Ys))
    c = [c c]'
    c = repeat([c...], inner=2)

    return c
end