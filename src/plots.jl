"""
    plot_vector_field(xs, ys, Σ::ControlAffineSystem; scale=0.15, lw=1)
    plot_vector_field(xs, ys, Σ::ControlAffineSystem, k::FeedbackController; scale=0.15, lw=1)

Plot the vector field of a open or closed-loop `ControlAffineSystem`.
"""
function plot_vector_field(xs, ys, Σ::ControlAffineSystem; scale=0.15, lw=1)
    f(x1, x2) = Σ.f([x1, x2])

    return VectorFieldPlots.plot_vector_field(xs, ys, f, scale=scale, lw=lw)
end

function plot_vector_field(xs, ys, Σ::ControlAffineSystem, k::FeedbackController; scale=0.15, lw=1)
    f(x1, x2) = Σ.f([x1, x2]) + Σ.g([x1, x2])*k([x1, x2])

    return VectorFieldPlots.plot_vector_field(xs, ys, f, scale=scale, lw=lw)
end

"""
    plot_vector_field!(xs, ys, Σ::ControlAffineSystem; scale=0.15, lw=1)
    plot_vector_field!(xs, ys, Σ::ControlAffineSystem, k::FeedbackController; scale=0.15, lw=1)

Plot the vector field of a open or closed-loop `ControlAffineSystem` on an existing figure.
"""
function plot_vector_field!(xs, ys, Σ::ControlAffineSystem; scale=0.15, lw=1)
    f(x1, x2) = Σ.f([x1, x2])

    return VectorFieldPlots.plot_vector_field!(xs, ys, f, scale=scale, lw=lw)
end

function  plot_vector_field!(xs, ys, Σ::ControlAffineSystem, k::FeedbackController; scale=0.15, lw=1)
    f(x1, x2) = Σ.f([x1, x2]) + Σ.g([x1, x2])*k([x1, x2])

    return VectorFieldPlots.plot_vector_field!(xs, ys, f, scale=scale, lw=lw)
end

"""
    vector_field_colors(Xs, Ys, f::Function)

Map the magnitude of a vector to a color.

To make this work you have to do some weird stuff. I don't know why, and the answer is taken
from https://discourse.julialang.org/t/quiver-plot-plots-pyplot-with-different-colors-depending-on-the-length-of-the-arrow/59577/5
"""
# function vector_field_colors(Xs, Ys, f::Function)
#     c = norm.(f.(Xs, Ys))
#     c = [c c]'
#     c = repeat([c...], inner=2)

#     return c
# end

"""
    plot_circle(x, y, r; samples=500, lw=1.0, c=:black, fillalpha=0.2)
    plot_circle!(x, y, r; samples=500, lw=1.0, c=:gray70, fillalpha=1.0)

Plot a circle with center (x,y) and radius r.
"""
function plot_circle(x, y, r; samples=500, lw=1.0, c=:black, fillalpha=0.2)
    xs, ys = circle_shape(x, y, r, samples=samples)
    fig = plot(xs, ys, seriestype=[:shape], lw=lw, c=c, linecolor=:black, fillalpha=fillalpha)

    return fig
end

function plot_circle!(x, y, r; samples=500, lw=1.0, c=:gray70, fillalpha=1.0)
    xs, ys = circle_shape(x, y, r, samples=samples)
    fig = Plots.current()
    plot!(fig, xs, ys, seriestype=[:shape], lw=lw, c=c, linecolor=:black, fillalpha=fillalpha)
    
    return fig
end

"""
    circle_shape(x, y, r; samples=500)

Create coordinates of circle for plotting.
"""
function circle_shape(x, y, r; samples=500)
    θ = LinRange(0.0, 2*π, samples)
    xs = x .+ r*cos.(θ)
    ys = y .+ r*sin.(θ)

    return xs, ys
end