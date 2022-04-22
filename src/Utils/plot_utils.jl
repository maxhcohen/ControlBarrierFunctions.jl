# Various utility functions used to make nice plots

# Custom colors
mycolors = [
    PlotThemes.RGB255(49, 114, 174),  # blue
    PlotThemes.RGB255(68, 156, 118),  # green
    PlotThemes.RGB255(117, 112, 173), # purple
    PlotThemes.RGB255(224, 107, 97),  # red
    PlotThemes.RGB255(221, 162, 66),  # yellow
    PlotThemes.RGB255(202, 98, 159),  # magenta
    PlotThemes.RGB255(114, 179, 224), # cyan
]

# Custom theme
mytheme = PlotThemes.PlotTheme(Dict([
            :grid => false,
            :fontfamily => "Computer Modern",
            :label => "",
            :guidefontsize => 12,
            :titlefontsize => 12,
            :legendfontsize => 12,
            :tickfontsize => 8,
            :framestyle => :box,
            :linewidth => 3,
            :palette => mycolors,
            :colorgradient => :viridis,
            :foregroundcolorlegend => nothing,
            :backgroundcolorlegend => nothing,
            ]))


function custom_theme()
    return mytheme
end

function custom_plots()
    add_theme(:mytheme, custom_theme())
    theme(:mytheme)
end

"""
    circle_shape(x, y, r)
    
Constucts circle object to be used in plots.

Example usage: 
plot!(circle_shape(x, y, r), seriestype=[:shape], fillcolor=:red, fillalpha=0.2, 
        linecolor=:black, lw=3, edgecolor=:black, label="")
"""
function circle_shape(x, y, r)
    θ = LinRange(0, 2 * π, 500)
    return x .+ r * cos.(θ), y .+ r * sin.(θ)
end

function circle_shape(O::CircularObstacle)
    θ = LinRange(0, 2 * π, 500)
    return O.c[1] .+ O.r * cos.(θ), O.c[2] .+ O.r * sin.(θ)
end
