# Misc utility functions

"""
    circle_shape(x, y, r)
    
Constucts circle object to be used in plots.
Example: 
plot!(circleshape(x, y, r), seriestype=[:shape], fillcolor=:red, fillalpha=0.2, 
        linecolor=:black, lw=3, edgecolor=:black, label="")
"""
function circle_shape(x, y, r)
    θ = LinRange(0, 2*π, 500)
    x .+ r*cos.(θ), y .+ r*sin.(θ)
end