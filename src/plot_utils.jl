# Various utility functions used to make nice plots

"""
	custom_colors()

Construct custom color scheme based on the color scheme used in 

M. J. Kochenderfer, T. A. Wheeler, and K. H. Wray, "Algorithms for Decision Making."

Returns a vector of the following RGB values:

myBlue = RGB{Float64}(49/255, 114/255, 174/255)
myRed = RGB{Float64}(224/255, 107/255, 97/255)
myGreen = RGB{Float64}(68/255, 156/255, 118/255)
myPurple = RGB{Float64}(117/255, 112/255, 173/255)
myYellow = RGB{Float64}(221/255, 162/255, 66/255)
myMagenta = RGB{Float64}(202/255, 98/255, 159/255)
myCyan = RGB{Float64}(114/255, 179/255, 224/255)
"""
function custom_colors()
	myBlue = RGB{Float64}(49/255, 114/255, 174/255)
	myRed = RGB{Float64}(224/255, 107/255, 97/255)
	myGreen = RGB{Float64}(68/255, 156/255, 118/255)
	myPurple = RGB{Float64}(117/255, 112/255, 173/255)
	myYellow = RGB{Float64}(221/255, 162/255, 66/255)
	myMagenta = RGB{Float64}(202/255, 98/255, 159/255)
	myCyan = RGB{Float64}(114/255, 179/255, 224/255)
	mycolor_palette = [
		myBlue,
		myRed,
		myGreen,
		myPurple,
		myYellow,
		myMagenta,
		myCyan,
	]

return mycolor_palette
end

"""
	custom_plots()

Set default values for all plot settings. It performs the following actions
	- turns gridlines OFF
	- sets default linewidth to 2.5
	- sets the default color palette to a custom one
	- sets the default font to LaTeX style font
	- sets the frame style to :box
	- makes the default label of any line default to nothing
"""
function custom_plots()
	default(
	grid = false,
	linewidth = 3.0,
	guidefontsize = 12.0,
	legendfontsize = 10.0,
	color_palette = custom_colors(),
	fontfamily = "Computer Modern",
	framestyle = :box,
	label = "",
	)
end

"""
    circle_shape(x, y, r)
    
Constucts circle object to be used in plots.

Example usage: 
plot!(circleshape(x, y, r), seriestype=[:shape], fillcolor=:red, fillalpha=0.2, 
        linecolor=:black, lw=3, edgecolor=:black, label="")
"""
function circle_shape(x, y, r)
    θ = LinRange(0, 2*π, 500)
    x .+ r*cos.(θ), y .+ r*sin.(θ)
end
