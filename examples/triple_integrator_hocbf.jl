# Load in packages
using Revise
using ControlBarrierFunctions
using LinearAlgebra
using Plots
using MatrixEquations

# Create a 2D triple integrator
n = 6
m = 2
f(x) = [x[3], x[4], x[5], x[6], 0.0, 0.0]
g(x) = [0.0 0.0; 0.0 0.0; 0.0 0.0; 0.0 0.0; 1.0 0.0; 0.0 1.0]
Σ = ControlAffineSystem("triple_integrator_2d", n, m, f, g)

# Make an LQR controller for getting to a goal location
A = vcat([0.0 0.0 1.0 0.0 0.0 0.0], [0.0 0.0 0.0 1.0 0.0 0.0], [0.0 0.0 0.0 0.0 1.0 0.0], [0.0 0.0 0.0 0.0 0.0 1.0], zeros(1, 6), zeros(1, 6))
B = g(zeros(n))
Q = diagm([1.0, 1.0, 0.4, 0.4, 0.2, 0.2])
R = diagm(ones(2))
_, _, K = arec(A, B, R, Q)

# Desired controller
xd = [0.0, 0.0]
e(x) = vcat(x[1:2] - xd, x[3:6])
kd(x) = -K * e(x)

# Obstacle
xo = [-1.0, 1.0]
ro = 0.4
h(x) = norm(x[1:2] - xo)^2 - ro^2

# Make a HOCBF
r = 3
α = [1.0, 1.0, 1.0]
hocbf = HighOrderCBF(h, Σ, r, α);

# Make safety filter from HOCBF
k = ExplicitSafetyFilter(hocbf, Σ, kd);

# Run sim
T = 15.0
x0 = [-2.1, 2.0, 0.0, 0.0, 0.0, 0.0]
sol = simulate(Σ, k, x0, T)

# Set up plots
default(;
	fontfamily = "Computer Modern",
	palette = :tab10,
	framestyle = :box,
	grid = false,
	lw = 2,
	guidefont = 12,
	tickfont = 10,
	legendfont = 10,
	size = (500, 500),
	legend_background_color = nothing,
	legend_foreground_color = nothing,
);

# Plot trajectory
plot(sol; idxs = (1, 2), label = "", axis_equal = true)
contour!(
	-1.5:0.01:-0.5,
	0.5:0.01:1.5,
	(x, y) -> h([x, y]);
	levels = [0.0],
	colorbar = false,
	c = "black",
)
