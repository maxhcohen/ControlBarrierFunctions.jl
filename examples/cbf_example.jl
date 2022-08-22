# Import packages
using Revise
using CBFToolbox
using LinearAlgebra
using Plots
using LaTeXStrings
default(grid=false, framestyle=:box, fontfamily="Computer Modern", label="")

# First we need to define a control affine system

n = 2 # State dimension
m = 2 # Control dimension
f(x) = zeros(2) # Drift dynamics
g(x) = diagm(ones(2)) # Control directions
Σ = ControlAffineSystem(n, m, f, g) # Construct ControlAffineSystem

# Next we need to define CBFs - we'll consider CBFs for two circular obstacles

# CBF for first obstacle
xo = [-1.5, 1.5] # Center of obstacle
ro = 0.4 # Obstacle radius
h(x) = norm(x - xo)^2 - ro^2 # Function defining the CBF
α(s) = s^3 # Extended class K function
CBF = ControlBarrierFunction(h, α) # Construct a Control Barrier function

# Repeat same steps for the other obstacle
xo2 = [-0.7, -0.2]
ro2 = 0.4
h2(x) = norm(x - xo2)^2 - ro2^2
CBF2 = ControlBarrierFunction(h2, α)

# To reach the goal we define a CLF
V(x) = 0.5x'x # Lyapunov candidate
γ(x) = V(x) # Negative definite function defining the rate of CLF decay V̇(x) ≤ -γ(x)
CLF = ControlLyapunovFunction(V, γ) # Construct a ControlLyapunovFunction

# Now we can use the CBF and CLF to define different control policies
k0 = CLFQuadProg(Σ, CLF) # CLF-QP
k = CBFQuadProg(Σ, [CBF, CBF2], k0) # CBF-QP using the CLF-QP as a nominal policy

# Start plotting some stuff

# Vector field coordinates
xx = -3:0.2:1
yy = -1:0.2:3

# Initial conditions for phase portrait
xx_phase = -3.0:1.0:1.0
yy_phase = -1.0:1.0:3.0
T = 20.0

# Plot vector field and phase portrait
begin
    fig = plot(xlabel=L"x_1", ylabel=L"x_2", dpi=200)
    plot_phase_portrait!(xx_phase, yy_phase, Σ, k, T, lw=2)
    plot_vector_field!(xx, yy, Σ, k)
    plot_circle!(xo[1], xo[2], ro)
    plot_circle!(xo2[1], xo2[2], ro2)
    xlims!(-3.1, 1.0)
    display(fig)
end
# savefig(fig, "cbf_example.png")