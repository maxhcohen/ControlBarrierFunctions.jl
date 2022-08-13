using Revise
using CBFToolbox
using LinearAlgebra
using Plots
using LaTeXStrings
default(grid=false, framestyle=:box, fontfamily="Computer Modern", label="")

# Define system
n = 2
m = 1
f(x) = [-0.6x[1] - x[2], x[1]^3]
g(x) = [0.0, x[2]]
Σ = ControlAffineSystem(n, m, f, g)

# Define CLF
V(x) = 0.25*x[1]^4 + 0.5*x[2]^2
γ(x) = V(x)
CLF = ControlLyapunovFunction(V, γ)

# Define safe set
h(x) = 1 - x[1] - x[2]^2
α(s) = s
CBF = ControlBarrierFunction(h, α)

# QP Controllers
kCLF = CLFQuadProg(Σ, CLF)
kCBF = CBFQuadProg(Σ, CBF)
kCBFCLF = CBFQuadProg(Σ, CBF, kCLF)

# Coordinates for contours
xx_contour = -4:0.1:2
yy_contour = -4:0.1:3

# Coordinates for vector field
xx_quiver = -4.0:0.3:2.0
yy_quiver = -4.0:0.3:3.0

# Coordinates for phase portraits
xx_phase = -3.0:1.5:1.0
yy_phase = -2.0:1.5:2.0
T = 10.0

# Take a look at the open-loop vector field
fig = plot_vector_field(xx_quiver, yy_quiver, Σ, scale=0.25)
plot_phase_portrait!(xx_phase, yy_phase, Σ, T, lw=2)
plot_safe_set!(xx_contour, yy_contour, CBF, lw=2)
xlabel!(L"x_1")
ylabel!(L"x_2")
title!("Open loop vector field")
xlims!(-4, 2)
ylims!(-4, 3)
display(fig)

# Take a look at the closed-loop vector field with CLF controller
fig = plot_vector_field(xx_quiver, yy_quiver, Σ, kCLF, scale=0.25)
plot_phase_portrait!(xx_phase, yy_phase, Σ, kCLF, T, lw=2)
plot_safe_set!(xx_contour, yy_contour, CBF, lw=2)
xlabel!(L"x_1")
ylabel!(L"x_2")
title!("CLF vector field")
xlims!(-4, 2)
ylims!(-4, 3)
display(fig)

# Take a look at the closed-loop vector field with CBF controller
fig = plot_vector_field(xx_quiver, yy_quiver, Σ, kCBF, scale=0.25)
plot_phase_portrait!(xx_phase, yy_phase, Σ, kCBF, T, lw=2)
plot_safe_set!(xx_contour, yy_contour, CBF, lw=2)
xlabel!(L"x_1")
ylabel!(L"x_2")
title!("CBF vector field")
xlims!(-4, 2)
ylims!(-4, 3)
display(fig)

# Take a look at the closed-loop vector field with CBF-CLF controller
fig = plot_vector_field(xx_quiver, yy_quiver, Σ, kCBFCLF, scale=0.25)
plot_phase_portrait!(xx_phase, yy_phase, Σ, kCBFCLF, T, lw=2)
plot_safe_set!(xx_contour, yy_contour, CBF, lw=2)
xlabel!(L"x_1")
ylabel!(L"x_2")
title!("CBF + CLF vector field")
xlims!(-4, 2)
ylims!(-4, 3)
display(fig)