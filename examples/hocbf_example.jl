using Revise
using CBFToolbox
using LinearAlgebra
using Plots
using LaTeXStrings
julia_palette = deleteat!(distinguishable_colors(10, [c for c in palette(:julia)]), 5:6)
default(fontfamily="Computer Modern", grid=false, framestyle=:box, lw=2, label="", palette=julia_palette)

# Define system
n = 4
m = 2
f(x) = vcat(x[3:4], zeros(2))
g(x) = vcat(zeros(2,2), diagm(ones(2)))
Σ = ControlAffineSystem(n, m, f, g)

# Define nominal controller
Kp = 1.0
Kd = 1.0
pd_controller(x) = -Kp*x[1:2] - Kd*x[3:4]
k0 = StateFeedbackController(pd_controller)

# Define safety constraint 
xo = [-1.0, 1.0]
ro = 0.4
h(x) = norm(x[1:2] - xo)^2 - ro^2

# Construct HOCBF
α1(s) = s
α2(s) = 5s^3
HOCBF = SecondOrderCBF(Σ, h, α1, α2)
k = CBFQuadProg(Σ, HOCBF, k0)

# Run simulation
x0 = vcat([-2.1, 2.0], zeros(2))
T = 15.0
S = Simulation(T)
x = S(Σ, k, x0)

# Plot results
fig = plot(x, idxs=(1,2), lw=2, label="")
plot_circle!(xo[1], xo[2], ro)
xlabel!(L"x_1")
ylabel!(L"x_2")