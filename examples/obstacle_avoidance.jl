# Load in packages
using Revise
using CBFToolbox
using LinearAlgebra
using Plots

# Create a single integrator
n = 2
m = 2
f(x) = zeros(2)
g(x) = diagm(ones(2))
Σ = ControlAffineSystem("single integrator", n, m, f, g)

# Create a CBF for an obstacle
xo = [-1.0, 1.0]
ro = 0.4
h(x) = norm(x - xo)^2 - ro^2
α(r) = r
cbf = ControlBarrierFunction(h, Σ, α);

# Nominal controller
kd(x) = -x

# Create different safety filters
umin = -0.5 * ones(m)
umax = 0.5 * ones(m)
kE = ExplicitSafetyFilter(cbf, Σ, kd);
kQP = QPSafetyFilter(cbf, Σ, kd, umin, umax);
kS = SmoothSafetyFilter(cbf, Σ, kd; σ=0.01);

# Run a simulation
x0 = [-2.1, 2.0]
T = 15.0
solE = simulate(Σ, kE, x0, T)
solQP = simulate(Σ, kQP, x0, T)
solS = simulate(Σ, kS, x0, T)

# Set up plots
default(; fontfamily="Computer Modern", palette=:tab10, framestyle=:box, grid=false, lw=2)

# Plot trajectory under different controllers
plot(solE; idxs=(1, 2), label="Explicit Filter")
plot!(solQP; idxs=(1, 2), label="QP Filter")
plot!(solS; idxs=(1, 2), label="Smooth Filter")
contour!(
    -1.5:0.01:-0.5,
    0.5:0.01:1.5,
    (x, y) -> h([x, y]);
    levels=[0.0],
    colorbar=false,
    c="black",
)
