"""
Cannonical adaptive cruise control (ACC) problem for CBFs. Vehicle parameters are based on those in:

    A. D. Ames, J. W. Grizzle, and P. Tabuada, "Control barrier function based quadratic programs with application to adaptive cruise control," IEEE Conference on Decision and Control, 2014.

See the above paper for more details.
"""

# Load in packages
using Revise
using ControlBarrierFunctions
using LinearAlgebra
using Plots

# Parameters of adaptive cruise control dynamics
grav = 9.81 # Acceleration due to gravity [m/s^2]
M = 1650.0 # Mass of vehicle [kg]
f0 = 0.1 # Friction coefficient [N]
f1 = 5.0 # Friction coefficient [Nm/s]
f2 = 0.25 # Friction coefficient [Nm/s^2]
vd = 24.0 # Desired velocity [m/s]
v0 = 13.89 # Velocity of lead vehicle [m/s]

# ACC dynamics: x = [v, z] = [velocity of ego, distance to lead vehicle]
f(x) = [-(1 / M) * (f0 + f1 * x[1] + f2 * x[1]^2), v0 - x[1]]
g(x) = [1 / M, 0.0]
Σ = ControlAffineSystem("acc", 2, 1, f, g);

# Desired controller: drive to desired velocity
kd(x) = -1.0 * M * (x[1] - vd)

# Barrier function
h(x) = x[2] - 1.8 * x[1]
α(r) = r
cbf = ControlBarrierFunction(h, Σ, α);

# Construct safety filter
k = ExplicitSafetyFilter(cbf, Σ, kd);

# Simulation parameters
x0 = [20.0, 100.0]
T = 25.0

# Run Simulation
sol = simulate(Σ, k, x0, T)

# Set up plots
default(;
    fontfamily="Computer Modern",
    palette=:tab10,
    framestyle=:box,
    grid=false,
    lw=2,
    guidefont=12,
    tickfont=10,
    legendfont=10,
    size=(500, 500),
    legend_background_color=nothing,
    legend_foreground_color=nothing,
);

# Plot velocity
fig1 = plot(
    sol;
    idxs=1,
    xlabel=raw"$t$",
    ylabel=raw"$v(t)$",
    label="Velocity",
    xlims=(0, T),
    ylims=(13.0, 25.0),
)
hline!([vd]; c=:black, ls=:dash, label="Desired Velocity")

# Plot value of barrier function
ts = 0.0:0.02:T
fig2 = plot(
    ts,
    h.(sol.(ts));
    c=2,
    xlabel=raw"$t$",
    ylabel="Barrier Function",
    label=raw"$h(\mathbf{x}(t))$",
    xlims=(0, T),
)
hline!([0.0]; c=:black, ls=:dash, label=raw"$h(\mathbf{x})=0$")

# Make into one plot
fig = plot(fig1, fig2; layout=(2, 1))