"""
Planar segway example from:

    T. G. Molnar, R. K. Cosner, A. W. Singletary, W. Ubellacker, 
    and A. D. Ames, "Model-free safety-critical control for robotic systems," 
    IEEE Robotics and Automation Letters, 2022.

See the above paper for more details.
"""

# Load in packages
using Revise
using ControlBarrierFunctions
using LinearAlgebra
using Plots

# Segway parameters
grav = 9.81
R = 0.195
M = 2 * 2.485
Jc = 2 * 0.0559
L = 0.169
m = 44.798
Jg = 3.836
m0 = 52.710
J0 = 5.108
Km = 2 * 1.262
bt = 2 * 1.225

# Dynamics of Segway in Euler-Lagrange form
D(q) = [m0 m*L*cos(q[2]); m*L*cos(q[2]) J0]
function H(q, q̇)
    return [
        -m * L * sin(q[2]) * q̇[2] + bt * (q̇[1] - R * q̇[2]) / R,
        -m * grav * L * sin(q[2]) - bt * (q̇[1] - R * q̇[2]),
    ]
end
B(q) = [Km / R, -Km]

# Convert to control affine form
function f(x)
    q, q̇ = x[1:2], x[3:4]
    return [q̇; -D(q) \ H(q, q̇)]
end
function g(x)
    q, q̇ = x[1:2], x[3:4]
    return [zeros(2); D(q) \ B(q)]
end
Σ = ControlAffineSystem("segway", 4, 1, f, g)

# Safety constraint: avoid wall located at pmax
pmax = 2.0
h(x) = pmax - x[1]

# Create single integrator reduced-order model
Σ0 = ControlAffineSystem("single integrator", 1, 1, x -> 0.0, x -> 1.0);

# Desired input for reduced-order model: move forward at 1 m/s
kd(x) = 1.0

# Safety filter for reduced-order model
cbf = ControlBarrierFunction(h, Σ0, s -> 0.5 * s)
k0 = ExplicitSafetyFilter(cbf, Σ0, kd);

# Tracking controller for full-order model
k(x) = 50.0 * (x[3] - k0(x[1])) + 150.0 * x[2] + 40.0 * x[4]

# Simulate full-order model
x0 = [0.0, -0.138, 0.0, 0.0]
T = 10.0
ts = 0.0:0.01:T
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

# Plot position of segway
fig1 = plot(sol; idxs=1, c=1, xlabel=raw"$t$", ylabel=raw"$p(t)$", label="")
hline!([pmax]; c=:black, ls=:dash, label=raw"$p_{\max}$")

# Plot pitch of segway
fig2 = plot(sol; idxs=2, c=2, xlabel=raw"$t$", ylabel=raw"$\varphi(t)$", label="")

# Plot velocity of segway
fig3 = plot(sol; idxs=3, c=3, xlabel=raw"$t$", ylabel=raw"$\dot{p}(t)$", label="Actual")
plot!(ts, k0.(sol.(ts, idxs=1)); c=3, ls=:dash, label="Commanded", alpha=0.5)

# Plot pitch rate
fig4 = plot(sol; idxs=4, c=4, xlabel=raw"$t$", ylabel=raw"$\dot{\varphi}(t)$", label="")

# Put all the plots together
fig = plot(fig1, fig2, fig3, fig4)