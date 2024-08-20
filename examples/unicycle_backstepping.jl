"""
Using backstepping to generate a CBF for a unicycle. Example is based on that in:

    A. J. Taylor, P. Ong, T. G. Molnar, and A. D. Ames, "Safe Backstepping with Control Barrier Functions," IEEE Conference on Decision and Control, 2022.

See the above paper for more details.
"""

# Load in packages
using Revise
using ControlBarrierFunctions
using LinearAlgebra
using Plots

# Create a unicycle
n = 3
m = 2
f(x) = zeros(3)
g(x) = [cos(x[3]) 0.0; sin(x[3]) 0.0; 0.0 1.0]
Σ = ControlAffineSystem("unicycle", n, m, f, g);

# Define safety constraint
xo = [-1.0, 1.0]
ro = 0.4
h0(x) = norm(x[1:2] - xo)^2 - ro^2

# Use a single integrator to generate safe velocities for unicycle
Σ0 = ControlAffineSystem("single integrator", 2, 2, x -> zeros(2), x -> diagm(ones(2)));

# Make CBF for single integrator reduced-order model
α(r) = r
cbf0 = ControlBarrierFunction(h0, Σ0, α);

# Get smooth safety filter for reduced-order model
kd0(x) = -x
k0 = SmoothSafetyFilter(cbf0, Σ0, kd0);

# CBF for unicycle: attempts to align unicycle's heading with safe velocity
h(x) = h0(x) - 0.5 * norm([cos(x[3]), sin(x[3])] - k0(x[1:2]))^2
cbf = ControlBarrierFunction(h, Σ, α);

# Desired controller for full-order system
Kp = 0.2
Kψ = 3
qd = [0.0, 0.0]
k0norm(x) = normalize(k0(x)) # Convert safe velocity into safe speed
ψ0(x) = atan(k0norm(x)[2], k0norm(x)[1]) # Convert safe velocity into safe heading
kd(x) = [Kp * norm(x[1:2] - qd), -Kψ * (sin(x[3]) - sin(ψ0(x[1:2])))]

# Safety filter for unicycle
k = ExplicitSafetyFilter(cbf, Σ, kd);

# Run sim from different initial headings
q0 = [-2.1, 2.1] # Initial position
θs = range(-2 * π / 3, 2 * π / 3, 30) # Range of initial headings
T = 20.0
sols = []
for θ0 in θs
    push!(sols, simulate(Σ, k, [q0; θ0], T))
end

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

# Plot trajectory under different controllers
begin
    plot()
    for sol in sols
        plot!(sol; idxs=(1, 2), label=raw"", xlabel=raw"$x$", ylabel=raw"$y$")
    end
    scatter!([q0[1]], [q0[2]]; label=raw"$x_0$", c=2, ms=6)
    scatter!([0], [0]; label=raw"Goal", c=3, m=:square, ms=6)
    contour!(
        -1.5:0.01:-0.5,
        0.5:0.01:1.5,
        (x, y) -> h0([x, y]);
        levels=[0.0],
        colorbar=false,
        c="black",
    )
    xlims!(-2.5, 0.5)
    ylims!(-0.5, 2.5)
    plot!(; aspect_ratio=:equal)
end
