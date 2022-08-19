# Import packages
using Revise
using CBFToolbox
using LinearAlgebra
using Plots
using LaTeXStrings
default(grid=false, framestyle=:box, fontfamily="Computer Modern", label="")

# Define system dynamics: inverted pendulum
mass = 2.0
l = 1.0
grav = 10.0
f(x) = [x[2], (grav/l)*x[1]]
g(x) = [0.0, 1/(mass*l^2)]
Σ = ControlAffineSystem(2, 1, f, g)

# Define a CLF via backstepping
cp = 3.0
k0(x) = -cp*x[1]
V(x) = x[1]^2 + 0.5*(x[2] - k0(x))^2
γ(x) = V(x)
CLF = ControlLyapunovFunction(V, γ)

# Define a standard CLF controller for comparison
k = CLFQuadProg(Σ, CLF)

# Define a colection of ISS-CLF controllers with different damping factors
εs = [5.0, 1.0, 0.5]
ISSks = [ISSCLFQuadProg(Σ, CLF, ε) for ε in εs]

# Plot corresponding vector fields
xx = LinRange(-π, π, 20)
yy = LinRange(-π, π, 20)
xx_phase = LinRange(-π, π, 5)
yy_phase = LinRange(-π, π, 5)
begin
    fig = plot_vector_field(xx, yy, Σ, k, scale=0.35)
    plot_phase_portrait!(xx_phase, yy_phase, Σ, k, 20.0, lw=3)
    xlims!(-π, π)
    xlabel!(L"x_1")
    ylabel!(L"x_2")
    title!("CLF-QP")
    display(fig)
end

begin
    for (i, k) in enumerate(ISSks)
        fig = plot_vector_field(xx, yy, Σ, k, scale=0.35)
        plot_phase_portrait!(xx_phase, yy_phase, Σ, k, 20.0, lw=3)
        xlims!(-π, π)
        xlabel!(L"x_1")
        ylabel!(L"x_2")
        title!("ISS-CLF-QP "*L"\varepsilon="*string(εs[i]))
        display(fig)
    end
end
