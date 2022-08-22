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
knom(x) = -cp*x[1]
V(x) = x[1]^2 + 0.5*(x[2] - knom(x))^2
γ(x) = V(x)
CLF = ControlLyapunovFunction(V, γ)

# Define a standard CLF controller and ISS-CLF controller for comparison
kCLF = CLFQuadProg(Σ, CLF)
εCLF = 1.0
kISS = ISSCLFQuadProg(Σ, CLF, εCLF)

# Define safety constraints
h1(x) = x[1] + π/4
h2(x) = -x[1] + π/4

# Construct HOCBFs
α1(s) = s
α2(s) = s
HOCBF1 = SecondOrderCBF(Σ, h1, α1, α2)
HOCBF2 = SecondOrderCBF(Σ, h2, α1, α2)
HOCBFs = [HOCBF1, HOCBF2]

# Define an ISSf controller
ε0 = 1.0
λ = 0.1
kISSf = ISSfCBFQuadProg(Σ, HOCBFs, kISS, ε0, λ)

# Simulate system
x0 = [-π/6, 0.0]
T = 15.0
sim = Simulation(T)
x = sim(Σ, kISSf, x0)

# Plot corresponding vector fields
xx = LinRange(-π/2, π/2, 20)
yy = LinRange(-π/2, π/2, 20)
xx_contour = LinRange(-π/2, π/2, 40)
yy_contour = LinRange(-π/2, π/2, 40)
begin
    fig = plot_vector_field(xx, yy, Σ, kISSf, scale=0.15)
    plot_safe_set!(xx_contour, yy_contour, HOCBF1)
    plot_safe_set!(xx_contour, yy_contour, HOCBF2)
    plot!(x, idxs=(1,2), label="", c=1, lw=4)
    xlabel!(L"x_1")
    ylabel!(L"x_2")
    xlims!(-π/2, π/2)
    ylims!(-π/2, π/2)
    display(fig)
end
