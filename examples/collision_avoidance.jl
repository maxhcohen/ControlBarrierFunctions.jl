## Import packages
using CBFToolbox
using Plots
using LaTeXStrings

## Define a control affine system
n = 2
m = 2
f(x) = [0.0, 0.0]
g(x) = [1.0 0.0; 0.0 1.0]
Σ = ControlAffineSystem(n, m, f, g)

## Control bounds
U1 = [-2.0, 2.0]
U2 = [-2.0, 2.0]
U = [U1, U2]

## Define CLF
V(x) = 0.5x'x
γ(s) = s
CLF = ControlLyapunovFunction(V, γ)

## Define CBFQP
O1 = CircularObstacle([-1.0, 1.0], 0.4)
O2 = CircularObstacle([-0.5, -0.5], 0.4)
α(s) = s^3
CBF1 = ControlBarrierFunction(O1, α)
CBF2 = ControlBarrierFunction(O2, α)
CBFs = [CBF1, CBF2]
κ = CBFQP(Σ, CBFs, CLF, U)

## Construct simulation object
t0 = 0.0
tf = 10.0
dt = 0.01
sim = Simulation(t0, tf, dt)

## Run closed-loop sim
x0 = [-2.2, 2.0]
x = sim(Σ, κ, x0)

## Plot results
latexify_plots()
fig = plot(t, x', xlabel=L"t", ylabel=L"x(t)")
Plots.display(fig)

## Phase portrait
fig = plot(x[1,:], x[2,:], xlabel=L"x_1", ylabel=L"x_2")
plot!(circle_shape(O1), seriestype=[:shape], fillcolor=:red, fillalpha=0.2,
        linecolor=:black, lw=2, edgecolor=:black, label="")
plot!(circle_shape(O2), seriestype=[:shape], fillcolor=:red, fillalpha=0.2,
linecolor=:black, lw=2, edgecolor=:black, label="")
Plots.display(fig)
