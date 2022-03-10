## Import packages
using Revise
using CBFToolbox
using Plots;
latexify_plots();
using LaTeXStrings

## Construct a control affine system
n = 4
m = 2
f(x) = [x[4] * cos(x[3]), x[4] * sin(x[3]), 0.0, 0.0]
g(x) = [0.0 0.0; 0.0 0.0; 1.0 0.0; 0.0 1.0]
Σ = ControlAffineSystem(n, m, f, g)

## Define CLF
xd = 45.0
yd = 25.1
vd = 2.0
V(x) = (x[3] - atan((yd - x[2]) / (xd - x[1])))^2 + (x[4] - vd)^2
γ(s) = s
CLF = ControlLyapunovFunction(V, γ)
κCLF = CLFQP(Σ, CLF)

## Construct CBF
O = CircularObstacle([32.0, 25.0], 3.0)
α(s) = 0.5s
CBF = ControlBarrierFunction(O, α)
HOCBF = HighOrderCBF(CBF, Σ, 2)

## Construct our control policy
κ = CBFQP(Σ, HOCBF, κCLF)

## Construct simulation object
t0 = 0.0
tf = 30.0
dt = 0.01
sim = Simulation(t0, tf, dt)

## Run simulation
x0 = [5.0, 25.0, 0.0, 0.0]
T = sim(Σ, κ, x0)

## Plot results
fig1 = plot(; xlabel=L"t", ylabel=L"x(t)")
plot!(T.t, T.x')
Plots.display(fig1)

fig2 = plot(T.x[1, :], T.x[2, :]; xlabel=L"x_1", ylabel=L"x_2", aspect_ratio=:equal)
plot!(
    circle_shape(O);
    seriestype=[:shape],
    fillcolor=:red,
    fillalpha=0.2,
    linecolor=:black,
    lw=2,
    edgecolor=:black,
    label="",
)
scatter!([xd], [yd]; ms=8, c=3)
Plots.display(fig2)
