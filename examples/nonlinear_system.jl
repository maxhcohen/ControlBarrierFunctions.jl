## Import necessary packages
using Revise
using CBFToolbox
using Plots; latexify_plots()
using LaTeXStrings

## Define system dynamics
n = 2
m = 1
f(x) = [-0.6x[1] - x[2], x[1]^3]
g(x) = [0.0, x[2]]
Σ = ControlAffineSystem(n, m, f, g)

## Define CLF
V(x) = (1/4)x[1]^4 + (1/2)x[2]^2
γ(s) = s
CLF = ControlLyapunovFunction(V, γ)

## Define safe set and CBF
h(x) = 1 - x[1] - x[2]^2
α(s) = s^3
CBF = ControlBarrierFunction(h, α)

## Define our control policy
κ = CBFQP(Σ, CBF, CLF)

## Construct simulation object
t0 = 0.0
tf = 10.0
dt = 0.005
sim = Simulation(t0, tf, dt)

## Run simulation
x0 = [-4.0, 1.0]
T = sim(Σ, κ, x0)

## States
fig1 = plot(xlabel=L"t", ylabel=L"x(t)")
plot!(T.t, T.x')
Plots.display(fig1)

## CBF
fig2 = plot(xlabel=L"t", ylabel=L"h(x(t))")
plot!(T.t, [h(T.x[:,i]) for i in 1:length(sim)])
hline!([0.0], ls=:dot, c=:black, label=L"h(x)=0")
Plots.display(fig2)

## Phase portrait
fig3 = plot(xlabel=L"x_1", ylabel=L"x_2")
plot!(T.x[1,:], T.x[2,:])
h(x1, x2) = 1 - x1 - x2^2
contour!(-4.5:0.1:1, -3:0.1:3, h, levels=[0.0], colorbar=false, c=:black)
Plots.display(fig3)

## Save data to a CSV file if you'd like
# df = DataFrame(t=t, x1=x[1,:], x2=x[2,:], h=[h(x[:,i]) for i in 1:length(t)])
# CSV.write("nonlinear.csv", df)
