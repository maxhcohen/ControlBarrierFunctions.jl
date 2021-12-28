## Import necessary packages
using Revise
using CBFToolbox
using LinearAlgebra
using Plots
using LaTeXStrings

## Define system dynamics - simplified adaptive cruise control
n = 2
m = 1
v0 = 13.89
f(x) = [v0 - x[2], 0.0]
g(x) = [0.0, 1.0]
Σ = ControlAffineSystem(n, m, f, g)

## Define safe set and CBF
h(x) = x[1] - 1.8x[2]
CBF = ControlBarrierFunction(h, α=r->0.5r)

## Define nominal control policy
vd = 24
k(x) = -(x[2] - vd)

## Run simulation
t0 = 0.0
tf = 25.0
dt = 0.01
x0 = [100.0, 20.0]
t, x = run_sim(t0, tf, dt, x0, k, Σ, CBF)

## Plot results
custom_plots()

## Velocity
fig1 = plot(xlabel=L"t", ylabel=L"v(t)")
plot!(t, x[2,:])
hline!([vd], ls=:dot, c=:black, label=L"v_d")
Plots.display(fig1)

## CBF
fig2 = plot(xlabel=L"t", ylabel=L"h(x(t))")
plot!(t, [h(x[:,i]) for i in 1:length(t)])
hline!([0.0], ls=:dot, c=:black, label=L"h(x)=0")
Plots.display(fig2)
