# Simulations 

## Background
We provide a shorthand function for running simulations of Control Affine Systems equipped with feedback controllers:
```julia
sol = simulate(Σ, k, x0, T)
```
where `Σ` is a `ControlAffineSystem`, `k` is a feedback controller of the form `k(x,t)` or `k(x)`, `x0` is the initial state, and `T` is the length of the simulation. This function simply calls a solver from `DifferentialEquations.jl` to simulate the closed-loop system and returns an ODE solution object.

## Implementation
```@docs
simulate(Σ::ControlAffineSystem, x0, T::Float64)
simulate(Σ::ControlAffineSystem, k::Function, x0, T::Float64)
```