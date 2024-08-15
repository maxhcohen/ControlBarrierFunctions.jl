[![Build status (Github Actions)](https://github.com/maxhcohen/CBFToolbox.jl/workflows/CI/badge.svg)](https://github.com/maxhcohen/CBFToolbox.jl/actions)
[![codecov.io](http://codecov.io/github/maxhcohen/CBFToolbox.jl/coverage.svg?branch=main)](http://codecov.io/github/maxhcohen/CBFToolbox.jl?branch=main)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://maxhcohen.github.io/CBFToolbox.jl/dev)

# CBFToolbox.jl
A package for implementing control barrier functions (CBFs) in Julia.

## Overview
This toolbox provides utilities to construct nonlinear systems and control policies based on control barrier functions (CBFs).

## Installation
To download this package open the Julia REPL, enter the package manager (type `]` into the REPL) and run

```julia
using Pkg
Pkg.add("https://github.com/maxhcohen/CBFToolbox.jl")
```

## Usage
```julia
# Load in packages
using CBFToolbox
using LinearAlgebra
using Plots

# Create a single integrator
n = 2
m = 2
f(x) = zeros(2)
g(x) = diagm(ones(2))
Σ = ControlAffineSystem("single integrator", n, m, f, g)

# Create a CBF for an obstacle
xo = [-1.0, 1.0]
ro = 0.4
h(x) = norm(x - xo)^2 - ro^2
α(r) = r
cbf = ControlBarrierFunction(h, Σ, α);

# Nominal controller
kd(x) = -x

# Create a safety filter
k = ExplicitSafetyFilter(cbf, Σ, kd);

# Run a simulation
x0 = [-2.1, 2.0]
T = 15.0
sol = simulate(Σ, k, x0, T)

# Set up plots
default(fontfamily="Computer Modern", palette=:tab10, framestyle=:box, grid=false, lw=2)

# Plot trajectory
plot(sol, idxs=(1,2), label="")
contour!(-1.5:0.01:-0.5, 0.5:0.01:1.5, (x,y) -> h([x,y]), levels=[0.0], colorbar=false, c="black")
```
 
