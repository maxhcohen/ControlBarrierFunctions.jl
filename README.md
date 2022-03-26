# CBFToolbox.jl
Control barrier functions (CBFs) and control Lyapunov functions (CLFs) written in Julia.

## Installation
To download this package open the Julia REPL, enter the package manager (type `]` into the REPL) and run

    add https://github.com/maxhcohen/CBFToolbox.jl.git

## Overview
This toolbox provides utilities to construct nonlinear systems and control policies based on control barrier functions and control Lyapunov functions. The utilities in this toolbox make heavy use of Julia's multiple dispatch functionality and are intended to provide a lightweight base for more complex packages that leverage CBFs and CLFs. This package is under active development, so things may change somewhat frequently.

The utilities in this package are built into four main types:
1. `ControlAffineSystem`
2. `ControlLyapunovFunction`
3. `ControlBarrierFunction`
4. `Policy`

The `ControlAffineSystem` type represents a nonlinear control affine system of the form

    ẋ = f(x) + g(x)u

where `x ∈ ℝⁿ` is the system state, `u ∈ ℝᵐ` is the control input, `f(x)` captures the drift dynamics,  and `g(x)` captures the control directions. A `ControlAffineSystem` can be constructed as

    Σ = ControlAffineSystem(n, m, f, g)

where `n` is an integer representing the state dimension, `m` is an integer representing the control dimension, `f(x)` is a function representing the drift dynamics, and `g(x)` is a function representing the control directions. The control input for a system can be computed by constructing a `Policy` from either a `ControlLyapunovFunction`, a `ControlBarrierFunction`, or both. A `ControlLyapunovFunction` can be constructed as 

    CLF = ControlLyapunovFunction(V, γ)

where `V(x)` is a function representing the Lyapunov function candidate and `γ(s)` is a class K function that specifies the rate of convergence of the Lyapunov function. A `ControlBarrierFunction` can be constructed in a similar fashion as

    CBF = ControlBarrierFunction(h, α)

where `h(x)` is a function defining the safe set `C  ` as `C = {x ∈ ℝⁿ | h(x) ≥ 0}` and `α` is an extended class K function that specifies how quickly the system may approach the boundary of the safe set. The `Policy` type makes heavy use of multiple dispatch. For example, the commands

    k = CLFQP(Σ, CLF)
    k = CLFQP(Σ, CLF, A, b)

construct a `Policy` represented as a CLF quadratic program without any control constraints and with control constraints specified by `Au <= b`, respectively. Similarly, the following commands

    k = CBFQP(Σ, CBF)
    k = CBFQP(Σ, CBF, CLF)

construct a `Policy` represented as a CBF quadratic program without and with an additional CLF constraint, respectively. Once a `Policy` is constructed, it can be evaluated by calling `u = k(x)`, which will solve the corresponding quadratic program for the control input. Simulating a system under a policy can be performed by first constructing a `Simulation` object as

    sim = Simulation(t0, tf, dt)

where `t0`, `tf`, and `dt` represent the start time, stop time, and time-step, respectively. A simulation can then be run by first setting the initial condition of the system as `Σ.x0 = x0` and calling

    sim(Σ, k)

which populates the field `Σ.xs` with the resulting system trajectory.

## Simple Example
```julia
# Import CBFToolbox.jl
using CBFToolbox

# Define system: single integrator
n = 2
m = 2
f(x) = [0.0, 0.0]
g(x) = [1.0 0.0; 0.0 1.0]
Σ = ControlAffineSystem(n, m, f, g)

# Define Control Lyapunov Function
V(x) = 0.5x'x
γ(s) = s
CLF = ControlLyapunovFunction(V, γ)

# Define CBF: avoid a circular obstacle
xo = [-0.5, 0.5]
ro = 0.4
h(x) = (x[1] - xo[1])^2 + (x[2]- xo[2])^2 - ro^2
α(s) = s^3
CBF = ControlLyapunovFunction(h, α)

# Construct control policy: CBF-CLF-QP
k = CBFQP(Σ, CBF, CLF)

# Create simulation object
sim = Simulation(20.0)

# Set initial condition and run
Σ.x0 = [-1.5, 1.0]
sim(Σ, k)
```

## Various to-dos
- Generate formal documentation.
- Add more concrete constructions of common systems and CBFs.
- Provide implementations of more complex examples.
- Should `ControlAffineSystem` be an abstract type?

## Questions and Contributions
If you have any questions about the toolbox, have suggestions for improvements, or would like to make your own contribution to the toolbox feel free to reach out to the repo's owner at maxcohen@bu.edu.
 
