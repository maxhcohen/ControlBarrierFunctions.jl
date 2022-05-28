"""
    ControlAffineSystem <: System

Abstract ControlAffineSystem type. This abstract type will be used to derive specific 
instances of various control affine systems. At a minimum, every type that inherets from 
ControlAffineSystem should have the following methods:

 - state_dim(::ControlAffineSystem)
 - control_dim(::ControlAffineSystem)
 - drift(::ControlAffineSystem, x)
 - actuation(::ControlAffineSystem, x)
 - closed_loop_dynamics(::ControlAffineSystem, x, u)
"""
abstract type ControlAffineSystem <: System end

# Define some shorthand functions for ControlAffineSystem
_f(Σ::ControlAffineSystem, x) = drift(Σ::ControlAffineSystem, x)
_g(Σ::ControlAffineSystem, x) = actuation(Σ::ControlAffineSystem, x) 
closed_loop_dynamics(Σ::ControlAffineSystem, x, u) = _f(Σ, x) + _g(Σ, x)*u
dynamics(Σ::ControlAffineSystem, x, u) = closed_loop_dynamics(Σ::ControlAffineSystem, x, u)

function dynamics(Σ::ControlAffineSystem, x)
    u = control_dim(Σ) == 1 ? 0.0 : zeros(control_dim(Σ))
    return closed_loop_dynamics(Σ::ControlAffineSystem, x, u)
end

function open_loop_dynamics(Σ::ControlAffineSystem, x)
    u = control_dim(Σ) == 1 ? 0.0 : zeros(control_dim(Σ))
    return closed_loop_dynamics(Σ::ControlAffineSystem, x, u)
end

function integrate(Σ::ControlAffineSystem, x, u, ts)
    rhs(x, u, t) = dynamics(Σ, x, u)
    prob = ODEProblem(rhs, x, ts, u)
    sol = solve(prob, Tsit5())

    return state_dim(Σ) == 1 ? sol[end] : sol[:,end]
end

function simulate(Σ::ControlAffineSystem, x, ts)
    rhs(x, u, t) = dynamics(Σ, x)
    prob = ODEProblem(rhs, x, ts)
    sol = solve(prob, Tsit5())

    return sol
end

function linearize(Σ::ControlAffineSystem, x)
    A = jacobian(x -> _f(Σ, x), x)[1]
    B = _g(Σ, x)

    return A, B
end

# """
#     ControlAffineSystem

# Models a nonlinear control affine system.

# # Fields
# - `n::Int`: state dimension
# - `m::Int`: control dimension
# - `f`: function f(x) that models the drift dynamics
# - `g`: function g(x) that models the control directions
# - `x0`: initial condition
# - `x`: current state of system
# - `xs`: system trajectory
# """
# mutable struct ControlAffineSystem <: System
#     n::Int
#     m::Int
#     f
#     g
#     x0
#     x
#     xs
# end

# """
#     ControlAffineSystem(n::Int, m::Int, f, g)
#     ControlAffineSystem(n::Int, m::Int, f, g, x0)

# Construct a nonlinear control affine system.
# """
# ControlAffineSystem(n::Int, m::Int, f, g) = ControlAffineSystem(n, m, f, g, zeros(n), missing, missing)
# ControlAffineSystem(n::Int, m::Int, f, g, x0) = ControlAffineSystem(n, m, f, g, x0, missing, missing)

# """
#     initialize!(Σ::ControlAffineSystem)

# Set current state to initial condition and initialize system trajectory.
# """
# function initialize!(Σ::ControlAffineSystem)
#     Σ.x = Σ.x0
#     Σ.xs = [Σ.x]

#     return Σ
# end

# """
#     step!(Σ::ControlAffineSystem, u, dt::Float64)

# Integrate the dynamics of the system under control input u for dt seconds.
# """
# function step!(Σ::ControlAffineSystem, u, dt::Float64)
#     rhs(x, u, t) = Σ.f(x) + Σ.g(x)u
#     sol = solve(ODEProblem(rhs, Σ.x, (0, dt), u))
#     Σ.x = Σ.n == 1 ? sol[end] : sol[:, end]
#     push!(Σ.xs, Σ.x)

#     return Σ
# end
