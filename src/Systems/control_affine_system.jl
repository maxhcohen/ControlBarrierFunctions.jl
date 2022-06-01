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
_f(Σ::ControlAffineSystem, x) = drift(Σ, x)
_g(Σ::ControlAffineSystem, x) = actuation(Σ, x)
_F(Σ::ControlAffineSystem, x) = regressor(Σ, x)
_φ(Σ::ControlAffineSystem, x) = matched_regressor(Σ, x)
_f0(Σ::ControlAffineSystem, x) = nominal_drift(Σ, x)
closed_loop_dynamics(Σ::ControlAffineSystem, x, u) = _f(Σ, x) + _g(Σ, x)*u
dynamics(Σ::ControlAffineSystem, x, u) = closed_loop_dynamics(Σ, x, u)
fully_actuated(Σ::ControlAffineSystem) = degrees_of_freedom(Σ) == control_dim(Σ)

function dynamics(Σ::ControlAffineSystem, x)
    u = control_dim(Σ) == 1 ? 0.0 : zeros(control_dim(Σ))
    return closed_loop_dynamics(Σ, x, u)
end

function open_loop_dynamics(Σ::ControlAffineSystem, x)
    u = control_dim(Σ) == 1 ? 0.0 : zeros(control_dim(Σ))
    return closed_loop_dynamics(Σ, x, u)
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

function drift_jacobian(Σ::ControlAffineSystem, x)
    return jacobian(x -> _f(Σ, x), x)[1]
end

function linearize(Σ::ControlAffineSystem, x)
    A = jacobian(x -> _f(Σ, x), x)[1]
    B = _g(Σ, x)

    return A, B
end