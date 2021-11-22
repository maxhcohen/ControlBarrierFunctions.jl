"Abstract type used to derive special classes of control systems"
abstract type ControlSystem end

"""
    ControlAffineSystem

Models a nonlinear control affine system of the form 
    x' = f(x) + g(x)u
# Fields
- `n::Int`: state dimension
- `m::Int`: control dimension
- `f`: function f(x) that models the drift dynamics
- `g`: function g(x) that models the control directions
- `rhs`: function rhs(x,u,t) that returns the right hand side f(x) + g(x)u at time t
"""
struct ControlAffineSystem <: ControlSystem
    n::Int
    m::Int
    f
    g
    rhs
end

"""
    ControlAffineSystem(n::Int, m::Int, f, g)

Construct a ControlAffineSystem given system dimension and dynamics.
"""
function ControlAffineSystem(n::Int, m::Int, f, g)
    rhs(x,u,t) = f(x) + g(x)*u
    return ControlAffineSystem(n, m, f, g, rhs)
end

"""
    simulate(x, u, ts, Σ::ControlSystem)
    
Integrate the dynamics of a ControlSystem starting at state x under control u over ts.
This will be used to simulate the system forward one time-step within a larger sim.
"""
function simulate(x, u, ts, Σ::ControlSystem)
    sol = DifferentialEquations.solve(ODEProblem(Σ.rhs, x, (ts[1], ts[end]), u), saveat=ts)
    return sol[:,end]
end