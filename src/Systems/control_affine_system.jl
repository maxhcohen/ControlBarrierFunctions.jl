"""
    ControlAffineSystem

Models a nonlinear control affine system.

# Fields
- `n::Int`: state dimension
- `m::Int`: control dimension
- `f`: function f(x) that models the drift dynamics
- `g`: function g(x) that models the control directions
- `x0`: initial condition
- `x`: current state of system
- `xs`: system trajectory
"""
mutable struct ControlAffineSystem <: System
    n::Int
    m::Int
    f
    g
    x0
    x
    xs
end

"""
    ControlAffineSystem(n::Int, m::Int, f, g)
    ControlAffineSystem(n::Int, m::Int, f, g, x0)

Construct a nonlinear control affine system.
"""
ControlAffineSystem(n::Int, m::Int, f, g) = ControlAffineSystem(n, m, f, g, zeros(n), missing, missing)
ControlAffineSystem(n::Int, m::Int, f, g, x0) = ControlAffineSystem(n, m, f, g, x0, missing, missing)

"""
    initialize!(Σ::ControlAffineSystem)

Set current state to initial condition and initialize system trajectory.
"""
function initialize!(Σ::ControlAffineSystem)
    Σ.x = Σ.x0
    Σ.xs = [Σ.x]

    return Σ
end

"""
    step!(Σ::ControlAffineSystem, u, dt::Float64)

Integrate the dynamics of the system under control input u for dt seconds.
"""
function step!(Σ::ControlAffineSystem, u, dt::Float64)
    rhs(x, u, t) = Σ.f(x) + Σ.g(x)u
    sol = solve(ODEProblem(rhs, Σ.x, (0, dt), u))
    Σ.x = Σ.n == 1 ? sol[end] : sol[:, end]
    push!(Σ.xs, Σ.x)

    return Σ
end
