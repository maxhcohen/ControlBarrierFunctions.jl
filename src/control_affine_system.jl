"""
    ControlAffineSystem <: System

Models a nonlinear control affine system.

# Fields
- `n::Int`: state dimension
- `m::Int`: control dimension
- `f`: function f(x) that models the drift dynamics
- `g`: function g(x) that models the control directions
- `A`: defines control constraints of the form A*u <= b
- `b`: defines control constraints of the form A*u <= b
"""
mutable struct ControlAffineSystem <: System
    n::Int
    m::Int
    f::Function
    g::Function
    A
    b
end

"""
    ControlAffineSystem(n::Int, m::Int, f::Function, g::Function)
    ControlAffineSystem(n::Int, m::Int, f::Function, g::Function, U)

Construct a nonlinear control affine system. If no control bounds are passed in, set the
bounds on each control to infinity. If passed in a list of bounds on the individual control
inputs, convert them to the form of A*u <= b.
"""
function ControlAffineSystem(n::Int, m::Int, f::Function, g::Function)
    # Set unbounded control constraints
    Nconstraints = 2*m
    A = zeros(Nconstraints, m)
    b = Inf*ones(Nconstraints)

    # Populate entries of A matrix
    a = [1.0, -1.0]
    for i in 1:2:Nconstraints
        A[i:i+1, i] = a
    end

    return ControlAffineSystem(n, m, f, g, A, b)
end

function ControlAffineSystem(n::Int, m::Int, f::Function, g::Function, U)
    # Convert list of control constraints to standard inequality contraints
    Nconstraints = 2*m
    A = zeros(Nconstraints, m)
    b = zeros(Nconstraints)

    # Populate entries of A matrix and b vector
    a = [1.0, -1.0]
    j = 1
    for i in 1:2:Nconstraints
        A[i:i+1, i] = a
        ubounds = U[j]
        b[j:j+1] = [ubounds[2], ubounds[1]]
    end

    return ControlAffineSystem(n, m, f, g, A, b)
end