"""
    ControlAffineSystem <: System

Models a nonlinear control affine system.

# Fields
- `n::Int`: state dimension
- `m::Int`: control dimension
- `f::Function`: function f(x) that models the drift dynamics
- `g::Function`: function g(x) that models the control directions
- `A`: defines control constraints of the form A*u <= b
- `b`: defines control constraints of the form A*u <= b
"""
struct ControlAffineSystem <: System
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

Construct a nonlinear control affine system from dimensions and dynamics functions. 

If no control bounds are passed in, the bounds on each control are set to infinity. If 
passed in a list of bounds on the individual control inputs, convert them 
to the form of `A*u <= b`.

# Examples
```julia-repl
julia> n = 2
2

julia> m = 1
1

julia> f(x) = [x[2], 0.0]
f (generic function with 1 method)

julia> g(x) = [0.0, 1.0]
g (generic function with 1 method)

julia> Σ = ControlAffineSystem(n, m, f, g)
ControlAffineSystem(2, 1, f, g, [1.0; -1.0;;], [Inf, Inf])
```
"""
function ControlAffineSystem(n::Int, m::Int, f::Function, g::Function)
    # Set unbounded control constraints
    Nconstraints = 2*m
    A = zeros(Nconstraints, m)
    b = Inf*ones(Nconstraints)

    # Populate entries of A matrix
    a = [1.0, -1.0]
    j = 1
    for i in 1:2:Nconstraints
        A[i:i+1, j] = a
        j += 1
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
        A[i:i+1, j] = a
        ubounds = U[j]
        b[i:i+1] = [ubounds[2], -ubounds[1]]
        j += 1
    end

    return ControlAffineSystem(n, m, f, g, A, b)
end

############################################################################################
#                                      Simulations                                         #
############################################################################################

"""
    (S::Simulation)(Σ::ControlAffineSystem, x::Union{Float64, Vector{Float64}})
    (S::Simulation)(Σ::ControlAffineSystem, x::SVector)

Run open-loop simulation of control affine system from initial state x.
"""
function (S::Simulation)(Σ::ControlAffineSystem, x::Union{Float64, Vector{Float64}})
     # Make in-place and out of place rhs functions
    rhs(x, p, t) = Σ.f(x)

    function rhs!(dx, x, p, t)
        dx .= Σ.f(x)
        nothing
    end

    problem = S.inplace ? ODEProblem(rhs!, x, S.tf) :  ODEProblem(rhs, x, S.tf)
    trajectory = solve(problem, Tsit5())

    return trajectory
end

function (S::Simulation)(Σ::ControlAffineSystem, x::SVector)
    right_hand_side(x, p, t) = Σ.f(x)
    problem = ODEProblem(right_hand_side, x, S.tf)
    trajectory = solve(problem, Tsit5())

    return trajectory
end

"""
    (S::Simulation)(Σ::ControlAffineSystem, k::FeedbackController, x::Union{Float64, Vector{Float64}})
    (S::Simulation)(Σ::ControlAffineSystem, k::FeedbackController, x::SVector)

Run closed-loop simulation of control affine system from initial state x under feedback
control policy.
"""
function (S::Simulation)(
    Σ::ControlAffineSystem,
    k::FeedbackController,
    x::Union{Float64, Vector{Float64}}
    )

    # Make in-place and out of place rhs functions
    rhs(x, p, t) = Σ.f(x) + Σ.g(x)*k(x)

    function rhs!(dx, x, p, t)
        dx .= Σ.f(x) + Σ.g(x)*k(x)
        nothing
    end

    problem = S.inplace ? ODEProblem(rhs!, x, S.tf) :  ODEProblem(rhs, x, S.tf)
    trajectory = solve(problem, Tsit5())

    return trajectory
end

function (S::Simulation)(
    Σ::ControlAffineSystem,
    k::FeedbackController,
    x::SVector,
    )
    right_hand_side(x, p, t) = SVector{Σ.n}(Σ.f(x) + Σ.g(x)*k(x))
    problem = ODEProblem(right_hand_side, x, S.tf)
    trajectory = solve(problem, Tsit5())

    return trajectory
end

"""
    (S::Simulation)(Σ::ControlAffineSystem, X::Vector{Vector{Float64}})

Run multiple open-loop simulations of control affine system from initial states X.
"""
function (S::Simulation)(Σ::ControlAffineSystem, X::Vector{Vector{Float64}})
    trajectories = [S(Σ, x) for x in X]

    return trajectories
end

"""
    (S::Simulation)(Σ::ControlAffineSystem, k::FeedbackController, X::Vector{Vector{Float64}})

Run multiple closed-loop simulations of a control affine system from initial states X.
"""
function (S::Simulation)(
    Σ::ControlAffineSystem,
    k::FeedbackController,
    X::Vector{Vector{Float64}}
    )
    trajectories = [S(Σ, k, x) for x in X]

    return trajectories
end