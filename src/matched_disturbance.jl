"""
    MatchedDisturbance <: Disturbance

# Fields
- d::Function : function `d(x,t)` defining the disturbance
- dmax::Float64 : maximum bound on disturbance of the form `||d|| <= dmax`
- A : captures halfspace constraints on the disturbance of the form `A*d <= b`
- b : captures halfspace constraints on the disturbance of the form `A*d <= b`
"""
struct MatchedDisturbance <: Disturbance
    d::Function
    dmax::Float64
    A
    b
end

(D::MatchedDisturbance)(x, t) = D.d(x, t)

function MatchedDisturbance(Σ::ControlAffineSystem, d::Function)
    N = 2*Σ.m
    dmax = Inf
    A = zeros(N, Σ.m)
    b = Inf*ones(N)

    # Populate entries of A matrix
    a = [1.0, -1.0]
    j = 1
    for i in 1:2:N
        A[i:i+1, j] = a
        j += 1
    end

    return MatchedDisturbance(d, dmax, A, b)
end

function (S::Simulation)(
    Σ::ControlAffineSystem, 
    d::MatchedDisturbance, 
    x::Union{Float64, Vector{Float64}}
    )
    right_hand_side(x, p, t) = Σ.f(x) + Σ.g(x)*d(x,t)
    problem = ODEProblem(right_hand_side, x, [S.t0, S.tf])
    trajectory = solve(problem)

    return trajectory
end

function (S::Simulation)(
    Σ::ControlAffineSystem,
    k::FeedbackController,
    d::MatchedDisturbance, 
    x::Union{Float64, Vector{Float64}}
    )
    right_hand_side(x, p, t) = Σ.f(x) + Σ.g(x)*(k(x) + d(x,t))
    problem = ODEProblem(right_hand_side, x, [S.t0, S.tf])
    trajectory = solve(problem)

    return trajectory
end

function (S::Simulation)(
    Σ::ControlAffineSystem, 
    d::MatchedDisturbance, 
    X::Vector{Vector{Float64}}
    )

    return [S(Σ, d, x) for x in X]
end

function (S::Simulation)(
    Σ::ControlAffineSystem,
    k::FeedbackController,
    d::MatchedDisturbance,
    X::Vector{Vector{Float64}}
    )
    
    return [S(Σ, k, d, x) for x in X]
end


