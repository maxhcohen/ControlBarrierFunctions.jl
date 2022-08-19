"""
    ISSCLFQuadProg <: FeedbackController

Input-to-state stabilizing control Lyapunov function quadratic program.

# Fields
- `solve` : function that solves the QP for the control input
- `H` : quadratic weight in QP objective
- `F` : linear weight in QP objective
"""
struct ISSCLFQuadProg <: FeedbackController
    solve::Function
    H::Union{Float64, Matrix{Float64}, Function}
    F::Union{Float64, Vector{Float64}, Function}
end

"""
    (k::ISSCLFQuadProg)(x) = k.solve(x)

Solve ISS-CLF-QP at state `x`.
"""
(k::ISSCLFQuadProg)(x) = k.solve(x)

"""
    ISSCLFQuadProg(Σ::ControlAffineSystem, CLF::ControlLyapunovFunction, ε::Float64)

Construct an `ISSCLFQuadProg` where `ε` is the nonlinear damping coefficient for ISS.
"""
function ISSCLFQuadProg(Σ::ControlAffineSystem, CLF::ControlLyapunovFunction, ε::Float64)
    # Set parameters for objective function
    H = Σ.m == 1 ? 1.0 : Matrix(1.0I, Σ.m, Σ.m)
    F = Σ.m == 1 ? 0.0 : zeros(Σ.m)

    # Construct quadratic program
    function solve(x)
        # Build QP and instantiate control decision variable
        model = Model(OSQP.Optimizer)
        set_silent(model)
        Σ.m == 1 ? @variable(model, u) : @variable(model, u[1:Σ.m])

        # Compute Lie derivatives
        LfV = drift_lie_derivative(CLF, Σ, x)
        LgV = control_lie_derivative(CLF, Σ, x)
        γ = CLF.α(x) + (1/ε)*LgV*LgV'

        # Check if we're relaxing the CLF constraint
        if CLF.relax
            @variable(model, δ)
            @constraint(model, LfV + LgV*u <= -γ + δ)
            @objective(model, Min, 0.5*u'*H*u + F'*u + CLF.p*δ^2)
        else
            @constraint(model, LfV + LgV*u <= -γ)
            @objective(model, Min, 0.5*u'*H*u + F'*u)
        end

        # Add control bounds on system - recall these default to unbounded controls
        if ~(Inf in Σ.b)
            @constraint(model, Σ.A * u .<= Σ.b)
        end

        # Solve QP
        optimize!(model)

        return Σ.m == 1 ? value(u) : value.(u)
    end

    return ISSCLFQuadProg(solve, H, F)
end