"""
    ISSfCBFQuadProg <: FeedbackController

Input-to-state safe control barrier function quadratic program.

# Fields
- `solve` : function that solves the QP for the control input
- `H` : quadratic weight in QP objective
- `F` : linear weight in QP objective
"""
struct ISSfCBFQuadProg <: FeedbackController
    solve::Function
    H::Union{Float64, Matrix{Float64}, Function}
    F::Union{Float64, Vector{Float64}, Function}
    ε::Function
end

"""
    (k::ISSfCBFQuadProg)(x) = ISSfCBFQuadProg.solve(x)

Solve ISSf-CBF-QP at state `x`.
"""
(k::ISSfCBFQuadProg)(x) = k.solve(x)

"""
    ISSfCBFQuadProg(Σ::ControlAffineSystem, CBFs::Vector{ControlBarrierFunction}, ε::Function)
    ISSfCBFQuadProg(Σ::ControlAffineSystem, CBFs::Vector{ControlBarrierFunction}, ε0::Float64)
    ISSfCBFQuadProg(
        Σ::ControlAffineSystem, 
        CBFs::Vector{ControlBarrierFunction}, 
        ε0::Float64, 
        λ::Float64
        )
    ISSfCBFQuadProg(
        Σ::ControlAffineSystem, 
        CBFs::Vector{ControlBarrierFunction}, 
        k::FeedbackController, 
        ε::Function
        )
    ISSfCBFQuadProg(
        Σ::ControlAffineSystem, 
        CBFs::Vector{ControlBarrierFunction},
        k::FeedbackController, 
        ε0::Float64
        )
    ISSfCBFQuadProg(
        Σ::ControlAffineSystem, 
        CBFs::Vector{ControlBarrierFunction}, 
        k::FeedbackController, 
        ε0::Float64, 
        λ::Float64
        )

Construct an `ISSfCBFQuadProg` object.
"""
function ISSfCBFQuadProg(Σ::ControlAffineSystem, CBFs::Vector{ControlBarrierFunction}, ε::Function)
    # Set parameters for objective function
    H = Σ.m == 1 ? 1.0 : Matrix(1.0I, Σ.m, Σ.m)
    F = Σ.m == 1 ? 0.0 : zeros(Σ.m)

    # Construct quadratic program
    function solve(x)
        # Build QP and instantiate control decision variable
        model = Model(OSQP.Optimizer)
        set_silent(model)
        Σ.m == 1 ? @variable(model, u) : @variable(model, u[1:Σ.m])

        # Set CBF constraint and objective
        for CBF in CBFs
            Lfh = drift_lie_derivative(CBF, Σ, x)
            Lgh = control_lie_derivative(CBF, Σ, x)
            α = CBF.α(CBF.h(x)) - (1/ε(CBF.h(x)))*Lgh*Lgh'
            @constraint(model, Lfh + Lgh*u >= -α)
        end
        @objective(model, Min, 0.5*u'*H*u + F'*u)

        # Add control bounds on system - recall these default to unbounded controls
        if ~(Inf in Σ.b)
            @constraint(model, Σ.A * u .<= Σ.b)
        end

        # Solve QP
        optimize!(model)

        return Σ.m == 1 ? value(u) : value.(u)
    end

    return ISSfCBFQuadProg(solve, H, F, ε)
end

function ISSfCBFQuadProg(Σ::ControlAffineSystem, CBFs::Vector{ControlBarrierFunction}, ε0::Float64)
    # Construct damping function
    ε(s) = ε0

    return ISSfCBFQuadProg(Σ, CBFs, ε)
end

function ISSfCBFQuadProg(
    Σ::ControlAffineSystem, 
    CBFs::Vector{ControlBarrierFunction}, 
    ε0::Float64, 
    λ::Float64
    )
    # Construct damping function
    ε(s) = ε0*exp(λ*s)

    return ISSfCBFQuadProg(Σ, CBFs, ε)
end

function ISSfCBFQuadProg(
    Σ::ControlAffineSystem, 
    CBFs::Vector{ControlBarrierFunction}, 
    k::FeedbackController, 
    ε::Function
    )
    # Set parameters for objective function
    H = Σ.m == 1 ? 1.0 : Matrix(1.0I, Σ.m, Σ.m)
    F(x) = -H*k(x)

    # Construct quadratic program
    function solve(x)
        # Build QP and instantiate control decision variable
        model = Model(OSQP.Optimizer)
        set_silent(model)
        Σ.m == 1 ? @variable(model, u) : @variable(model, u[1:Σ.m])

        # Set CBF constraint and objective
        for CBF in CBFs
            Lfh = drift_lie_derivative(CBF, Σ, x)
            Lgh = control_lie_derivative(CBF, Σ, x)
            α = CBF.α(CBF.h(x)) - (1/ε(CBF.h(x)))*Lgh*Lgh'
            @constraint(model, Lfh + Lgh*u >= -α)
        end
        @objective(model, Min, 0.5*u'*H*u + F(x)'*u)

        # Add control bounds on system - recall these default to unbounded controls
        if ~(Inf in Σ.b)
            @constraint(model, Σ.A * u .<= Σ.b)
        end

        # Solve QP
        optimize!(model)

        return Σ.m == 1 ? value(u) : value.(u)
    end

    return ISSfCBFQuadProg(solve, H, F, ε)
end

function ISSfCBFQuadProg(
    Σ::ControlAffineSystem, 
    CBFs::Vector{ControlBarrierFunction},
    k::FeedbackController, 
    ε0::Float64
    )
    # Construct damping function
    ε(s) = ε0

    return ISSfCBFQuadProg(Σ, CBFs, k, ε)
end

function ISSfCBFQuadProg(
    Σ::ControlAffineSystem, 
    CBFs::Vector{ControlBarrierFunction}, 
    k::FeedbackController, 
    ε0::Float64, 
    λ::Float64
    )
    # Construct damping function
    ε(s) = ε0*exp(λ*s)

    return ISSfCBFQuadProg(Σ, CBFs, k, ε)
end

function ISSfCBFQuadProg(Σ::ControlAffineSystem, CBF::ControlBarrierFunction, ε::Function)
    return ISSfCBFQuadProg(Σ, [CBF], ε)
end

function ISSfCBFQuadProg(Σ::ControlAffineSystem, CBF::ControlBarrierFunction, ε::Float64)
    return ISSfCBFQuadProg(Σ, [CBF], ε)
end

function ISSfCBFQuadProg(
    Σ::ControlAffineSystem,
     CBF::ControlBarrierFunction,
      ε::Float64, 
      λ::Float64
      )
    return ISSfCBFQuadProg(Σ, [CBF], ε, λ)
end

function ISSfCBFQuadProg(
    Σ::ControlAffineSystem, 
    CBF::ControlBarrierFunction, 
    k::FeedbackController, 
    ε::Function
    )
    return ISSfCBFQuadProg(Σ, [CBF], k, ε)
end

function ISSfCBFQuadProg(
    Σ::ControlAffineSystem, 
    CBF::ControlBarrierFunction, 
    k::FeedbackController, 
    ε::Float64
    )
    return ISSfCBFQuadProg(Σ, [CBF], k, ε)
end

function ISSfCBFQuadProg(
    Σ::ControlAffineSystem, 
    CBF::ControlBarrierFunction,
    k::FeedbackController, 
    ε::Float64, 
    λ::Float64
    )
    return ISSfCBFQuadProg(Σ, [CBF], k, ε, λ)
end