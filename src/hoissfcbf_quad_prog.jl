function ISSfCBFQuadProg(Σ::ControlAffineSystem, HOCBFs::Vector{SecondOrderCBF}, ε::Function)
    # Set parameters for objective function
    H = Σ.m == 1 ? 1.0 : Matrix(1.0I, Σ.m, Σ.m)
    F = Σ.m == 1 ? 0.0 : zeros(Σ.m)

    # Construct quadratic program
    function solve(x)
        # Build QP and instantiate control decision variable
        model = Model(OSQP.Optimizer)
        set_silent(model)
        Σ.m == 1 ? @variable(model, u) : @variable(model, u[1:Σ.m])

        # Set HOCBF constraint and objective
        for HOCBF in HOCBFs
            Lfψ = drift_lie_derivative(HOCBF, Σ, x)
            Lgψ = control_lie_derivative(HOCBF, Σ, x)
            α = HOCBF.α2(HOCBF.ψ1(x)) - (1/ε(HOCBF.ψ1(x)))*Lgψ*Lgψ'
            @constraint(model, Lfψ + Lgψ*u >= -α)
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

function ISSfCBFQuadProg(Σ::ControlAffineSystem, HOCBFs::Vector{SecondOrderCBF}, ε0::Float64)
    # Construct damping function
    ε(s) = ε0

    return ISSfCBFQuadProg(Σ, HOCBFs, ε)
end

function ISSfCBFQuadProg(
    Σ::ControlAffineSystem, 
    HOCBFs::Vector{SecondOrderCBF},
    ε0::Float64, 
    λ::Float64
    )
    # Construct damping function
    ε(s) = ε0*exp(λ*s)

    return ISSfCBFQuadProg(Σ, HOCBFs, ε)
end

function ISSfCBFQuadProg(
    Σ::ControlAffineSystem, 
    HOCBFs::Vector{SecondOrderCBF},
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

        # Set HOCBF constraint and objective
        for HOCBF in HOCBFs
            Lfψ = drift_lie_derivative(HOCBF, Σ, x)
            Lgψ = control_lie_derivative(HOCBF, Σ, x)
            α = HOCBF.α2(HOCBF.ψ1(x)) - (1/ε(HOCBF.ψ1(x)))*Lgψ*Lgψ'
            @constraint(model, Lfψ + Lgψ*u >= -α)
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
    HOCBFs::Vector{SecondOrderCBF},
    k::FeedbackController, 
    ε0::Float64
    )
    # Construct damping function
    ε(s) = ε0

    return ISSfCBFQuadProg(Σ, HOCBFs, k, ε)
end

function ISSfCBFQuadProg(
    Σ::ControlAffineSystem, 
    HOCBFs::Vector{SecondOrderCBF},
    k::FeedbackController, 
    ε0::Float64, 
    λ::Float64
    )
    # Construct damping function
    ε(s) = ε0*exp(λ*s)

    return ISSfCBFQuadProg(Σ, HOCBFs, k, ε)
end

function ISSfCBFQuadProg(Σ::ControlAffineSystem, HOCBF::SecondOrderCBF, ε::Function)
    return ISSfCBFQuadProg(Σ, [HOCBF], ε)
end

function ISSfCBFQuadProg(Σ::ControlAffineSystem, HOCBF::SecondOrderCBF, ε::Float64)
    return ISSfCBFQuadProg(Σ, [HOCBF], ε)
end

function ISSfCBFQuadProg(
    Σ::ControlAffineSystem,
    HOCBF::SecondOrderCBF,
    ε::Float64, 
    λ::Float64
    )
    return ISSfCBFQuadProg(Σ, [HOCBF], ε, λ)
end

function ISSfCBFQuadProg(
    Σ::ControlAffineSystem, 
    HOCBF::SecondOrderCBF, 
    k::FeedbackController, 
    ε::Function
    )
    return ISSfCBFQuadProg(Σ, [HOCBF], k, ε)
end

function ISSfCBFQuadProg(
    Σ::ControlAffineSystem, 
    HOCBF::SecondOrderCBF, 
    k::FeedbackController, 
    ε::Float64
    )
    return ISSfCBFQuadProg(Σ, [HOCBF], k, ε)
end

function ISSfCBFQuadProg(
    Σ::ControlAffineSystem, 
    HOCBF::SecondOrderCBF,
    k::FeedbackController, 
    ε::Float64, 
    λ::Float64
    )
    return ISSfCBFQuadProg(Σ, [HOCBF], k, ε, λ)
end