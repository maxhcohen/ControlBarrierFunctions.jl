"""
    CBFQuadProg <: FeedbackController

Control Barrier Function (CBF)-based quadratic program (QP) to compute control inputs
for a control affine system

# Fields
- `solve` : function that solves the QP for the control input
- `H` : quadratic weight in QP objective
- `F` : linear weight in QP objective
"""
struct CBFQuadProg <: FeedbackController
    solve::Function
    H::Union{Float64, Matrix{Float64}, Function}
    F::Union{Float64, Vector{Float64}, Function}
end

"Solve `CBFQuadProg` at state `x`."
(k::CBFQuadProg)(x) = k.solve(x)

function CBFQuadProg(Σ::ControlAffineSystem, CBFs::Vector{ControlBarrierFunction})
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
            α = CBF.α(CBF.h(x))
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

    return CBFQuadProg(solve, H, F)
end

function CBFQuadProg(
    Σ::ControlAffineSystem, 
    CBFs::Vector{ControlBarrierFunction}, 
    k::FeedbackController
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
            α = CBF.α(CBF.h(x))
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

    return CBFQuadProg(solve, H, F)
end

function CBFQuadProg(
    Σ::ControlAffineSystem, 
    CBFs::Vector{ControlBarrierFunction}, 
    k::FeedbackController,
    H::Union{Float64, Matrix{Float64}}
    )
    # Set parameters for objective function
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
            α = CBF.α(CBF.h(x))
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

    return CBFQuadProg(solve, H, F)
end

function CBFQuadProg(
    Σ::ControlAffineSystem, 
    CBFs::Vector{ControlBarrierFunction}, 
    k::FeedbackController,
    H::Function
    )
    # Set parameters for objective function
    F(x) = -H(x)*k(x)

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
            α = CBF.α(CBF.h(x))
            @constraint(model, Lfh + Lgh*u >= -α)
        end
        @objective(model, Min, 0.5*u'*H(x)*u + F(x)'*u)

        # Add control bounds on system - recall these default to unbounded controls
        if ~(Inf in Σ.b)
            @constraint(model, Σ.A * u .<= Σ.b)
        end

        # Solve QP
        optimize!(model)

        return Σ.m == 1 ? value(u) : value.(u)
    end

    return CBFQuadProg(solve, H, F)
end

function CBFQuadProg(
    Σ::ControlAffineSystem, 
    CBFs::Vector{ControlBarrierFunction}, 
    H::Function,
    F::Function
    )
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
            α = CBF.α(CBF.h(x))
            @constraint(model, Lfh + Lgh*u >= -α)
        end
        @objective(model, Min, 0.5*u'*H(x)*u + F(x)'*u)

        # Add control bounds on system - recall these default to unbounded controls
        if ~(Inf in Σ.b)
            @constraint(model, Σ.A * u .<= Σ.b)
        end

        # Solve QP
        optimize!(model)

        return Σ.m == 1 ? value(u) : value.(u)
    end

    return CBFQuadProg(solve, H, F)
end

function CBFQuadProg(
    Σ::ControlAffineSystem, 
    CBFs::Vector{ControlBarrierFunction}, 
    H::Union{Float64, Matrix{Float64}},
    F::Union{Float64, Vector{Float64}}
    )
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
            α = CBF.α(CBF.h(x))
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

    return CBFQuadProg(solve, H, F)
end

function CBFQuadProg(
    Σ::ControlAffineSystem, 
    CBFs::Vector{ControlBarrierFunction}, 
    H::Function,
    F::Union{Float64, Vector{Float64}}
    )
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
            α = CBF.α(CBF.h(x))
            @constraint(model, Lfh + Lgh*u >= -α)
        end
        @objective(model, Min, 0.5*u'*H(x)*u + F'*u)

        # Add control bounds on system - recall these default to unbounded controls
        if ~(Inf in Σ.b)
            @constraint(model, Σ.A * u .<= Σ.b)
        end

        # Solve QP
        optimize!(model)

        return Σ.m == 1 ? value(u) : value.(u)
    end

    return CBFQuadProg(solve, H, F)
end

function CBFQuadProg(
    Σ::ControlAffineSystem, 
    CBFs::Vector{ControlBarrierFunction}, 
    H::Union{Float64, Matrix{Float64}},
    F::Function
    )
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
            α = CBF.α(CBF.h(x))
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

    return CBFQuadProg(solve, H, F)
end

# If a single CBF passed in pack it into a vector and then call constructor
CBFQuadProg(Σ::ControlAffineSystem, CBF::ControlBarrierFunction) = CBFQuadProg(Σ, [CBF])

function CBFQuadProg(
    Σ::ControlAffineSystem,
    CBF::ControlBarrierFunction, 
    k::FeedbackController
    )
    return CBFQuadProg(Σ, [CBF], k)
end

function CBFQuadProg(
    Σ::ControlAffineSystem, 
    CBF::ControlBarrierFunction, 
    k::FeedbackController,
    H::Union{Float64, Matrix{Float64}}
    )
    return CBFQuadProg(Σ, [CBF], k, H)
end

function CBFQuadProg(
    Σ::ControlAffineSystem, 
    CBF::ControlBarrierFunction,
    k::FeedbackController,
    H::Function
    )
    return CBFQuadProg(Σ, [CBF], k, H)
end

function CBFQuadProg(
    Σ::ControlAffineSystem, 
    CBF::ControlBarrierFunction,
    H::Function,
    F::Function
    )
    return CBFQuadProg(Σ, [CBF], H, F)
end

function CBFQuadProg(
    Σ::ControlAffineSystem, 
    CBF::ControlBarrierFunction,
    H::Union{Float64, Matrix{Float64}},
    F::Union{Float64, Vector{Float64}}
    )
    return CBFQuadProg(Σ, [CBF], H, F)
end

function CBFQuadProg(
    Σ::ControlAffineSystem, 
    CBF::ControlBarrierFunction, 
    H::Function,
    F::Union{Float64, Vector{Float64}}
    )
    return CBFQuadProg(Σ, [CBF], H, F)
end

function CBFQuadProg(
    Σ::ControlAffineSystem, 
    CBF::ControlBarrierFunction,
    H::Union{Float64, Matrix{Float64}},
    F::Function
    )
    return CBFQuadProg(Σ, [CBF], H, F)
end