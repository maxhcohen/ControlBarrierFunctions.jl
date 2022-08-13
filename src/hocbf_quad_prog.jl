function CBFQuadProg(Σ::ControlAffineSystem, HOCBFs::Vector{SecondOrderCBF})
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
            α = HOCBF.α2(HOCBF.ψ1(x))
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

    return CBFQuadProg(solve, H, F)
end

function CBFQuadProg(
    Σ::ControlAffineSystem, 
    HOCBFs::Vector{SecondOrderCBF}, 
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

        # Set HOCBF constraint and objective
        for HOCBF in HOCBFs
            Lfψ = drift_lie_derivative(HOCBF, Σ, x)
            Lgψ = control_lie_derivative(HOCBF, Σ, x)
            α = HOCBF.α2(HOCBF.ψ1(x))
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

    return CBFQuadProg(solve, H, F)
end

function CBFQuadProg(
    Σ::ControlAffineSystem, 
    HOCBFs::Vector{SecondOrderCBF}, 
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

        # Set HOCBF constraint and objective
        for HOCBF in HOCBFs
            Lfψ = drift_lie_derivative(HOCBF, Σ, x)
            Lgψ = control_lie_derivative(HOCBF, Σ, x)
            α = HOCBF.α2(HOCBF.ψ1(x))
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

    return CBFQuadProg(solve, H, F)
end

function CBFQuadProg(
    Σ::ControlAffineSystem, 
    HOCBFs::Vector{SecondOrderCBF}, 
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

        # Set HOCBF constraint and objective
        for HOCBF in HOCBFs
            Lfψ = drift_lie_derivative(HOCBF, Σ, x)
            Lgψ = control_lie_derivative(HOCBF, Σ, x)
            α = HOCBF.α2(HOCBF.ψ1(x))
            @constraint(model, Lfψ + Lgψ*u >= -α)
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
    HOCBFs::Vector{SecondOrderCBF}, 
    H::Function,
    F::Function
    )
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
            α = HOCBF.α2(HOCBF.ψ1(x))
            @constraint(model, Lfψ + Lgψ*u >= -α)
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
    HOCBFs::Vector{SecondOrderCBF}, 
    H::Union{Float64, Matrix{Float64}},
    F::Union{Float64, Vector{Float64}}
    )
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
            α = HOCBF.α2(HOCBF.ψ1(x))
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

    return CBFQuadProg(solve, H, F)
end

function CBFQuadProg(
    Σ::ControlAffineSystem, 
    HOCBFs::Vector{SecondOrderCBF}, 
    H::Function,
    F::Union{Float64, Vector{Float64}}
    )
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
            α = HOCBF.α2(HOCBF.ψ1(x))
            @constraint(model, Lfψ + Lgψ*u >= -α)
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
    HOCBFs::Vector{SecondOrderCBF},
    H::Union{Float64, Matrix{Float64}},
    F::Function
    )
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
            α = HOCBF.α2(HOCBF.ψ1(x))
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

    return CBFQuadProg(solve, H, F)
end

# If a single HOCBF passed in pack it into a vector and then call constructor
CBFQuadProg(Σ::ControlAffineSystem, HOCBF::SecondOrderCBF) = CBFQuadProg(Σ, [HOCBF])

function CBFQuadProg(
    Σ::ControlAffineSystem,
    HOCBF::SecondOrderCBF, 
    k::FeedbackController
    )
    return CBFQuadProg(Σ, [HOCBF], k)
end

function CBFQuadProg(
    Σ::ControlAffineSystem, 
    HOCBF::SecondOrderCBF, 
    k::FeedbackController,
    H::Union{Float64, Matrix{Float64}}
    )
    return CBFQuadProg(Σ, [HOCBF], k, H)
end

function CBFQuadProg(
    Σ::ControlAffineSystem, 
    HOCBF::SecondOrderCBF, 
    k::FeedbackController,
    H::Function
    )
    return CBFQuadProg(Σ, [HOCBF], k, H)
end

function CBFQuadProg(
    Σ::ControlAffineSystem, 
    HOCBF::SecondOrderCBF, 
    H::Function,
    F::Function
    )
    return CBFQuadProg(Σ, [HOCBF], H, F)
end

function CBFQuadProg(
    Σ::ControlAffineSystem, 
    HOCBF::SecondOrderCBF, 
    H::Union{Float64, Matrix{Float64}},
    F::Union{Float64, Vector{Float64}}
    )
    return CBFQuadProg(Σ, [HOCBF], H, F)
end

function CBFQuadProg(
    Σ::ControlAffineSystem, 
    HOCBF::SecondOrderCBF, 
    H::Function,
    F::Union{Float64, Vector{Float64}}
    )
    return CBFQuadProg(Σ, [HOCBF], H, F)
end

function CBFQuadProg(
    Σ::ControlAffineSystem, 
    HOCBF::SecondOrderCBF, 
    H::Union{Float64, Matrix{Float64}},
    F::Function
    )
    return CBFQuadProg(Σ, [HOCBF], H, F)
end