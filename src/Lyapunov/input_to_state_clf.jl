"""
    InputToStateCLF <: CertificateFunction

# Fields
- `CLF::ControlLyapunovFunction` : nominal CLF for a `ControlAffineSystem`
- `ε::Float64` : nonlinear damping coefficient for ISS term
"""
struct InputToStateCLF <: CertificateFunction
    CLF::ControlLyapunovFunction
    ε::Float64
end

InputToStateCLF(CLF::ControlLyapunovFunction) = InputToStateCLF(CLF, 1.0)

function CLFQuadProg(Σ::ControlAffineSystem, ISSCLF::InputToStateCLF)
    # Set parameters for objective function
    H = Σ.m == 1 ? 1.0 : Matrix(1.0I, Σ.m, Σ.m)
    F = Σ.m == 1 ? 0.0 : zeros(Σ.m)

    # Pull out underlying CLF
    CLF = ISSCLF.CLF
    ε = ISSCLF.ε

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

    return CLFQuadProg(solve, H, F)
end

function CLFQuadProg(
    Σ::ControlAffineSystem,
    ISSCLF::InputToStateCLF,
    H::Union{Float64, Matrix{Float64}},
    F::Union{Float64, Vector{Float64}}
    )

    # Pull out underlying CLF
    CLF = ISSCLF.CLF
    ε = ISSCLF.ε

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

    return CLFQuadProg(solve, H, F)
end

function CLFQuadProg(
    Σ::ControlAffineSystem,
    ISSCLF::InputToStateCLF,
    H::Function,
    F::Function
    )

    # Pull out underlying CLF
    CLF = ISSCLF.CLF
    ε = ISSCLF.ε

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
            @objective(model, Min, 0.5*u'*H(x)*u + F(x)'*u + CLF.p*δ^2)
        else
            @constraint(model, LfV + LgV*u <= -γ)
            @objective(model, Min, 0.5*u'*H(x)*u + F(x)'*u)
        end

        # Add control bounds on system - recall these default to unbounded controls
        if ~(Inf in Σ.b)
            @constraint(model, Σ.A * u .<= Σ.b)
        end

        # Solve QP
        optimize!(model)

        return Σ.m == 1 ? value(u) : value.(u)
    end

    return CLFQuadProg(solve, H, F)
end

function CLFQuadProg(
    Σ::ControlAffineSystem,
    ISSCLF::InputToStateCLF,
    H::Function,
    F::Union{Float64, Vector{Float64}}
    )
    # Pull out underlying CLF
    CLF = ISSCLF.CLF
    ε = ISSCLF.ε

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
            @objective(model, Min, 0.5*u'*H(x)*u + F'*u + CLF.p*δ^2)
        else
            @constraint(model, LfV + LgV*u <= -γ)
            @objective(model, Min, 0.5*u'*H(x)*u + F'*u)
        end

        # Add control bounds on system - recall these default to unbounded controls
        if ~(Inf in Σ.b)
            @constraint(model, Σ.A * u .<= Σ.b)
        end

        # Solve QP
        optimize!(model)

        return Σ.m == 1 ? value(u) : value.(u)
    end

    return CLFQuadProg(solve, H, F)
end

function CLFQuadProg(
    Σ::ControlAffineSystem,
    ISSCLF::InputToStateCLF,
    H::Union{Float64, Matrix{Float64}},
    F::Function
    )

    # Pull out underlying CLF
    CLF = ISSCLF.CLF
    ε = ISSCLF.ε

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
            @objective(model, Min, 0.5*u'*H*u + F(x)'*u + CLF.p*δ^2)
        else
            @constraint(model, LfV + LgV*u <= -γ)
            @objective(model, Min, 0.5*u'*H*u + F(x)'*u)
        end

        # Add control bounds on system - recall these default to unbounded controls
        if ~(Inf in Σ.b)
            @constraint(model, Σ.A * u .<= Σ.b)
        end

        # Solve QP
        optimize!(model)

        return Σ.m == 1 ? value(u) : value.(u)
    end

    return CLFQuadProg(solve, H, F)
end

"Compute gradient of a ISSCLF evaluated at `x`."
gradient(ISSCLF::InputToStateCLF, x) = gradient(ISSCLF.CLF, x)

"Compute lie derivatives of ISSCLF at `x`."
function lie_derivatives(ISSCLF::InputToStateCLF, Σ::ControlAffineSystem, x)
    return lie_derivatives(ISSCLF.CLF, Σ, x)
end

"Compute Lie derivative of ISSCLF along drift dynamics ``LfV(x) = ∇V(x) * f(x)``."
function drift_lie_derivative(ISSCLF::InputToStateCLF, Σ::ControlAffineSystem, x)
    return gradient(ISSCLF, x) * Σ.f(x)
end

"Compute Lie derivative of ISSCLF along control directions ``LgV(x) = ∇V(x) * g(x)``."
function control_lie_derivative(ISSCLF::InputToStateCLF, Σ::ControlAffineSystem, x)
    return gradient(ISSCLF, x) * Σ.g(x)
end