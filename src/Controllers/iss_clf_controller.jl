struct ISSCLFController <: FeedbackController
    V::ControlLyapunovFunction
    ε
    A
    b
    p
end

ISSCLFController(V::ControlLyapunovFunction, ε) = ISSCLFController(V, ε, nothing, nothing, nothing)
ISSCLFController(V::ControlLyapunovFunction, ε, p) = ISSCLFController(V, ε, nothing, nothing, p)
ISSCLFController(V::ControlLyapunovFunction, ε, A, b) = ISSCLFController(V, ε, A, b, nothing)

function iss_clf_condition(V::ControlLyapunovFunction, Σ::ControlAffineSystem, ε, x, u)
    LfV = _Lf0V(V, Σ, x)
    LgV = _LgV(V, Σ, x)
    φ = _φ(Σ, x)

    return LfV + LgV*u + (1/ε)*norm(LgV*φ)^2  + V.α(V(x))
end

function (k::ISSCLFController)(Σ::ControlAffineSystem, x)
    # Build QP and instantiate control decision variable
    model = Model(OSQP.Optimizer)
    set_silent(model)
    control_dim(Σ) == 1 ? @variable(model, u) : @variable(model, u[1:control_dim(Σ)])

    # Check if we want to add a relaxation variable
    if k.p === nothing
        @objective(model, Min, 0.5*u'u)
    else
        @variable(model, δ)
        @objective(model, Min, 0.5*u'u + k.p*δ^2)
    end

    # Add CLF constraint
    if k.p === nothing
        @constraint(model, iss_clf_condition(k.V, Σ, k.ε, x, u) <= 0.0)
    else
        @constraint(model, iss_clf_condition(k.V, Σ, k.ε, x, u) <= δ)
    end

    # Check if we should add control constraints
    if k.A === nothing
    else
        @constraint(model, k.A * u .<= k.b)
    end

    # Solve QP
    optimize!(model)

    return control_dim(Σ) == 1 ? value(u) : value.(u)
end