struct CLFController <: FeedbackController
    V::ControlLyapunovFunction
    A
    b
    p
end

CLFController(V::ControlLyapunovFunction) = CLFController(V, nothing, nothing, nothing)
CLFController(V::ControlLyapunovFunction, p) = CLFController(V, nothing, nothing, p)
CLFController(V::ControlLyapunovFunction, A, b) = CLFController(V, A, b, nothing)

function (k::CLFController)(Σ::ControlAffineSystem, x)
    # CLF stuff
    V = k.V
    α = k.V.α
    LfV, LgV = lie_derivatives(k.V, Σ, x)

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
    @constraint(model, LfV + LgV*u <= -α(V(x)))

    # Check if we should add control constraints
    if k.A === nothing
    else
        @constraint(model, k.A * u .<= k.b)
    end

    # Solve QP
    optimize!(model)

    return control_dim(Σ) == 1 ? value(u) : value.(u)
end