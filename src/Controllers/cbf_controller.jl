struct CBFController <: FeedbackController
    h::ControlBarrierFunction
    A
    b
end

CBFController(h::ControlBarrierFunction) = CBFController(h, nothing, nothing)

function (k::CBFController)(Σ::ControlAffineSystem, x)
    # CBF stuff
    h = k.h
    α = k.h.α
    Lfh, Lgh = lie_derivatives(h, Σ, x)

    # Build QP and instantiate control decision variable
    model = Model(OSQP.Optimizer)
    set_silent(model)
    control_dim(Σ) == 1 ? @variable(model, u) : @variable(model, u[1:control_dim(Σ)])
    @objective(model, Min, 0.5*u'u)

    # Add CBF constraint
    @constraint(model, Lfh + Lgh*u >= -α(h(x)))

    # Check if we should add control constraints
    if k.A === nothing
    else
        @constraint(model, k.A * u .<= k.b)
    end

    # Solve QP
    optimize!(model)

    return control_dim(Σ) == 1 ? value(u) : value.(u)
end

function (k::CBFController)(Σ::ControlAffineSystem, x, k0)
    # CBF stuff
    h = k.h
    α = k.h.α
    Lfh, Lgh = lie_derivatives(h, Σ, x)

    # Build QP and instantiate control decision variable
    model = Model(OSQP.Optimizer)
    set_silent(model)
    control_dim(Σ) == 1 ? @variable(model, u) : @variable(model, u[1:control_dim(Σ)])
    @objective(model, Min, 0.5*u'u - k0'u)

    # Add CBF constraint
    @constraint(model, Lfh + Lgh*u >= -α(h(x)))

    # Check if we should add control constraints
    if k.A === nothing
    else
        @constraint(model, k.A * u .<= k.b)
    end

    # Solve QP
    optimize!(model)

    return control_dim(Σ) == 1 ? value(u) : value.(u)
end