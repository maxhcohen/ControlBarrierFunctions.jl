struct CBFController <: FeedbackController
    h::Vector{ControlBarrierFunction}
    A
    b
end

CBFController(h::ControlBarrierFunction) = CBFController([h], nothing, nothing)
CBFController(h::Vector{ControlBarrierFunction}) = CBFController(h, nothing, nothing)

function (k::CBFController)(Σ::ControlAffineSystem, x)
    # Build QP and instantiate control decision variable
    model = Model(OSQP.Optimizer)
    set_silent(model)
    control_dim(Σ) == 1 ? @variable(model, u) : @variable(model, u[1:control_dim(Σ)])
    @objective(model, Min, 0.5*u'u)

    # Add CBF constraints
    for h in k.h
        @constraint(model, cbf_condition(h, Σ, x, u) >= 0.0)
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

function (k::CBFController)(Σ::ControlAffineSystem, x, k0)
    # Build QP and instantiate control decision variable
    model = Model(OSQP.Optimizer)
    set_silent(model)
    control_dim(Σ) == 1 ? @variable(model, u) : @variable(model, u[1:control_dim(Σ)])
    @objective(model, Min, 0.5*u'u - k0'u)

    # Add CBF constraints
    for h in k.h
        @constraint(model, cbf_condition(h, Σ, x, u) >= 0.0)
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