struct ISSCBFController <: FeedbackController
    h::Vector{ControlBarrierFunction}
    ε
    A
    b
end

ISSCBFController(h::ControlBarrierFunction, ε) = ISSCBFController([h], ε, nothing, nothing)
ISSCBFController(h::Vector{ControlBarrierFunction}, ε) = ISSCBFController(h, ε, nothing, nothing)

function iss_cbf_condition(h::ControlBarrierFunction, Σ::ControlAffineSystem, ε, x, u)
    r = h.relative_degree
    if r == 1
        Lfh = _Lfh(h, Σ, x)
        Lgh = _Lgh(h, Σ, x)
        φ = _φ(Σ, x)
        out = Lfh + Lgh*u + h.α(h(x)) - (1/ε)*norm(Lgh*φ)^2
    elseif r == 2
        φ = _φ(Σ, x)
        LgLfh = _LgLfh(h, Σ, x)
        out = _dψ1(h, Σ, x, u) + h.α(_ψ1(h, Σ, x)) - (1/ε)*norm(LgLfh*φ)^2
    end

    return out
end

function (k::ISSCBFController)(Σ::ControlAffineSystem, x)
    # Build QP and instantiate control decision variable
    model = Model(OSQP.Optimizer)
    set_silent(model)
    control_dim(Σ) == 1 ? @variable(model, u) : @variable(model, u[1:control_dim(Σ)])
    @objective(model, Min, 0.5*u'u)

    # Add CBF constraints
    for h in k.h
        @constraint(model, iss_cbf_condition(h, Σ, k.ε, x, u) >= 0.0)
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

function (k::ISSCBFController)(Σ::ControlAffineSystem, x, k0)
    # Build QP and instantiate control decision variable
    model = Model(OSQP.Optimizer)
    set_silent(model)
    control_dim(Σ) == 1 ? @variable(model, u) : @variable(model, u[1:control_dim(Σ)])
    @objective(model, Min, 0.5*u'u - k0'u)

    # Add CBF constraints
    for h in k.h
        @constraint(model, iss_cbf_condition(h, Σ, k.ε, x, u) >= 0.0)
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