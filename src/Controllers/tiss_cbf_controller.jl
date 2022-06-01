struct TISSCBFController <: FeedbackController
    h::Vector{ControlBarrierFunction}
    ε
    λ
    A
    b
end

TISSCBFController(h::ControlBarrierFunction, ε, λ) = TISSCBFController([h], ε, λ, nothing, nothing)
TISSCBFController(h::Vector{ControlBarrierFunction}, ε, λ) = TISSCBFController(h, ε, λ, nothing, nothing)

function tiss_cbf_condition(h::ControlBarrierFunction, Σ::ControlAffineSystem, ε0, λ, x, u)
    r = h.relative_degree
    if r == 1
        Lfh = _Lf0h(h, Σ, x)
        Lgh = _Lgh(h, Σ, x)
        φ = _φ(Σ, x)
        ε = ε0*exp(λ * h(x))
        out = Lfh + Lgh*u + h.α(h(x)) - (1/ε)*norm(Lgh*φ)^2
    elseif r == 2
        φ = _φ(Σ, x)
        LgLfh = _LgLf0h(h, Σ, x)
        ε = ε0*exp(λ * _ψ1(h, Σ, x))
        out = _dψ10(h, Σ, x, u) + h.α(_ψ10(h, Σ, x)) - (1/ε)*norm(LgLfh*φ)^2
    end

    return out
end

function (k::TISSCBFController)(Σ::ControlAffineSystem, x)
    # Build QP and instantiate control decision variable
    model = Model(OSQP.Optimizer)
    set_silent(model)
    control_dim(Σ) == 1 ? @variable(model, u) : @variable(model, u[1:control_dim(Σ)])
    @objective(model, Min, 0.5*u'u)

    # Add CBF constraints
    for h in k.h
        @constraint(model, tiss_cbf_condition(h, Σ, k.ε, k.λ, x, u) >= 0.0)
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

function (k::TISSCBFController)(Σ::ControlAffineSystem, x, k0)
    # Build QP and instantiate control decision variable
    model = Model(OSQP.Optimizer)
    set_silent(model)
    control_dim(Σ) == 1 ? @variable(model, u) : @variable(model, u[1:control_dim(Σ)])
    @objective(model, Min, 0.5*u'u - k0'u)

    # Add CBF constraints
    for h in k.h
        @constraint(model, tiss_cbf_condition(h, Σ, k.ε, k.λ, x, u) >= 0.0)
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