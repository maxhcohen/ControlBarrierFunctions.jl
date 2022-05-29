struct CBFCLFController <: FeedbackController
    h::Vector{ControlBarrierFunction}
    V::ControlLyapunovFunction
    A
    b
    p
end

function CBFCLFController(h::ControlBarrierFunction, V::ControlLyapunovFunction)
    return CBFCLFController([h], V, nothing, nothing, nothing)
end

function CBFCLFController(h::ControlBarrierFunction, V::ControlLyapunovFunction, p) 
    return CBFCLFController([h], V, nothing, nothing, p)
end

function CBFCLFController(h::ControlBarrierFunction, V::ControlLyapunovFunction, A, b) 
    return CBFCLFController([h], V, A, b, nothing)
end

function CBFCLFController(h::Vector{ControlBarrierFunction}, V::ControlLyapunovFunction)
    return CBFCLFController(h, V, nothing, nothing, nothing)
end

function CBFCLFController(h::Vector{ControlBarrierFunction}, V::ControlLyapunovFunction, p) 
    return CBFCLFController(h, V, nothing, nothing, p)
end

function CBFCLFController(h::Vector{ControlBarrierFunction}, V::ControlLyapunovFunction, A, b) 
    return CBFCLFController(h, V, A, b, nothing)
end

function (k::CBFCLFController)(Σ::ControlAffineSystem, x)
    # Build QP and instantiate control decision variable
    model = Model(OSQP.Optimizer)
    set_silent(model)
    control_dim(Σ) == 1 ? @variable(model, u) : @variable(model, u[1:control_dim(Σ)])

    # Check if we want to add a relaxation variable to objective
    if k.p === nothing
        @objective(model, Min, 0.5*u'u)
    else
        @variable(model, δ)
        @objective(model, Min, 0.5*u'u + k.p*δ^2)
    end

    # Add CBF constraints
    for h in k.h
        α = h.α
        Lfh, Lgh = lie_derivatives(h, Σ, x)
        @constraint(model, Lfh + Lgh*u >= -α(h(x)))
    end

    # Add (relaxed) CLF constraint
    V = k.V
    γ = k.V.α
    LfV, LgV = lie_derivatives(V, Σ, x)
    if k.p === nothing
        @constraint(model, LfV + LgV*u <= -γ(V(x)))
    else
        @constraint(model, LfV + LgV*u <= -γ(V(x)) + δ)
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