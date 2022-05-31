struct ISSCBFCLFController <: FeedbackController
    h::Vector{ControlBarrierFunction}
    V::ControlLyapunovFunction
    εCBF
    εCLF
    A
    b
    p
end

function ISSCBFCLFController(h::ControlBarrierFunction, V::ControlLyapunovFunction, εCBF, εCLF)
    return ISSCBFCLFController([h], V, εCBF, εCLF, nothing, nothing, nothing)
end

function ISSCBFCLFController(h::ControlBarrierFunction, V::ControlLyapunovFunction, εCBF, εCLF, p)
    return ISSCBFCLFController([h], V, εCBF, εCLF, nothing, nothing, p)
end

function ISSCBFCLFController(h::ControlBarrierFunction, V::ControlLyapunovFunction, εCBF, εCLF, A, b)
    return ISSCBFCLFController([h], V, εCBF, εCLF, A, b, nothing)
end

function ISSCBFCLFController(h::Vector{ControlBarrierFunction}, V::ControlLyapunovFunction, εCBF, εCLF)
    return ISSCBFCLFController(h, V, εCBF, εCLF, nothing, nothing, nothing)
end

function ISSCBFCLFController(h::Vector{ControlBarrierFunction}, V::ControlLyapunovFunction, εCBF, εCLF, p)
    return ISSCBFCLFController(h, V, εCBF, εCLF, nothing, nothing, p)
end

function ISSCBFCLFController(h::Vector{ControlBarrierFunction}, V::ControlLyapunovFunction, εCBF, εCLF, A, b)
    return ISSCBFCLFController(h, V, εCBF, εCLF, A, b, nothing)
end

function (k::ISSCBFCLFController)(Σ::ControlAffineSystem, x)
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
    for (i, h) in enumerate(k.h)
        @constraint(model, iss_cbf_condition(h, Σ, k.εCBF[i], x, u) >= 0.0)
    end

    # Add (relaxed) CLF constraint
    if k.p === nothing
        @constraint(model, iss_clf_condition(k.V, Σ, k.εCLF, x, u) <= 0.0)
    else
        @constraint(model, iss_clf_condition(k.V, Σ, k.εCLF, x, u) <= δ)
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