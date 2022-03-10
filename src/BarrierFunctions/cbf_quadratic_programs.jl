"""
	CBFQP <: Policy

Struct representing a Control Barrier Function based quadratic program (CBFQP).

# Fields
- `control`: function u = control(x) describing the control law.
"""
struct CBFQP <: Policy
    control
end

"""
	(κ::CBFQP)(x)

Solve CBF-QP at state x.
"""
function (κ::CBFQP)(x)
    return κ.control(x)
end

"""
	CBFQP(Σ::ControlAffineSystem, CBF::ControlBarrierFunction)
	CBFQP(Σ::ControlAffineSystem, CBF::ControlBarrierFunction, A, b)
    CBFQP(Σ::ControlAffineSystem, CBF::ControlBarrierFunction, κ::FeedbackPolicy)
	CBFQP(Σ::ControlAffineSystem, CBF::ControlBarrierFunction, κ::FeedbackPolicy, A, b)
	CBFQP(Σ::ControlAffineSystem, CBF::ControlBarrierFunction, κ::CLFQP)
	CBFQP(Σ::ControlAffineSystem, CBF::ControlBarrierFunction, κ::CLFQP, A, b)
	CBFQP(Σ::ControlAffineSystem, CBF::ControlBarrierFunction, CLF::ControlLyapunovFunction, p::Float64=10e2)
	CBFQP(Σ::ControlAffineSystem, CBF::ControlBarrierFunction, CLF::ControlLyapunovFunction, A, b, p::Float64=10e2)
	CBFQP(Σ::ControlAffineSystem, CBFs::Vector{ControlBarrierFunction})
	CBFQP(Σ::ControlAffineSystem, CBFs::Vector{ControlBarrierFunction}, A, b)
	CBFQP(Σ::ControlAffineSystem, CBFs::Vector{ControlBarrierFunction}, κ::FeedbackPolicy)
	CBFQP(Σ::ControlAffineSystem, CBFs::Vector{ControlBarrierFunction}, κ::FeedbackPolicy, A, b)
	CBFQP(Σ::ControlAffineSystem, CBFs::Vector{ControlBarrierFunction}, CLF::ControlLyapunovFunction, p::Float64=10e2)
	CBFQP(Σ::ControlAffineSystem, CBFs::Vector{ControlBarrierFunction}, CLF::ControlLyapunovFunction, A, b, p::Float64=10e2)

Construct a CBF-QP from a ControlAffineSystem and ControlBarrierFunction.
"""
function CBFQP(Σ::ControlAffineSystem, CBF::ControlBarrierFunction)
    function control(x)
        model = Model(OSQP.Optimizer)
        set_silent(model)
        Σ.m == 1 ? @variable(model, u) : @variable(model, u[1:(Σ.m)])
        @objective(model, Min, (1 / 2)u'u)
        @constraint(model, cbf, CBF.∇h(x) * (Σ.f(x) + Σ.g(x) * u) >= -CBF.α(CBF(x)))
        optimize!(model)

        return Σ.m == 1 ? value(u) : value.(u)
    end

    return CBFQP(control)
end

function CBFQP(Σ::ControlAffineSystem, CBF::ControlBarrierFunction, A, b)
    function control(x)
        model = Model(OSQP.Optimizer)
        set_silent(model)
        Σ.m == 1 ? @variable(model, u) : @variable(model, u[1:(Σ.m)])
        @objective(model, Min, (1 / 2)u'u)
        @constraint(model, cbf, CBF.∇h(x) * (Σ.f(x) + Σ.g(x) * u) >= -CBF.α(CBF(x)))
        @constraint(model, A * u .<= b)
        optimize!(model)

        return Σ.m == 1 ? value(u) : value.(u)
    end

    return CBFQP(control)
end

function CBFQP(Σ::ControlAffineSystem, CBF::ControlBarrierFunction, κ::FeedbackPolicy)
    function control(x)
        model = Model(OSQP.Optimizer)
        set_silent(model)
        Σ.m == 1 ? @variable(model, u) : @variable(model, u[1:(Σ.m)])
        @objective(model, Min, (1 / 2)u'u - κ(x)'u)
        @constraint(model, cbf, CBF.∇h(x) * (Σ.f(x) + Σ.g(x) * u) >= -CBF.α(CBF(x)))
        optimize!(model)

        return Σ.m == 1 ? value(u) : value.(u)
    end

    return CBFQP(control)
end

function CBFQP(Σ::ControlAffineSystem, CBF::ControlBarrierFunction, κ::FeedbackPolicy, A, b)
    function control(x)
        model = Model(OSQP.Optimizer)
        set_silent(model)
        Σ.m == 1 ? @variable(model, u) : @variable(model, u[1:(Σ.m)])
        @objective(model, Min, (1 / 2)u'u - κ(x)'u)
        @constraint(model, cbf, CBF.∇h(x) * (Σ.f(x) + Σ.g(x) * u) >= -CBF.α(CBF(x)))
        @constraint(model, A * u .<= b)
        optimize!(model)

        return Σ.m == 1 ? value(u) : value.(u)
    end

    return CBFQP(control)
end

function CBFQP(Σ::ControlAffineSystem, CBF::ControlBarrierFunction, κ::CLFQP)
    function control(x)
        model = Model(OSQP.Optimizer)
        set_silent(model)
        Σ.m == 1 ? @variable(model, u) : @variable(model, u[1:(Σ.m)])
        @objective(model, Min, (1 / 2)u'u - κ(x)'u)
        @constraint(model, cbf, CBF.∇h(x) * (Σ.f(x) + Σ.g(x) * u) >= -CBF.α(CBF(x)))
        optimize!(model)

        return Σ.m == 1 ? value(u) : value.(u)
    end

    return CBFQP(control)
end

function CBFQP(Σ::ControlAffineSystem, CBF::ControlBarrierFunction, κ::CLFQP, A, b)
    function control(x)
        model = Model(OSQP.Optimizer)
        set_silent(model)
        Σ.m == 1 ? @variable(model, u) : @variable(model, u[1:(Σ.m)])
        @objective(model, Min, (1 / 2)u'u - κ(x)'u)
        @constraint(model, cbf, CBF.∇h(x) * (Σ.f(x) + Σ.g(x) * u) >= -CBF.α(CBF(x)))
        @constraint(model, A * u .<= b)
        optimize!(model)

        return Σ.m == 1 ? value(u) : value.(u)
    end

    return CBFQP(control)
end

function CBFQP(
    Σ::ControlAffineSystem,
    CBF::ControlBarrierFunction,
    CLF::ControlLyapunovFunction,
    p::Float64=10e2,
)
    function control(x)
        model = Model(OSQP.Optimizer)
        set_silent(model)
        Σ.m == 1 ? @variable(model, u) : @variable(model, u[1:(Σ.m)])
        @variable(model, δ)
        @objective(model, Min, (1 / 2)u'u + p * δ^2)
        @constraint(model, cbf, CBF.∇h(x) * (Σ.f(x) + Σ.g(x) * u) >= -CBF.α(CBF(x)))
        @constraint(model, clf, CLF.∇V(x) * (Σ.f(x) + Σ.g(x) * u) <= -CLF.α(CLF(x)) + δ)
        optimize!(model)

        return Σ.m == 1 ? value(u) : value.(u)
    end

    return CBFQP(control)
end

function CBFQP(
    Σ::ControlAffineSystem,
    CBF::ControlBarrierFunction,
    CLF::ControlLyapunovFunction,
    A,
    b,
    p::Float64=10e2,
)
    function control(x)
        model = Model(OSQP.Optimizer)
        set_silent(model)
        Σ.m == 1 ? @variable(model, u) : @variable(model, u[1:(Σ.m)])
        @variable(model, δ)
        @objective(model, Min, (1 / 2)u'u + p * δ^2)
        @constraint(model, cbf, CBF.∇h(x) * (Σ.f(x) + Σ.g(x) * u) >= -CBF.α(CBF(x)))
        @constraint(model, clf, CLF.∇V(x) * (Σ.f(x) + Σ.g(x) * u) <= -CLF.α(CLF(x)) + δ)
        @constraint(model, A * u .<= b)
        optimize!(model)

        return Σ.m == 1 ? value(u) : value.(u)
    end

    return CBFQP(control)
end

function CBFQP(Σ::ControlAffineSystem, CBFs::Vector{ControlBarrierFunction})
    function control(x)
        model = Model(OSQP.Optimizer)
        set_silent(model)
        Σ.m == 1 ? @variable(model, u) : @variable(model, u[1:(Σ.m)])
        @objective(model, Min, (1 / 2)u'u)
        for CBF in CBFs
            @constraint(model, CBF.∇h(x) * (Σ.f(x) + Σ.g(x) * u) >= -CBF.α(CBF(x)))
        end
        optimize!(model)

        return Σ.m == 1 ? value(u) : value.(u)
    end

    return CBFQP(control)
end

function CBFQP(Σ::ControlAffineSystem, CBFs::Vector{ControlBarrierFunction}, A, b)
    function control(x)
        model = Model(OSQP.Optimizer)
        set_silent(model)
        Σ.m == 1 ? @variable(model, u) : @variable(model, u[1:(Σ.m)])
        @objective(model, Min, (1 / 2)u'u)
        for CBF in CBFs
            @constraint(model, CBF.∇h(x) * (Σ.f(x) + Σ.g(x) * u) >= -CBF.α(CBF(x)))
        end
        @constraint(model, A * u .<= b)
        optimize!(model)

        return Σ.m == 1 ? value(u) : value.(u)
    end

    return CBFQP(control)
end

function CBFQP(
    Σ::ControlAffineSystem, CBFs::Vector{ControlBarrierFunction}, κ::FeedbackPolicy
)
    function control(x)
        model = Model(OSQP.Optimizer)
        set_silent(model)
        Σ.m == 1 ? @variable(model, u) : @variable(model, u[1:(Σ.m)])
        @objective(model, Min, (1 / 2)u'u - κ(x)'u)
        for CBF in CBFs
            @constraint(model, CBF.∇h(x) * (Σ.f(x) + Σ.g(x) * u) >= -CBF.α(CBF(x)))
        end
        optimize!(model)

        return Σ.m == 1 ? value(u) : value.(u)
    end

    return CBFQP(control)
end

function CBFQP(
    Σ::ControlAffineSystem, CBFs::Vector{ControlBarrierFunction}, κ::FeedbackPolicy, A, b
)
    function control(x)
        model = Model(OSQP.Optimizer)
        set_silent(model)
        Σ.m == 1 ? @variable(model, u) : @variable(model, u[1:(Σ.m)])
        @objective(model, Min, (1 / 2)u'u - κ(x)'u)
        for CBF in CBFs
            @constraint(model, CBF.∇h(x) * (Σ.f(x) + Σ.g(x) * u) >= -CBF.α(CBF(x)))
        end
        @constraint(model, A * u .<= b)
        optimize!(model)

        return Σ.m == 1 ? value(u) : value.(u)
    end

    return CBFQP(control)
end

function CBFQP(
    Σ::ControlAffineSystem,
    CBFs::Vector{ControlBarrierFunction},
    CLF::ControlLyapunovFunction,
    p::Float64=10e2,
)
    function control(x)
        model = Model(OSQP.Optimizer)
        set_silent(model)
        Σ.m == 1 ? @variable(model, u) : @variable(model, u[1:(Σ.m)])
        @variable(model, δ)
        @objective(model, Min, (1 / 2)u'u + p * δ^2)
        for CBF in CBFs
            @constraint(model, CBF.∇h(x) * (Σ.f(x) + Σ.g(x) * u) >= -CBF.α(CBF(x)))
        end
        @constraint(model, clf, CLF.∇V(x) * (Σ.f(x) + Σ.g(x) * u) <= -CLF.α(CLF(x)) + δ)
        optimize!(model)

        return Σ.m == 1 ? value(u) : value.(u)
    end

    return CBFQP(control)
end

function CBFQP(
    Σ::ControlAffineSystem,
    CBFs::Vector{ControlBarrierFunction},
    CLF::ControlLyapunovFunction,
    A,
    b,
    p::Float64=10e2,
)
    function control(x)
        model = Model(OSQP.Optimizer)
        set_silent(model)
        Σ.m == 1 ? @variable(model, u) : @variable(model, u[1:(Σ.m)])
        @variable(model, δ)
        @objective(model, Min, (1 / 2)u'u + p * δ^2)
        for CBF in CBFs
            @constraint(model, CBF.∇h(x) * (Σ.f(x) + Σ.g(x) * u) >= -CBF.α(CBF(x)))
        end
        @constraint(model, clf, CLF.∇V(x) * (Σ.f(x) + Σ.g(x) * u) <= -CLF.α(CLF(x)) + δ)
        @constraint(model, A * u .<= b)
        optimize!(model)

        return Σ.m == 1 ? value(u) : value.(u)
    end

    return CBFQP(control)
end
