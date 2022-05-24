"""
	TimeVaryingCLFQP <: Policy

Struct representing a Control Lyapunov Function based quadratic program (CLFQP).

# Fields
- `control`: function u = control(x, t) describing the control law.
"""
struct TimeVaryingCLFQP <: Policy
    control
end

"""
	(k::TimeVaryingCLFQP)(x, t)

Solve CLF-QP for control input at state x.
"""
function (κ::TimeVaryingCLFQP)(x, t)
    return κ.control(x, t)
end

"""
    TimeVaryingCLFQP(Σ::ControlAffineSystem, CLF::TimeVaryingCLF)
    TimeVaryingCLFQP(Σ::ControlAffineSystem, CLF::TimeVaryingCLF, A, b)
    TimeVaryingCLFQP(Σ::ControlAffineSystem, CLF::TimeVaryingCLF, A, b, p::Float64)

Construct a CLF-QP from a ControlAffineSystem and ControlLyapunovFunction.
"""
function TimeVaryingCLFQP(Σ::ControlAffineSystem, CLF::TimeVaryingCLF)
    function control(x, t)
        model = Model(OSQP.Optimizer)
        set_silent(model)
        Σ.m == 1 ? @variable(model, u) : @variable(model, u[1:(Σ.m)])
        @objective(model, Min, (1 / 2)u'u)
        @constraint(model, CLF.∇V(x, t) * (Σ.f(x) + Σ.g(x) * u) + CLF.∂ₜV(x, t) <= -CLF.α(CLF(x, t)))
        optimize!(model)

        return Σ.m == 1 ? value(u) : value.(u)
    end

    return TimeVaryingCLFQP(control)
end

function TimeVaryingCLFQP(Σ::ControlAffineSystem, CLF::TimeVaryingCLF, A, b)
    function control(x, t)
        model = Model(OSQP.Optimizer)
        set_silent(model)
        Σ.m == 1 ? @variable(model, u) : @variable(model, u[1:(Σ.m)])
        @objective(model, Min, (1 / 2)u'u)
        @constraint(model, CLF.∇V(x, t) * (Σ.f(x) + Σ.g(x) * u) + CLF.∂ₜV(x, t) <= -CLF.α(CLF(x, t)))
        @constraint(model, A * u .<= b)
        optimize!(model)

        return Σ.m == 1 ? value(u) : value.(u)
    end

    return TimeVaryingCLFQP(control)
end

function TimeVaryingCLFQP(Σ::ControlAffineSystem, CLF::TimeVaryingCLF, A, b, p::Float64)
    function control(x, t)
        model = Model(OSQP.Optimizer)
        set_silent(model)
        Σ.m == 1 ? @variable(model, u) : @variable(model, u[1:(Σ.m)])
        @variable(model, δ)
        @objective(model, Min, (1 / 2)u'u + p*δ^2)
        @constraint(model, CLF.∇V(x, t) * (Σ.f(x) + Σ.g(x) * u) + CLF.∂ₜV(x, t) <= -CLF.α(CLF(x, t)))
        @constraint(model, A * u .<= b)
        optimize!(model)

        return Σ.m == 1 ? value(u) : value.(u)
    end

    return TimeVaryingCLFQP(control)
end
