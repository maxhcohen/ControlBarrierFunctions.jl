"""
	TimeVaryingCBFQP <: Policy

Struct representing a Control Barrier Function based quadratic program (CBFQP).

# Fields
- `control`: function u = control(x) describing the control law.
"""
struct TimeVaryingCBFQP <: Policy
    control
end

"""
	(κ::TimeVaryingCBFQP)(x, t)

Solve CBF-QP at state x and time t.
"""
function (κ::TimeVaryingCBFQP)(x, t)
    return κ.control(x, t)
end

"""
    TimeVaryingCBFQP(Σ::ControlAffineSystem, CBF::ControlBarrierFunction, κ::TimeVaryingCLFQP)
    TimeVaryingCBFQP(Σ::ControlAffineSystem, CBF::ControlBarrierFunction, κ::TimeVaryingCLFQP, A, b)
    TimeVaryingCBFQP(Σ::ControlAffineSystem, HOCBF::HighOrderCBF, κ::TimeVaryingCLFQP)
    TimeVaryingCBFQP(Σ::ControlAffineSystem, HOCBF::HighOrderCBF, κ::TimeVaryingCLFQP, A, b)
    TimeVaryingCBFQP(Σ::ControlAffineSystem, CBFs::Vector{ControlBarrierFunction}, κ::TimeVaryingCLFQP)
    TimeVaryingCBFQP(Σ::ControlAffineSystem, CBFs::Vector{ControlBarrierFunction}, κ::TimeVaryingCLFQP, A, b)
    TimeVaryingCBFQP(Σ::ControlAffineSystem, HOCBFs::Vector{HighOrderCBF}, κ::TimeVaryingCLFQP)
    TimeVaryingCBFQP(Σ::ControlAffineSystem, HOCBFs::Vector{HighOrderCBF}, κ::TimeVaryingCLFQP, A, b)

Construct a TimeVaryingCBFQP.
"""
function TimeVaryingCBFQP(Σ::ControlAffineSystem, CBF::ControlBarrierFunction, κ::TimeVaryingCLFQP)
    function control(x, t)
        model = Model(OSQP.Optimizer)
        set_silent(model)
        Σ.m == 1 ? @variable(model, u) : @variable(model, u[1:(Σ.m)])
        @objective(model, Min, (1 / 2)u'u - κ(x,t)'u)
        @constraint(model, CBF.∇h(x) * (Σ.f(x) + Σ.g(x) * u) >= -CBF.α(CBF(x)))
        optimize!(model)

        return Σ.m == 1 ? value(u) : value.(u)
    end

    return TimeVaryingCBFQP(control)
end

function TimeVaryingCBFQP(Σ::ControlAffineSystem, CBF::ControlBarrierFunction, κ::TimeVaryingCLFQP, A, b)
    function control(x, t)
        model = Model(OSQP.Optimizer)
        set_silent(model)
        Σ.m == 1 ? @variable(model, u) : @variable(model, u[1:(Σ.m)])
        @objective(model, Min, (1 / 2)u'u - κ(x,t)'u)
        @constraint(model, CBF.∇h(x) * (Σ.f(x) + Σ.g(x) * u) >= -CBF.α(CBF(x)))
        @constraint(model, A * u .<= b)
        optimize!(model)

        return Σ.m == 1 ? value(u) : value.(u)
    end

    return TimeVaryingCBFQP(control)
end

function TimeVaryingCBFQP(Σ::ControlAffineSystem, HOCBF::HighOrderCBF, κ::TimeVaryingCLFQP)
    function control(x, t)
        model = Model(OSQP.Optimizer)
        set_silent(model)
        Σ.m == 1 ? @variable(model, u) : @variable(model, u[1:(Σ.m)])
        @objective(model, Min, (1 / 2)u'u - κ(x,t)'u)
        @constraint(model, HOCBF.∇ψ(x) * (Σ.f(x) + Σ.g(x) * u) >= -HOCBF.α(HOCBF(x)))
        optimize!(model)

        return Σ.m == 1 ? value(u) : value.(u)
    end

    return TimeVaryingCBFQP(control)
end

function TimeVaryingCBFQP(Σ::ControlAffineSystem, HOCBF::HighOrderCBF, κ::TimeVaryingCLFQP, A, b)
    function control(x, t)
        model = Model(OSQP.Optimizer)
        set_silent(model)
        Σ.m == 1 ? @variable(model, u) : @variable(model, u[1:(Σ.m)])
        @objective(model, Min, (1 / 2)u'u - κ(x,t)'u)
        @constraint(model, HOCBF.∇ψ(x) * (Σ.f(x) + Σ.g(x) * u) >= -HOCBF.α(HOCBF(x)))
        @constraint(model, A * u .<= b)
        optimize!(model)

        return Σ.m == 1 ? value(u) : value.(u)
    end

    return TimeVaryingCBFQP(control)
end

function TimeVaryingCBFQP(Σ::ControlAffineSystem, CBFs::Vector{ControlBarrierFunction}, κ::TimeVaryingCLFQP)
    function control(x, t)
        model = Model(OSQP.Optimizer)
        set_silent(model)
        Σ.m == 1 ? @variable(model, u) : @variable(model, u[1:(Σ.m)])
        @objective(model, Min, (1 / 2)u'u - κ(x,t)'u)
        for CBF in CBFs
            @constraint(model, CBF.∇h(x) * (Σ.f(x) + Σ.g(x) * u) >= -CBF.α(CBF(x)))
        end
        optimize!(model)

        return Σ.m == 1 ? value(u) : value.(u)
    end

    return TimeVaryingCBFQP(control)
end

function TimeVaryingCBFQP(Σ::ControlAffineSystem, CBFs::Vector{ControlBarrierFunction}, κ::TimeVaryingCLFQP, A, b)
    function control(x, t)
        model = Model(OSQP.Optimizer)
        set_silent(model)
        Σ.m == 1 ? @variable(model, u) : @variable(model, u[1:(Σ.m)])
        @objective(model, Min, (1 / 2)u'u - κ(x,t)'u)
        for CBF in CBFs
            @constraint(model, CBF.∇h(x) * (Σ.f(x) + Σ.g(x) * u) >= -CBF.α(CBF(x)))
        end
        @constraint(model, A * u .<= b)
        optimize!(model)

        return Σ.m == 1 ? value(u) : value.(u)
    end

    return TimeVaryingCBFQP(control)
end

function TimeVaryingCBFQP(Σ::ControlAffineSystem, HOCBFs::Vector{HighOrderCBF}, κ::TimeVaryingCLFQP)
    function control(x, t)
        model = Model(OSQP.Optimizer)
        set_silent(model)
        Σ.m == 1 ? @variable(model, u) : @variable(model, u[1:(Σ.m)])
        @objective(model, Min, (1 / 2)u'u - κ(x,t)'u)
        for HOCBF in HOCBFs 
            @constraint(model, HOCBF.∇ψ(x) * (Σ.f(x) + Σ.g(x) * u) >= -HOCBF.α(HOCBF(x)))
        end
        optimize!(model)

        return Σ.m == 1 ? value(u) : value.(u)
    end

    return TimeVaryingCBFQP(control)
end

function TimeVaryingCBFQP(Σ::ControlAffineSystem, HOCBFs::Vector{HighOrderCBF}, κ::TimeVaryingCLFQP, A, b)
    function control(x, t)
        model = Model(OSQP.Optimizer)
        set_silent(model)
        Σ.m == 1 ? @variable(model, u) : @variable(model, u[1:(Σ.m)])
        @objective(model, Min, (1 / 2)u'u - κ(x,t)'u)
        for HOCBF in HOCBFs 
            @constraint(model, HOCBF.∇ψ(x) * (Σ.f(x) + Σ.g(x) * u) >= -HOCBF.α(HOCBF(x)))
        end
        @constraint(model, A * u .<= b)
        optimize!(model)

        return Σ.m == 1 ? value(u) : value.(u)
    end

    return TimeVaryingCBFQP(control)
end