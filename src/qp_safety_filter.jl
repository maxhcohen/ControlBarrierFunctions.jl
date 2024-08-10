"""
    QPSafetyFilter <: SafetyFilter

Controller that solves a control barrier function-based quadratic program (CBF-QP).

Uses OSQP to solve the corresponding QP.

# Fields
- `k::Function` : function that computes safe control actions
"""
struct QPSafetyFilter <: SafetyFilter
    k::Function
end

"""
    (k::QPSafetyFilter)(x)

Functors for evaluating QP-based safety filter
"""
(k::QPSafetyFilter)(x) = k.k(x)
(k::QPSafetyFilter)(x, t) = k.k(x, t)

"""
    QPSafetyFilter(cbfs::Vector{ControlBarrierFunction}, Σ::ControlAffineSystem, kd::Function)

Construct an QPSafetyFilter from a cbf and a desired controller.
"""
function QPSafetyFilter(
    cbfs::Vector{ControlBarrierFunction}, Σ::ControlAffineSystem, kd::Function
)
    try
        kd(Σ.n == 1 ? rand() : rand(Σ.n), 0.0) # See if desired controller is time-varying
    catch e
        if isa(e, MethodError) # If controller is not time-varying 
            return QPSafetyFilter(x -> solve_cbf_qp(x, Σ, cbfs, kd))
        else
            return e
        end
    else # If controller is time-varying
        return QPSafetyFilter((x, t) -> solve_time_varying_cbf_qp(x, t, Σ, cbfs, kd))
    end
end

"""
    QPSafetyFilter(
    cbfs::Vector{ControlBarrierFunction}, Σ::ControlAffineSystem, kd::Function, umin, umax
)

Construct an QPSafetyFilter from a cbf and a desired controller.
"""
function QPSafetyFilter(
    cbfs::Vector{ControlBarrierFunction}, Σ::ControlAffineSystem, kd::Function, umin, umax
)
    try
        kd(Σ.n == 1 ? rand() : rand(Σ.n), 0.0) # See if desired controller is time-varying
    catch e
        if isa(e, MethodError) # If controller is not time-varying 
            return QPSafetyFilter(x -> solve_cbf_qp(x, Σ, cbfs, kd, umin, umax))
        else
            return e
        end
    else # If controller is time-varying
        return QPSafetyFilter((x, t) -> solve_time_varying_cbf_qp(x, t, Σ, cbfs, kd))
    end
end

"""
    QPSafetyFilter(cbf::ControlBarrierFunction, Σ::ControlAffineSystem, kd::Function)
    QPSafetyFilter(cbf::ControlBarrierFunction, Σ::ControlAffineSystem, kd::Function, umin, umax)

Add ability to pass in single CBF.
"""
function QPSafetyFilter(cbf::ControlBarrierFunction, Σ::ControlAffineSystem, kd::Function)
    return QPSafetyFilter([cbf], Σ, kd)
end
function QPSafetyFilter(
    cbf::ControlBarrierFunction, Σ::ControlAffineSystem, kd::Function, umin, umax
)
    return QPSafetyFilter([cbf], Σ, kd, umin, umax)
end

"""
    solve_cbf_qp(x, Σ::ControlAffineSystem, cbfs::Vector{ControlBarrierFunction}, kd::Function)

Solve CBF-QP using OSQP and JuMP.
"""
function solve_cbf_qp(
    x, Σ::ControlAffineSystem, cbfs::Vector{ControlBarrierFunction}, kd::Function
)
    model = Model(OSQP.Optimizer)
    set_silent(model)
    u = Σ.m == 1 ? @variable(model, u) : @variable(model, u[1:(Σ.m)])
    @objective(model, Min, 0.5 * (u - kd(x))' * (u - kd(x)))
    for cbf in cbfs
        @constraint(model, cbf.Lfh(x) + cbf.Lgh(x) * u ≥ -cbf.α(cbf(x)))
    end
    optimize!(model)

    return Σ.m == 1 ? value(u) : value.(u)
end

"""
    solve_cbf_qp(
        x, Σ::ControlAffineSystem, cbfs::Vector{ControlBarrierFunction}, kd::Function, umin, umax
    )

Solve CBF-QP using OSQP and JuMP while satisfying input bounds `umin ≤ u ≤ umax`.
"""
function solve_cbf_qp(
    x,
    Σ::ControlAffineSystem,
    cbfs::Vector{ControlBarrierFunction},
    kd::Function,
    umin,
    umax,
)
    model = Model(OSQP.Optimizer)
    set_silent(model)
    u = Σ.m == 1 ? @variable(model, u) : @variable(model, u[1:(Σ.m)])
    Σ.m == 1 ? @constraint(model, u ≤ umax) : @constraint(model, u .≤ umax)
    Σ.m == 1 ? @constraint(model, u ≥ umin) : @constraint(model, u .≥ umin)
    @objective(model, Min, 0.5 * (u - kd(x))' * (u - kd(x)))
    for cbf in cbfs
        @constraint(model, cbf.Lfh(x) + cbf.Lgh(x) * u ≥ -cbf.α(cbf(x)))
    end
    optimize!(model)

    return Σ.m == 1 ? value(u) : value.(u)
end

"""
    solve_time_varying_cbf_qp(x, t, Σ::ControlAffineSystem, cbfs::Vector{ControlBarrierFunction}, kd::Function)

Solve CBF-QP where desired controller is time-varying
"""
function solve_time_varying_cbf_qp(
    x, t, Σ::ControlAffineSystem, cbfs::Vector{ControlBarrierFunction}, kd::Function
)
    model = Model(OSQP.Optimizer)
    set_silent(model)
    u = Σ.m == 1 ? @variable(model, u) : @variable(model, u[1:(Σ.m)])
    @objective(model, Min, 0.5 * (u - kd(x, t))' * (u - kd(x, t)))
    for cbf in cbfs
        @constraint(model, cbf.Lfh(x) + cbf.Lgh(x) * u ≥ -cbf.α(cbf(x)))
    end
    optimize!(model)

    return Σ.m == 1 ? value(u) : value.(u)
end

"""
    solve_time_varying_cbf_qp(x, t, Σ::ControlAffineSystem, cbfs::Vector{ControlBarrierFunction}, kd::Function, umin, umax)

Solve CBF-QP where desired controller is time-varying with input bounds.
"""
function solve_time_varying_cbf_qp(
    x,
    t,
    Σ::ControlAffineSystem,
    cbfs::Vector{ControlBarrierFunction},
    kd::Function,
    umin,
    umax,
)
    model = Model(OSQP.Optimizer)
    set_silent(model)
    u = Σ.m == 1 ? @variable(model, u) : @variable(model, u[1:(Σ.m)])
    Σ.m == 1 ? @constraint(model, u ≤ umax) : @constraint(model, u .≤ umax)
    Σ.m == 1 ? @constraint(model, u ≥ umin) : @constraint(model, u .≥ umin)
    @objective(model, Min, 0.5 * (u - kd(x, t))' * (u - kd(x, t)))
    for cbf in cbfs
        @constraint(model, cbf.Lfh(x) + cbf.Lgh(x) * u ≥ -cbf.α(cbf(x)))
    end
    optimize!(model)

    return Σ.m == 1 ? value(u) : value.(u)
end
