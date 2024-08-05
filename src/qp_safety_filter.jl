"""
    QPSafetyFilter

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
    QPSafetyFilter(cbfs::Vector{ControlBarrierFunction}, Σ::ControlAffineSystem, kd::Function; tunable=false)

Construct an QPSafetyFilter from a cbf and a desired controller.

# Keywork arguments
- `tunable::Bool` : boolean to decide if coefficients on extended class K functions should be decision variables
"""
function QPSafetyFilter(cbfs::Vector{ControlBarrierFunction}, Σ::ControlAffineSystem, kd::Function)
    try
        kd(Σ.n == 1 ? rand() : rand(Σ.n), 0.0) # See if desired controller is time-varying
    catch e
        if isa(e, MethodError) # If controller is not time-varying 
            return QPSafetyFilter(x -> solve_cbf_qp(x, Σ, cbfs, kd))
        else
           return e
        end
    else # If controller is time-varying
        return QPSafetyFilter((x,t) -> solve_time_varying_cbf_qp(x, t, Σ, cbfs, kd))
    end
end

"""
    QPSafetyFilter(cbf::ControlBarrierFunction, Σ::ControlAffineSystem, kd::Function; kwargs...)

Add ability to pass in single CBF.
"""
QPSafetyFilter(cbf::ControlBarrierFunction, Σ::ControlAffineSystem, kd::Function) = QPSafetyFilter([cbf], Σ, kd)

"""
    solve_cbf_qp(x, Σ::ControlAffineSystem, cbfs::Vector{ControlBarrierFunction}, kd::Function)

Solve CBF-QP
"""
function solve_cbf_qp(x, Σ::ControlAffineSystem, cbfs::Vector{ControlBarrierFunction}, kd::Function)
    model = Model(OSQP.Optimizer)
    set_silent(model)
    u = Σ.m == 1 ? @variable(model, u) : @variable(model, u[1:Σ.m])
    @objective(model, Min, 0.5*(u - kd(x))'*(u - kd(x)))
    for cbf in cbfs
        @constraint(model, cbf.Lfh(x) + cbf.Lgh(x)*u ≥ -cbf.α(cbf(x)))
    end
    optimize!(model)

    return Σ.m == 1 ? value(u) : value.(u)
end

"""
    solve_time_varying_cbf_qp(x, t, Σ::ControlAffineSystem, cbfs::Vector{ControlBarrierFunction}, kd::Function)

Solve CBF-QP where desired controller is time-varying
"""
function solve_time_varying_cbf_qp(x, t, Σ::ControlAffineSystem, cbfs::Vector{ControlBarrierFunction}, kd::Function)
    model = Model(OSQP.Optimizer)
    set_silent(model)
    u = Σ.m == 1 ? @variable(model, u) : @variable(model, u[1:Σ.m])
    @objective(model, Min, 0.5*(u - kd(x,t))'*(u - kd(x,t)))
    for cbf in cbfs
        @constraint(model, cbf.Lfh(x) + cbf.Lgh(x)*u ≥ -cbf.α(cbf(x)))
    end
    optimize!(model)

    return Σ.m == 1 ? value(u) : value.(u)
end