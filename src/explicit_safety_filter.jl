"""
    ExplicitSafetyFilter <: SafetyFilter

Controller that uses the closed-form solution to a control barrier function quadratic program.

# Fields
- `k::Function` : function that computes safe control actions
"""
struct ExplicitSafetyFilter{T} <: SafetyFilter
    k::T
end

"""
    (k::ExplicitSafetyFilter)(x)

Functors for evaluating explicit safety filter
"""
(k::ExplicitSafetyFilter)(args...) = k.k(args...)

"""
    ExplicitSafetyFilter(cbf::ControlBarrierFunction, Σ::ControlAffineSystem, kd::Function)

Construct an ExplicitSafetyFilter from a cbf and a desired controller.
"""
function ExplicitSafetyFilter(
    cbf::ControlBarrierFunction, Σ::ControlAffineSystem, kd::Function
)
    function k(x, args...)
        Lgh = cbf.Lgh(x)
        Lfh = cbf.Lfh(x)
        kdx = kd(x, args...)
        a = Lfh + Lgh * kdx + cbf.α(cbf(x))
        kdx + λQP(a, norm(Lgh)^2) * Lgh'
    end
    ExplicitSafetyFilter{typeof(k)}(k)
end

"""
    ExplicitSafetyFilter(cbf::ControlBarrierFunction, Σ::ControlAffineSystem)

If no desired controller passed in then default it to zero.
"""
function ExplicitSafetyFilter(cbf::ControlBarrierFunction, Σ::ControlAffineSystem)
    kd(x) = Σ.m == 1 ? 0.0 : zeros(Σ.m)

    return ExplicitSafetyFilter(cbf, Σ, kd)
end

# Some helper functions
ReLU(x::Float64) = max(0.0, x)
λQP(a::Float64, b::Float64) = b == 0.0 ? 0.0 : ReLU(-a / b)
