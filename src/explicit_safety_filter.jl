"""
    ExplicitSafetyFilter <: SafetyFilter

Controller that uses the closed-form solution to a control barrier function quadratic program.

# Fields
- `k::Function` : function that computes safe control actions
"""
struct ExplicitSafetyFilter <: SafetyFilter
    k::Function
end

"""
    (k::ExplicitSafetyFilter)(x)

Functors for evaluating explicit safety filter
"""
(k::ExplicitSafetyFilter)(x) = k.k(x)
(k::ExplicitSafetyFilter)(x, t) = k.k(x, t)

"""
    ExplicitSafetyFilter(cbf::ControlBarrierFunction, Σ::ControlAffineSystem, kd::Function)

Construct an ExplicitSafetyFilter from a cbf and a desired controller.
"""
function ExplicitSafetyFilter(
    cbf::ControlBarrierFunction, Σ::ControlAffineSystem, kd::Function
)
    try
        kd(Σ.n == 1 ? rand() : rand(Σ.n), 0.0)
    catch e
        if isa(e, MethodError)
            a(x) = cbf.Lfh(x) + cbf.Lgh(x) * kd(x) + cbf.α(cbf(x))
            k(x) = kd(x) + λQP(a(x), norm(cbf.Lgh(x))^2) * cbf.Lgh(x)'

            return ExplicitSafetyFilter(k)
        else
            return e
        end
    else
        a(x, t) = cbf.Lfh(x) + cbf.Lgh(x) * kd(x, t) + cbf.α(cbf(x))
        k(x, t) = kd(x, t) + λQP(a(x, t), norm(cbf.Lgh(x))^2) * cbf.Lgh(x)'

        return ExplicitSafetyFilter(k)
    end
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
