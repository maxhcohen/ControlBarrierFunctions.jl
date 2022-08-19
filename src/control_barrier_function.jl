"""
    ControlBarrierFunction <: CertificateFunction

# Fields
- `h::Function` : function `h(x)` representing the barrier function candidate
- `α::Function` : extended class K function `α(h(x))` defining the CBF
"""
struct ControlBarrierFunction <: CertificateFunction
    h::Function
    α::Function
end

"""
Construct a `ControlBarrierFunction` from  functions `h(x)` and `α(h(x))`.

If no extended class K function is passed in default it to `α(s)=s`.
"""
ControlBarrierFunction(h::Function) = ControlBarrierFunction(h, s -> s)

"Compute gradient of a CBF evaluated at `x`."
function gradient(CBF::ControlBarrierFunction, x)
    return ForwardDiff.gradient(CBF.h, x)'
end

"Compute Lie derivatives of CBF along system dynamics."
function lie_derivatives(CBF::ControlBarrierFunction, Σ::ControlAffineSystem, x)
    return [drift_lie_derivative(CBF, Σ, x), control_lie_derivative(CBF, Σ, x)]
end

"Compute Lie derivative of CBF along drift dynamics ``Lfh(x) = ∇h(x) * f(x)``."
function drift_lie_derivative(CBF::ControlBarrierFunction, Σ::ControlAffineSystem, x)
    return gradient(CBF, x) * Σ.f(x)
end

"Compute Lie derivative of CBF along control directions ``Lgh(x) = ∇h(x) * g(x)``."
function control_lie_derivative(CBF::ControlBarrierFunction, Σ::ControlAffineSystem, x)
    return gradient(CBF, x) * Σ.g(x)
end