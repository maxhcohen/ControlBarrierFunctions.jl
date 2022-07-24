"""
    ControlLyapunovFunction <: CertificateFunction

Control Lyapunov Function (CLF) for a control affine system.

# Fields
- `V`: function V(x) that represents the Lyapunov function candidate
- `α`: positive definite function α(x) that represents the desired rate of CLF decay
- `relax` : bool that indicates if we want to relax the CLF constraint in a QP
- `p`: relaxation penality on CLF constraints. Only relevant if relax = true
"""
struct ControlLyapunovFunction <: CertificateFunction
    V::Function
    α::Function
    relax::Bool
    p::Real
end

"Constructors for `ControlLyapunovFunction`. If `α` not passed in, set it to `α(x)=V(x)`."
ControlLyapunovFunction(V::Function, α::Function) = ControlLyapunovFunction(V, α, false, 0.0)
ControlLyapunovFunction(V::Function) = ControlLyapunovFunction(V, x -> V(x), false, 0.0)

"Set CLF stability margin to what we would obtain if using Sontag's formula."
function ControlLyapunovFunction(V::Function, Σ::ControlAffineSystem)
    CLF = ControlLyapunovFunction(V)
    a(x) = drift_lie_derivative(CLF, Σ, x)
    b(x) = control_lie_derivative(CLF, Σ, x)
    α(x) = sqrt(a(x)^2 + norm(b(x))^4)

    return ControlLyapunovFunction(V, α)
end
function ControlLyapunovFunction(V::Function, Σ::ControlAffineSystem, relax::Bool, p::Real)
    CLF = ControlLyapunovFunction(V)
    a(x) = drift_lie_derivative(CLF, Σ, x)
    b(x) = control_lie_derivative(CLF, Σ, x)
    α(x) = sqrt(a(x)^2 + norm(b(x))^4)

    return ControlLyapunovFunction(V, α, relax, p)
end

"Compute gradient of a CLF evaluated at `x`."
function gradient(CLF::ControlLyapunovFunction, x)
    return ForwardDiff.gradient(CLF.V, x)'
end

"Compute Lie derivatives of CLF along system dynamics."
function lie_derivatives(CLF::ControlLyapunovFunction, Σ::ControlAffineSystem, x)
    return [drift_lie_derivative(CLF, Σ, x), control_lie_derivative(CLF, Σ, x)]
end

"Compute Lie derivative of CLF along drift dynamics ``LfV(x) = ∇V(x) * f(x)``."
function drift_lie_derivative(CLF::ControlLyapunovFunction, Σ::ControlAffineSystem, x)
    return gradient(CLF, x) * Σ.f(x)
end

"Compute Lie derivative of CLF along control directions ``LgV(x) = ∇V(x) * g(x)``."
function control_lie_derivative(CLF::ControlLyapunovFunction, Σ::ControlAffineSystem, x)
    return gradient(CLF, x) * Σ.g(x)
end
