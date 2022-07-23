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

# Constructors
ControlLyapunovFunction(V::Function, α::Function) = ControlLyapunovFunction(V, α, false, 0.0)
ControlLyapunovFunction(V::Function) = ControlLyapunovFunction(V, x -> V(x), false, 0.0)

# Stability margin if we were to use Sontag's formula
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

# Compute gradient and Lie derivatives
function gradient(CLF::ControlLyapunovFunction, x)
    return ForwardDiff.gradient(CLF.V, x)'
end

function lie_derivatives(CLF::ControlLyapunovFunction, Σ::ControlAffineSystem, x)
    return [drift_lie_derivative(CLF, Σ, x), control_lie_derivative(CLF, Σ, x)]
end

function drift_lie_derivative(CLF::ControlLyapunovFunction, Σ::ControlAffineSystem, x)
    return gradient(CLF, x) * Σ.f(x)
end

function control_lie_derivative(CLF::ControlLyapunovFunction, Σ::ControlAffineSystem, x)
    return gradient(CLF, x) * Σ.g(x)
end
