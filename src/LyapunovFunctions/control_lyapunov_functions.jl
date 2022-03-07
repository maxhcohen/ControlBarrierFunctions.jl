"Abstract Lyapunov function type"
abstract type LyapunovFunction <: CertificateFunction end

"""
    ControlLyapunovFunction

Control Lyapunov function V for a control affine system
# Fields
- `V`: function V(x) that represents the CLF
- `∇V`: function ∇V(x) defining the gradient of the CLF
- `α`: class K function α(V(x))
"""
struct ControlLyapunovFunction <: LyapunovFunction
    V
    ∇V
    α
end

"""
    ControlLyapunovFunction(V, α)

Construct a CLF given V(x) and class K function α(V(x)).
"""
function ControlLyapunovFunction(V, α)
    ∇V(x) = gradient(V, x)[1]'
    return ControlLyapunovFunction(V, ∇V, α)
end

"""
	(V::ControlLyapunovFunction)(x)

Evaluate CLF at state x.
"""
function (V::ControlLyapunovFunction)(x)
	return V.V(x)
end
