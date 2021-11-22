"""
    ControlLyapunovFunction

Control Lyapunov function V for a control affine system
# Fields
- `V`: function V(x) that represents the CLF
- `∇V`: function ∇V(x) defining the gradient of the CLF
- `α`: class K function α(V(x))
"""
struct ControlLyapunovFunction
    V
    ∇V
    α
end

"Evaluate CLF at state x"
(clf::ControlLyapunovFunction)(x) = clf.V(x)

"""
    ControlLyapunovFunction(V; α)

Construct a CLF given V(x) and class K function α(V(x)).
Extended classK function α() defaults to the identity i.e., α(V)=V.
"""
function ControlLyapunovFunction(V; α=r->r)
    ∇V(x) = ForwardDiff.gradient(V,x)'
    return ControlBarrierFunction(V, ∇V, α)
end