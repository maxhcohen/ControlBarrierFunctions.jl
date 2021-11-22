"""
    ControlBarrierFunction

Control barrier function h(x) for a safe set C.
# Fields
- `h`: function h(x) that represents the CBF
- `∇h`: function ∇h(x) defining the gradient of the CBF
- `α`: extended class K function α(h(x))
"""
struct ControlBarrierFunction
    h
    ∇h
    α
end

"Evaluate CBF at state x"
(cbf::ControlBarrierFunction)(x) = cbf.h(x)

"""
    ControlBarrierFunction(h; α)

Construct a CBF given h(x) and extended class K function α(h(x)).
Extended classK function α() defaults to the identity i.e., α(h)=h.
"""
function ControlBarrierFunction(h; α=r->r)
    ∇h(x) = ForwardDiff.gradient(h,x)'
    return ControlBarrierFunction(h, ∇h, α)
end