"Abstract barrier function type"
abstract type BarrierFunction <: CertificateFunction end

"""
    ControlBarrierFunction

Control barrier function h for a control affine system.
# Fields
- `h`: function h(x) that represents the CBF
- `∇h`: function ∇h(x) defining the gradient of the CBF
- `α`: extended class K function α(h(x))
"""
struct ControlBarrierFunction <: BarrierFunction
    h
    ∇h
    α
end

"""
	ControlBarrierFunction(h, α)

Construct a CBF given h(x) and extended class K function α(h(x)).
"""
function ControlBarrierFunction(h, α)
    ∇h(x) = gradient(h, x)[1]'
    return ControlBarrierFunction(h, ∇h, α)
end

"""
	(h::ControlBarrierFunction)(x)

Evaluate CBF at state x.
"""
function (h::ControlBarrierFunction)(x)
    return h.h(x)
end
