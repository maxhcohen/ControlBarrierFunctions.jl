"""
    ControlBarrierFunction

Control barrier function (CBF) defining a safe set as its zero superlevel set.

# Fields
- `h::Function` : function defining the safe set `h(x) ≥ 0`
- `α::Function` : extended class K function `a(h(x))` for CBF
- `∇h::Function` : gradient of CBF
- `Lfh::Function` : Lie derivative of CBF along drift vector field `f`
- `Lgh::Function` : Lie derivative of CBF along control directions `g`
"""
struct ControlBarrierFunction
    h::Function
    α::Function
    ∇h::Function
    Lfh::Function
    Lgh::Function
end

"""
    (cbf::ControlBarrierFunction)(x)

Functor for evaluating CBF at state `x`.
"""
(cbf::ControlBarrierFunction)(x) = cbf.h(x)

"""
    ControlBarrierFunction(h::Function, Σ::ControlAffineSystem, α::Function)

Construct a CBF from a function `h`, a control affine system `Σ`, and an extended class K function `α`.
"""
function ControlBarrierFunction(h::Function, Σ::ControlAffineSystem, α::Function)
    ∇h(x) = Σ.n == 1 ? ForwardDiff.derivative(h, x) : ForwardDiff.gradient(h, x)
    Lfh(x) = ∇h(x)'*Σ.f(x)
    Lgh(x) = ∇h(x)'*Σ.g(x)

    return ControlBarrierFunction(h, α, ∇h, Lfh, Lgh)
end

"""
    ControlBarrierFunction(h::Function, Σ::ControlAffineSystem)

If no extended class K function provided, default to the identify function.
"""
ControlBarrierFunction(h::Function, Σ::ControlAffineSystem) = ControlBarrierFunction(h, Σ, r -> r)