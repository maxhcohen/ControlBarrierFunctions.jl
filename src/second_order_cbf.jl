"""
    SecondOrderCBF <: HighOrderCBF

High order control barrier function (HOCBF) with relative degree 2.

# Fields
- `h::function` : constraint function representing the constraint set
- `ψ1::Function` : dynamically extended constraint function
- `α1::Function` : extended class K function
- `α2::Function` : extended class K function
"""
struct SecondOrderCBF <: HighOrderCBF
    h::Function
    ψ1::Function
    α1::Function
    α2::Function
end

"SecondSecondOrderCBF constructors."
function SecondOrderCBF(Σ::ControlAffineSystem, h::Function, α1::Function, α2::Function)
    ψ1(x) = ForwardDiff.gradient(h, x)' * Σ.f(x) + α1(h(x))

    return SecondOrderCBF(h, ψ1, α1, α2)
end

function SecondOrderCBF(Σ::ControlAffineSystem, h::Function)
    α1(s) = s
    α2(s) = s

    return SecondOrderCBF(Σ, h, α1, α2)
end 

function SecondOrderCBF(Σ::ControlAffineSystem, h::Function, α::Function)
    α1 = α
    α2 = α

    return SecondOrderCBF(Σ, h, α1, α2)
end 

"Compute gradient of a HOCBF evaluated at `x`."
function gradient(HOCBF::SecondOrderCBF, x)
    return ForwardDiff.gradient(HOCBF.ψ1, x)'
end

"Compute Lie derivative of HOCBF along drift dynamics."
function drift_lie_derivative(HOCBF::SecondOrderCBF, Σ::ControlAffineSystem, x)
    return gradient(HOCBF, x) * Σ.f(x)
end

"Compute Lie derivative of HOCBF along control directions."
function control_lie_derivative(HOCBF::SecondOrderCBF, Σ::ControlAffineSystem, x)
    return gradient(HOCBF, x) * Σ.g(x)
end

