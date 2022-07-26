"""
    CLFSontag <: FeedbackController

Control Lyapunov Function (CLF)-bassed controller using Sontag's formula from

E. D. Sontag, "A `universal` construction of Artstein's theorem on nonlinear stabilization,"
Systems & Control Letters, vol. 13, no. 2, pp. 117-123, 1989.
"""
struct CLFSontag <: FeedbackController
    control_law::Function
end

"Evaluate a `CLFSontag` controller at state `x`."
(k::CLFSontag)(x) = k.control_law(x)

"""
    CLFSontag(Σ::ControlAffineSystem, CLF::ControlLyapunovFunction)
    CLFSontag(Σ::ControlAffineSystem, CLF::ControlLyapunovFunction, α::Function)

Construct a `FeedbackController` that uses Sontag's universal CLF formula for stabilization
with or without an additional stability margin `α`.
"""
function CLFSontag(Σ::ControlAffineSystem, CLF::ControlLyapunovFunction)
    function control_law(x)
        a = drift_lie_derivative(CLF, Σ, x)
        b = control_lie_derivative(CLF, Σ, x)'
        u = clf_universal_formula(a, b)

        return u
    end

    return CLFSontag(control_law)
end
function CLFSontag(Σ::ControlAffineSystem, CLF::ControlLyapunovFunction, α::Function)
    function control_law(x)
        a = drift_lie_derivative(CLF, Σ, x) + α(x)
        b = control_lie_derivative(CLF, Σ, x)'
        u = clf_universal_formula(a, b)

        return u
    end

    return CLFSontag(control_law)
end

"""
    clf_universal_formula(a::Float64, b::Union{Float64, Vector{Float64}})

Sontag's universal CLF formula parameterized by `(a,b)`.
"""
function clf_universal_formula(a::Float64, b::Union{Float64, Vector{Float64}})
    if b == zeros(length(b))
        u = zeros(length(b))
    else
        u = -((a + sqrt(a^2 + norm(b)^4))/norm(b)^2) * b
    end

    return u
end