"""
    TimeVaryingCLF

Time varying Control Lyapunov function V for a control affine system
# Fields
- `V`: function V(x, t) that represents the CLF
- `∇V`: function ∇V(x, t) defining the spatial gradient of the CLF
- `∂ₜV`: function ∂ₜV(x, t) defining the time derivative of the CLF
- `α`: class K function α(V(x, t))
"""
struct TimeVaryingCLF <: LyapunovFunction
    V
    ∇V
    ∂ₜV
    α
end

"""
    TimeVaryingCLF(V, α)

Construct a CLF given V and class K function α.
"""
function TimeVaryingCLF(V, α)
    ∇V(x, t) = gradient((r, s) -> V(r, s), x, t)[1]'
    ∂ₜV(x, t) = gradient((r, s) -> V(r, s), x, t)[2]
    return TimeVaryingCLF(V, ∇V, ∂ₜV, α)
end

"""
	(V::TimeVaryingCLF)(x, t)

Evaluate CLF at state x and time t.
"""
function (V::TimeVaryingCLF)(x, t)
    return V.V(x, t)
end