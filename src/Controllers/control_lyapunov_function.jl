struct ControlLyapunovFunction <: CertificateFunction
    V
    α
end

(V::ControlLyapunovFunction)(x) = V.V(x)

function ControlLyapunovFunction(Σ::ControlAffineSystem, Q, R, α)
    # Linearize ControlAffineSystem
    A, B = linearize(Σ, zeros(state_dim(Σ)))

    # Compute value function of LQR problem
    P = are(Continuous, A, B, Q, R)

    return ControlLyapunovFunction(x -> x'P*x, α)
end

function lie_derivatives(V::ControlLyapunovFunction, Σ::ControlAffineSystem, x)
    ∇V = gradient(x -> V(x), x)[1]'
    LfV = ∇V * _f(Σ, x)
    LgV = ∇V * _g(Σ, x)

    return LfV, LgV
end
