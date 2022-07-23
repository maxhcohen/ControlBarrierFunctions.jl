"""
    LQRController <: FeedbackController

# Fields
- `P` : positive definite matrix associated with the value function x'P*x
- `K` : optimal LQR gain
"""
struct LQRController <: FeedbackController
    Q
    R
    x0
    P
    K
end

(k::LQRController)(x) = size(k.K)[1] == 1 ? dot(-k.K, x - k.x0) : -k.K * (x - k.x0)

"""
    LQRController(Σ::ControlAffineSystem, Q, R)
    LQRController(Σ::ControlAffineSystem, Q, R, x0)

LQRController constructors.
"""
function LQRController(Σ::ControlAffineSystem, Q, R)
    x0 = zeros(Σ.n)
    A, B = linearize(Σ, x0)
    K = lqr(Continuous, A, B, Q, R)
    P = are(Continuous, A, B, Q, R)

    return LQRController(Q, R, x0, P, K)
end

function LQRController(Σ::ControlAffineSystem, Q, R, x0)
    A, B = linearize(Σ, x0)
    K = lqr(Continuous, A, B, Q, R)
    P = are(Continuous, A, B, Q, R)

    return LQRController(Q, R, x0, P, K)
end

"""
    linearize(Σ::ControlAffineSystem, x0)

Linearize control affine system about point x0.
"""
function linearize(Σ::ControlAffineSystem, x0)
    A = ForwardDiff.jacobian(Σ.f, x0)
    B = Σ.g(x0)

    return A, B
end

