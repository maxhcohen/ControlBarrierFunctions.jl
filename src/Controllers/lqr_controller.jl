struct LQRController <: FeedbackController
    Q
    R
    A
    B
    P
    K
end

function LQRController(Σ::ControlAffineSystem, Q, R)
    # Linearize ControlAffineSystem
    A, B = linearize(Σ, zeros(state_dim(Σ)))

    # Compute LQR feedback gain
    K = lqr(Continuous, A, B, Q, R)

    # Compute associated value function: V(x) = x'Px
    P = are(Continuous, A, B, Q, R)

    return LQRController(Q, R, A, B, P, K)
end

function LQRController(Σ::ControlAffineSystem)
    # Linearize ControlAffineSystem
    A, B = linearize(Σ, zeros(state_dim(Σ)))

    # Default Q and R matrices
    Q = Matrix(1.0I, state_dim(Σ), state_dim(Σ))
    R = Matrix(1.0I, control_dim(Σ), control_dim(Σ))

    # Compute LQR feedback gain
    K = lqr(Continuous, A, B, Q, R)


    # Compute associated value function: V(x) = x'Px
    P = are(Continuous, A, B, Q, R)

    return LQRController(Q, R, A, B, P, K)
end

(k::LQRController)(x) = size(k.K)[1] == 1 ? dot(-k.K, x) : -k.K * x