struct FLController <: FeedbackController
    K
end

function FLController(Σ::ControlAffineSystem)
    # Get dimensions of state and output
    n = state_dim(Σ)
    m = control_dim(Σ)

    # Linear system after feedback linearization
    A = [zeros(m, m) Matrix(1.0I, m, m); zeros(m, m) zeros(m, m)]
    B = vcat(zeros(m, m), Matrix(1.0I, m, m))

    # Default Q and R matrices
    Q = Matrix(1.0I, n, n)
    R = Matrix(1.0I, m, m)

    # Compute LQR feedback gain
    K = lqr(Continuous, A, B, Q, R)

    return FLController(K)
end

function FLController(Σ::ControlAffineSystem, Q, R)
    # Get dimensions of state and output
    n = state_dim(Σ)
    m = control_dim(Σ)

    # Linear system after feedback linearization
    A = [zeros(m, m) Matrix(1.0I, m, m); zeros(m, m) zeros(m, m)]
    B = vcat(zeros(m, m), Matrix(1.0I, m, m))

    # Compute LQR feedback gain
    K = lqr(Continuous, A, B, Q, R)

    return FLController(K)
end

function (k::FLController)(Σ::ControlAffineSystem, y::ConfigurationError, x)
    if y.relative_degree == 1
        Lfy, Lgy = lie_derivatives(y, Σ, x)
        ξ = y(Σ, x)
        ν = size(k.K)[1] == 1 ? dot(-k.K, ξ) : -k.K * ξ
        u = -inv(Lgy)*(Lfy - ν)
    elseif y.relative_degree == 2
        Lf²y, LfLgy = lie_derivatives(y, Σ, x)
        ξ = vcat(y(Σ, x), _Lfy(y, Σ, x))
        ν = size(k.K)[1] == 1 ? dot(-k.K, ξ) : -k.K * ξ
        u = -inv(LfLgy)*(Lf²y - ν)
    end

    return u
end