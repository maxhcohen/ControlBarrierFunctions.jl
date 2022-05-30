struct ControlLyapunovFunction <: CertificateFunction
    V
    α
end

(V::ControlLyapunovFunction)(x) = V.V(x)

function lie_derivatives(V::ControlLyapunovFunction, Σ::ControlAffineSystem, x)
    ∇V = gradient(x -> V(x), x)[1]'
    LfV = ∇V * _f(Σ, x)
    LgV = ∇V * _g(Σ, x)

    return LfV, LgV
end

function ControlLyapunovFunction(Σ::ControlAffineSystem, Q, R, α)
    # Linearize ControlAffineSystem
    A, B = linearize(Σ, zeros(state_dim(Σ)))

    # Compute value function of LQR problem
    P = are(Continuous, A, B, Q, R)

    return ControlLyapunovFunction(x -> x'P*x, α)
end

function ControlLyapunovFunction(Σ::ControlAffineSystem, Λ, α)
    # Check if ControlAffineSystem is an Euler-Lagrange system
    n = state_dim(Σ)
    N = degrees_of_freedom(Σ)
    if euler_lagrange(Σ::ControlAffineSystem)
        function lyap(x)
            q = n == 1 ? x[1] : x[1:N]
            q̇ = n == 1 ? x[2] : x[N+1:end]
            r = q̇ + Λ*q
            M = inertia_matrix(Σ, q)

            return 0.5 * r'M*r + 0.5 * q'q
        end
    else
        pritnln("Please make sure system is Euler-Lagrange.")
    end

    return ControlLyapunovFunction(x -> lyap(x), α)
end


