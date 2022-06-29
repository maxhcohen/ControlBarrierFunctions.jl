"""
    ControlLyapunovFunction <: CertificateFunction

Control Lyapunov Function type with Lyapunov function candidate V(x) and class K function α/
"""
struct ControlLyapunovFunction <: CertificateFunction
    V
    α
end

"""
    (V::ControlLyapunovFunction)(x)

Evaluate CLF at state x.
"""
(V::ControlLyapunovFunction)(x) = V.V(x)

"""
    _∇V(V::ControlLyapunovFunction, x)

Compute gradient of CLF at x.
"""
function _∇V(V::ControlLyapunovFunction, x)
    _V(x) = V(x)
    
    return ForwardDiff.gradient(_V, x)'
end

"""
    _LfV(V::ControlLyapunovFunction, Σ::ControlAffineSystem, x)

Compute the Lie derivative of V along the drift vector field.
"""
function _LfV(V::ControlLyapunovFunction, Σ::ControlAffineSystem, x)
    return _∇V(V, x) * _f(Σ, x)
end

"""
    _LgV(V::ControlLyapunovFunction, Σ::ControlAffineSystem, x)

Compute the Lie derivative of V along the vector field of control directions.
"""
function _LgV(V::ControlLyapunovFunction, Σ::ControlAffineSystem, x)
    return _∇V(V, x) * _g(Σ, x)
end

"""
    _Lf0V(V::ControlLyapunovFunction, Σ::ControlAffineSystem, x)

Compute the Lie derivative of V along the nominal drift dynamics. 
"""
function _Lf0V(V::ControlLyapunovFunction, Σ::ControlAffineSystem, x)
    _∇V(V, x) * _f0(Σ, x)
end

function _LFV(V::ControlLyapunovFunction, Σ::ControlAffineSystem, x)
    _∇V(V, x) * _F(Σ, x)
end

function lie_derivatives(V::ControlLyapunovFunction, Σ::ControlAffineSystem, x)
    return _LfV(V, Σ, x), _LgV(V, Σ, x)
end

function clf_condition(V::ControlLyapunovFunction, Σ::ControlAffineSystem, x, u)
    return _LfV(V, Σ, x) + _LgV(V, Σ, x)*u + V.α(V(x))
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


