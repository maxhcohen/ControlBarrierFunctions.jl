"""
    CLFSontag <: FeedbackController

Control Lyapunov Function (CLF)-bassed controller using Sontag's formula from

E. D. Sontag, "A `universal` construction of Artstein's theorem on nonlinear stabilization,"
Systems & Control Letters, vol. 13, no. 2, pp. 117-123, 1989.
"""
struct CLFSontag <: FeedbackController
    control_law::Function
end

(k::CLFSontag)(x) = k.control_law(x)

function CLFSontag(Σ::ControlAffineSystem, CLF::ControlLyapunovFunction)
    function control_law(x)
        a = drift_lie_derivative(CLF, Σ, x)
        b = control_lie_derivative(CLF, Σ, x)'
        if b == zeros(Σ.m)
            u = zeros(Σ.m)
        else
            u = -((a + sqrt(a^2 + norm(b)^4))/norm(b)^2) * b
        end

        return u
    end

    return CLFSontag(control_law)
end