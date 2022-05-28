struct ControlBarrierFunction <: CertificateFunction
    h
    α
end

(h::ControlBarrierFunction)(x) = h.h(x)

function lie_derivatives(h::ControlBarrierFunction, Σ::ControlAffineSystem, x)
    ∇h = gradient(x -> h(x), x)[1]'
    Lfh = ∇h * _f(Σ, x)
    Lgh = ∇h * _g(Σ, x)

    return Lfh, Lgh
end