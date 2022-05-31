struct ControlBarrierFunction <: CertificateFunction
    h
    α
    relative_degree
end

ControlBarrierFunction(h, α) = ControlBarrierFunction(h, α, 1)

(h::ControlBarrierFunction)(x) = h.h(x)

function _ψ0(h::ControlBarrierFunction, x)
    return h(x)
end

function _∇h(h::ControlBarrierFunction, x)
    _h(x) = h(x)

    return ForwardDiff.gradient(_h, x)'
end

function _dα(h::ControlBarrierFunction, r)
    _α(r) = h.α(r)

    return ForwardDiff.derivative(_α, r)
end

function _Lfh(h::ControlBarrierFunction, Σ::ControlAffineSystem, x)
    return _∇h(h, x) * _f(Σ, x)
end

function _Lgh(h::ControlBarrierFunction, Σ::ControlAffineSystem, x)
    return _∇h(h, x) * _g(Σ, x)
end

function _dψ0(h::ControlBarrierFunction, Σ::ControlAffineSystem, x)
    return _Lfh(h, Σ, x)
end

function _ψ1(h::ControlBarrierFunction, Σ::ControlAffineSystem, x)
    return _dψ0(h, Σ, x) + h.α(_ψ0(h, x))
end

function _dLfh(h::ControlBarrierFunction, Σ::ControlAffineSystem, x)
    Lfh(x) = _Lfh(h, Σ, x)

    return ForwardDiff.gradient(Lfh, x)'
end

function _Lf²h(h::ControlBarrierFunction, Σ::ControlAffineSystem, x)
    return _dLfh(h, Σ, x) * _f(Σ, x)
end

function _LgLfh(h::ControlBarrierFunction, Σ::ControlAffineSystem, x)
    return _dLfh(h, Σ, x) * _g(Σ, x)
end

function _dψ1(h::ControlBarrierFunction, Σ::ControlAffineSystem, x)
    return _Lf²h(h, Σ, x) + _dα(h, h(x))*_Lfh(h, Σ, x)
end

function _dψ1(h::ControlBarrierFunction, Σ::ControlAffineSystem, x, u)
    return _Lf²h(h, Σ, x) + _LgLfh(h, Σ, x)*u +  _dα(h, h(x))*_Lfh(h, Σ, x)
end

function lie_derivatives(h::ControlBarrierFunction, Σ::ControlAffineSystem, x)
    r = h.relative_degree
    if r == 1
        out1 = _Lfh(h, Σ, x)
        out2 = _Lgh(h, Σ, x)
    elseif r == 2
        out1 = _Lf²h(h, Σ, x)
        out2 = _LgLfh(h, Σ, x)
    end

    return out1, out2
end

function cbf_condition(h::ControlBarrierFunction, Σ::ControlAffineSystem, x, u)
    r = h.relative_degree
    if r == 1
        out = _Lfh(h, Σ, x) + _Lgh(h, Σ, x)*u + h.α(h(x))
    elseif r == 2
        out = _dψ1(h, Σ, x, u) + h.α(_ψ1(h, Σ, x))
    end

    return out
end