struct ConfigurationError <: Output
    qd
    relative_degree
end

function (y::ConfigurationError)(Σ::ControlAffineSystem, x)
    N = control_dim(Σ)
    if N == 1
        e = x[1] - y.qd
    else
        e = x[1:N] - y.qd
    end
    
    return e
end

function _dy(y::ConfigurationError, Σ::ControlAffineSystem, x)
    _y(x) = y(Σ, x)
    if control_dim(Σ) == 1
        out = ForwardDiff.gradient(_y, x)'
    else
        out = ForwardDiff.jacobian(_y, x)
    end

    return out
end

function _Lfy(y::ConfigurationError, Σ::ControlAffineSystem, x)
    return _dy(y, Σ, x) * _f(Σ, x)
end

function _dLfy(y::ConfigurationError, Σ::ControlAffineSystem, x)
    Lfy(x) = _Lfy(y, Σ, x)
    if control_dim(Σ) == 1
        out = ForwardDiff.gradient(Lfy, x)'
    else
        out = ForwardDiff.jacobian(Lfy, x)
    end
    
    return out
end

function _Lf²y(y::ConfigurationError, Σ::ControlAffineSystem, x)
    return _dLfy(y, Σ, x) * _f(Σ, x)
end

function _LgLfy(y::ConfigurationError, Σ::ControlAffineSystem, x)
    return _dLfy(y, Σ, x) * _g(Σ, x)
end

function lie_derivatives(y::ConfigurationError, Σ::ControlAffineSystem, x)
    if y.relative_degree == 2
        out1 = _Lf²y(y, Σ, x)
        out2 = _LgLfy(y, Σ, x)
    end

    return out1, out2
end


