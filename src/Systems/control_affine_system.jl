"""
    ControlAffineSystem

Models a nonlinear control affine system.

# Fields
- `n::Int`: state dimension
- `m::Int`: control dimension
- `f`: function f(x) that models the drift dynamics
- `g`: function g(x) that models the control directions
"""
struct ControlAffineSystem <: System
    n::Int
    m::Int
    f
    g
end

"""
	Base.step(x, u, t0, tf, Σ::ControlAffineSystem)

Integrate the dynamics of a ControlSystemAffine starting at state x under control u over
from t0 to tf.
"""
function Base.step(Σ::ControlAffineSystem, x, u, t0::Float64, tf::Float64)
	rhs(x, u, t) = Σ.f(x) + Σ.g(x)u
	sol = solve(ODEProblem(rhs, x, (t0, tf), u))
	xnew = Σ.n == 1 ? sol[end] : sol[:,end]

	return xnew
end
