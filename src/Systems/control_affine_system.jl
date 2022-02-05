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

"""
    (S::Simulation)(Σ::ControlAffineSystem, x0)

Simulate an open-loop trajectory of a control affine system from initial condition x0.
"""
function (S::Simulation)(Σ::ControlAffineSystem, x0)
    xs = Vector{typeof(x0)}(undef, length(S.ts))
	xs[1] = x0
	for i in 1:length(S.ts)-1
		t = S.ts[i]
		x = xs[i]
		u = Σ.m == 1 ? 0.0 : zeros(Σ.m)
		xs[i+1] = step(Σ, x, u, t, t + S.dt)
	end

	return vec2mat(xs)
end
