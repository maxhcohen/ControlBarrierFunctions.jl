"""
	CLFQP <: Policy

Struct representing a Control Lyapunov Function based quadratic program (CLFQP).

# Fields
- `control`: function u = control(x) describing the control law.
"""
struct CLFQP <: Policy
	control
end

"""
	(k::CLFQP)(x)

Solve CLF-QP for control input at state x.
"""
function (κ::CLFQP)(x)
	return κ.control(x)
end

"""
	simulate(Σ::ControlAffineSystem, k::CLFQP, x0, time::Tuple)

Simulate the closed-loop trajectory of a control affine system under a CLF-based policy 
from initial condition x0.
"""
function simulate(Σ::ControlAffineSystem, κ::CLFQP, x0, time::NamedTuple)
	ts = time.t0:time.dt:time.tf
	xs = Vector{typeof(x0)}(undef, length(ts))
	xs[1] = x0
	for i in 1:length(ts)-1
		t = ts[i]
		x = xs[i]
		u = κ(x)
		xs[i+1] = step(Σ, x, u, t, t + time.dt)
	end

	return ts, vec2mat(xs)
end

"""
	CLFQP(Σ::ControlAffineSystem, CLF::ControlLyapunovFunction)

Construct a CLF-QP from a ControlAffineSystem and ControlLyapunovFunction.
"""
function CLFQP(Σ::ControlAffineSystem, CLF::ControlLyapunovFunction)
	function control(x)
		model = Model(OSQP.Optimizer)
    	set_silent(model)
    	Σ.m == 1 ? @variable(model, u) : @variable(model, u[1:Σ.m])
		@objective(model, Min, (1/2)u'u)
		@constraint(model, clf, CLF.∇V(x)*(Σ.f(x) + Σ.g(x)*u) <= -CLF.α(CLF(x)))
		optimize!(model)
		
		return Σ.m == 1 ? value(u) : value.(u)
	end

	return CLFQP(control)
end

"""
	CLFQP(Σ::ControlAffineSystem, CLF::ControlLyapunovFunction, U)

Construct a CLF-QP subject to control constraints, where each component of U describes
the limits of each component of u.
"""
function CLFQP(Σ::ControlAffineSystem, CLF::ControlLyapunovFunction, U)
	function control(x)
		model = Model(OSQP.Optimizer)
    	set_silent(model)
    	Σ.m == 1 ? @variable(model, u) : @variable(model, u[1:Σ.m])
		@objective(model, Min, (1/2)u'u)
		@constraint(model, clf, CLF.∇V(x)*(Σ.f(x) + Σ.g(x)*u) <= -CLF.α(CLF(x)))
		if Σ.m == 1
			@constraint(model, U[1] <= u <= U[2])
		else
			for i in 1:Σ.m
				bound = U[i]
				@constraint(model, bound[1] <= u[i] <= bound[2])
			end
		end
		optimize!(model)
		
		return Σ.m == 1 ? value(u) : value.(u)
	end

	return CLFQP(control)
end