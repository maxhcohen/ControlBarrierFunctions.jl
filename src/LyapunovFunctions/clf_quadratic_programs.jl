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
	CLFQP(Σ::ControlAffineSystem, CLF::ControlLyapunovFunction)
    CLFQP(Σ::ControlAffineSystem, CLF::ControlLyapunovFunction, A, b)

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

function CLFQP(Σ::ControlAffineSystem, CLF::ControlLyapunovFunction, A, b)
	function control(x)
		model = Model(OSQP.Optimizer)
    	set_silent(model)
    	Σ.m == 1 ? @variable(model, u) : @variable(model, u[1:Σ.m])
		@objective(model, Min, (1/2)u'u)
		@constraint(model, clf, CLF.∇V(x)*(Σ.f(x) + Σ.g(x)*u) <= -CLF.α(CLF(x)))
        @constraint(model, A*u .<= b)
		optimize!(model)

	return Σ.m == 1 ? value(u) : value.(u)
	end

	return CLFQP(control)
end
