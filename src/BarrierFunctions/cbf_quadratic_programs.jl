"""
	CBFQP <: Policy

Struct representing a Control Barrier Function based quadratic program (CBFQP).

# Fields
- `control`: function u = control(x) describing the control law.
"""
struct CBFQP <: Policy
	control
end

"""
	(κ::CBFQP)(x)

Solve CBF-QP at state x.
"""
function (κ::CBFQP)(x)
	return κ.control(x)
end

# Main function used to simulate systems under a CBFQP policy.
"""
	simulate(Σ::ControlAffineSystem, κ::CBFQP, x0, time::Tuple)

Simulate the closed-loop trajectory of a control affine system under a CBF-based policy 
from initial condition x0.
"""
function simulate(Σ::ControlAffineSystem, κ::CBFQP, x0, time::NamedTuple)
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
	CBFQP(Σ::ControlAffineSystem, CBF::ControlBarrierFunction)

Construct a CBF-QP from a ControlAffineSystem and ControlBarrierFunction.
"""
function CBFQP(Σ::ControlAffineSystem, CBF::ControlBarrierFunction)
	function control(x)
		model = Model(OSQP.Optimizer)
    	set_silent(model)
    	Σ.m == 1 ? @variable(model, u) : @variable(model, u[1:Σ.m])
		@objective(model, Min, (1/2)u'u)
		@constraint(model, cbf, CBF.∇h(x)*(Σ.f(x) + Σ.g(x)*u) >= -CBF.α(CBF(x)))
		optimize!(model)
		
		return Σ.m == 1 ? value(u) : value.(u)
	end

	return CBFQP(control)
end


"""
	CBFQP(Σ::ControlAffineSystem, CBFs::Vector{ControlBarrierFunction})

Construct a CBF-QP from a ControlAffineSystem and ControlBarrierFunction.
"""
function CBFQP(Σ::ControlAffineSystem, CBFs::Vector{ControlBarrierFunction})
	function control(x)
		model = Model(OSQP.Optimizer)
    	set_silent(model)
    	Σ.m == 1 ? @variable(model, u) : @variable(model, u[1:Σ.m])
		@objective(model, Min, (1/2)u'u)
		for CBF in CBFs
			@constraint(model, CBF.∇h(x)*(Σ.f(x) + Σ.g(x)*u) >= -CBF.α(CBF(x)))
		end
		optimize!(model)
		
		return Σ.m == 1 ? value(u) : value.(u)
	end

	return CBFQP(control)
end

"""
	CBFQP(Σ::ControlAffineSystem, CBF::ControlBarrierFunction, U)

Construct a CBF-QP from a ControlAffineSystem and ControlBarrierFunction subject to
actuation constraints described by U.
"""
function CBFQP(Σ::ControlAffineSystem, CBF::ControlBarrierFunction, U)
	function control(x)
		model = Model(OSQP.Optimizer)
    	set_silent(model)
    	Σ.m == 1 ? @variable(model, u) : @variable(model, u[1:Σ.m])
		@objective(model, Min, (1/2)u'u)
		@constraint(model, cbf, CBF.∇h(x)*(Σ.f(x) + Σ.g(x)*u) >= -CBF.α(CBF(x)))
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

	return CBFQP(control)
end

"""
	CBFQP(Σ::ControlAffineSystem, CBFs::Vector{ControlBarrierFunction}, U)

Construct a CBF-QP from a ControlAffineSystem and ControlBarrierFunction subject to
actuation constraints described by U.
"""
function CBFQP(Σ::ControlAffineSystem, CBFs::Vector{ControlBarrierFunction}, U)
	function control(x)
		model = Model(OSQP.Optimizer)
    	set_silent(model)
    	Σ.m == 1 ? @variable(model, u) : @variable(model, u[1:Σ.m])
		@objective(model, Min, (1/2)u'u)
		for CBF in CBFs
			@constraint(model, CBF.∇h(x)*(Σ.f(x) + Σ.g(x)*u) >= -CBF.α(CBF(x)))
		end
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

	return CBFQP(control)
end

"""
	CBFQP(Σ::ControlAffineSystem, CBF::ControlBarrierFunction, κ::FeedbackPolicy)

Construct a CBF-QP as a safety filter for a nominal feedback control policy.
"""
function CBFQP(Σ::ControlAffineSystem, CBF::ControlBarrierFunction, κ::FeedbackPolicy)
	function control(x)
		model = Model(OSQP.Optimizer)
    	set_silent(model)
    	Σ.m == 1 ? @variable(model, u) : @variable(model, u[1:Σ.m])
		@objective(model, Min, (1/2)u'u - κ(x)'u)
		@constraint(model, cbf, CBF.∇h(x)*(Σ.f(x) + Σ.g(x)*u) >= -CBF.α(CBF(x)))
		optimize!(model)
		
		return Σ.m == 1 ? value(u) : value.(u)
	end

	return CBFQP(control)
end

"""
	CBFQP(Σ::ControlAffineSystem, CBFs::Vector{ControlBarrierFunction}, κ::FeedbackPolicy)

Construct a CBF-QP as a safety filter for a nominal feedback control policy.
"""
function CBFQP(Σ::ControlAffineSystem, CBFs::Vector{ControlBarrierFunction}, κ::FeedbackPolicy)
	function control(x)
		model = Model(OSQP.Optimizer)
    	set_silent(model)
    	Σ.m == 1 ? @variable(model, u) : @variable(model, u[1:Σ.m])
		@objective(model, Min, (1/2)u'u - κ(x)'u)
		for CBF in CBFs
			@constraint(model, CBF.∇h(x)*(Σ.f(x) + Σ.g(x)*u) >= -CBF.α(CBF(x)))
		end
		optimize!(model)
		
		return Σ.m == 1 ? value(u) : value.(u)
	end

	return CBFQP(control)
end

"""
	CBFQP(Σ::ControlAffineSystem, CBF::ControlBarrierFunction, κ::FeedbackPolicy, U)

Construct a CBF-QP as a safety filter for a nominal feedback control policy subject to
actuation constraints described by U.
"""
function CBFQP(Σ::ControlAffineSystem, CBF::ControlBarrierFunction, κ::FeedbackPolicy, U)
	function control(x)
		model = Model(OSQP.Optimizer)
    	set_silent(model)
    	Σ.m == 1 ? @variable(model, u) : @variable(model, u[1:Σ.m])
		@objective(model, Min, (1/2)u'u - κ(x)'u)
		@constraint(model, cbf, CBF.∇h(x)*(Σ.f(x) + Σ.g(x)*u) >= -CBF.α(CBF(x)))
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

	return CBFQP(control)
end

"""
	CBFQP(Σ::ControlAffineSystem, CBFs::Vector{ControlBarrierFunction}, κ::FeedbackPolicy, U)

Construct a CBF-QP as a safety filter for a nominal feedback control policy subject to
actuation constraints described by U.
"""
function CBFQP(Σ::ControlAffineSystem, CBFs::Vector{ControlBarrierFunction}, κ::FeedbackPolicy, U)
	function control(x)
		model = Model(OSQP.Optimizer)
    	set_silent(model)
    	Σ.m == 1 ? @variable(model, u) : @variable(model, u[1:Σ.m])
		@objective(model, Min, (1/2)u'u - κ(x)'u)
		for CBF in CBFs
			@constraint(model, CBF.∇h(x)*(Σ.f(x) + Σ.g(x)*u) >= -CBF.α(CBF(x)))
		end
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

	return CBFQP(control)
end

"""
	CBFQP(
		Σ::ControlAffineSystem, 
		CBF::ControlBarrierFunction, 
		CLF::ControlLyapunovFunction, 
		p::Float64=10e2
	)

Construct a CLF-CBF-QP, where p is the relaxation penalty on the CLF constraint.
"""
function CBFQP(
	Σ::ControlAffineSystem, 
	CBF::ControlBarrierFunction, 
	CLF::ControlLyapunovFunction, 
	p::Float64=10e2
)
	function control(x)
		model = Model(OSQP.Optimizer)
    	set_silent(model)
    	Σ.m == 1 ? @variable(model, u) : @variable(model, u[1:Σ.m])
		@variable(model, δ)
		@objective(model, Min, (1/2)u'u + p*δ^2)
		@constraint(model, cbf, CBF.∇h(x)*(Σ.f(x) + Σ.g(x)*u) >= -CBF.α(CBF(x)))
		@constraint(model, clf, CLF.∇V(x)*(Σ.f(x) + Σ.g(x)*u) <= -CLF.α(CLF(x)) + δ)
		optimize!(model)
		
		return Σ.m == 1 ? value(u) : value.(u)
	end

	return CBFQP(control)
end

"""
	CBFQP(
		Σ::ControlAffineSystem, 
		CBFs::Vector{ControlBarrierFunction},
		CLF::ControlLyapunovFunction, 
		p::Float64=10e2
	)

Construct a CLF-CBF-QP, where p is the relaxation penalty on the CLF constraint.
"""
function CBFQP(
	Σ::ControlAffineSystem, 
	CBFs::Vector{ControlBarrierFunction},
	CLF::ControlLyapunovFunction, 
	p::Float64=10e2
)
	function control(x)
		model = Model(OSQP.Optimizer)
    	set_silent(model)
    	Σ.m == 1 ? @variable(model, u) : @variable(model, u[1:Σ.m])
		@variable(model, δ)
		@objective(model, Min, (1/2)u'u + p*δ^2)
		for CBF in CBFs
			@constraint(model, CBF.∇h(x)*(Σ.f(x) + Σ.g(x)*u) >= -CBF.α(CBF(x)))
		end
		@constraint(model, clf, CLF.∇V(x)*(Σ.f(x) + Σ.g(x)*u) <= -CLF.α(CLF(x)) + δ)
		optimize!(model)
		
		return Σ.m == 1 ? value(u) : value.(u)
	end

	return CBFQP(control)
end

"""
	CBFQP(
		Σ::ControlAffineSystem, 
		CBF::ControlBarrierFunction, 
		CLF::ControlLyapunovFunction, 
		U,
		p::Float64=10e2
	)

Construct a CLF-CBF-QP, where p is the relaxation penalty on the CLF constraint.
"""
function CBFQP(
	Σ::ControlAffineSystem, 
	CBF::ControlBarrierFunction, 
	CLF::ControlLyapunovFunction, 
	U,
	p::Float64=10e2
)
	function control(x)
		model = Model(OSQP.Optimizer)
    	set_silent(model)
    	Σ.m == 1 ? @variable(model, u) : @variable(model, u[1:Σ.m])
		@variable(model, δ)
		@objective(model, Min, (1/2)u'u + p*δ^2)
		@constraint(model, cbf, CBF.∇h(x)*(Σ.f(x) + Σ.g(x)*u) >= -CBF.α(CBF(x)))
		@constraint(model, clf, CLF.∇V(x)*(Σ.f(x) + Σ.g(x)*u) <= -CLF.α(CLF(x)) + δ)
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

	return CBFQP(control)
end

"""
	CBFQP(
		Σ::ControlAffineSystem, 
		CBFs::Vector{ControlBarrierFunction}, 
		CLF::ControlLyapunovFunction, 
		U,
		p::Float64=10e2
	)

Construct a CLF-CBF-QP, where p is the relaxation penalty on the CLF constraint.
"""
function CBFQP(
	Σ::ControlAffineSystem, 
	CBFs::Vector{ControlBarrierFunction},
	CLF::ControlLyapunovFunction, 
	U,
	p::Float64=10e2
)
	function control(x)
		model = Model(OSQP.Optimizer)
    	set_silent(model)
    	Σ.m == 1 ? @variable(model, u) : @variable(model, u[1:Σ.m])
		@variable(model, δ)
		@objective(model, Min, (1/2)u'u + p*δ^2)
		for CBF in CBFs
			@constraint(model, CBF.∇h(x)*(Σ.f(x) + Σ.g(x)*u) >= -CBF.α(CBF(x)))
		end
		@constraint(model, clf, CLF.∇V(x)*(Σ.f(x) + Σ.g(x)*u) <= -CLF.α(CLF(x)) + δ)
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

	return CBFQP(control)
end