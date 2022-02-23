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

"""
    (sim::Simulation)((Σ::ControlAffineSystem, k::CBFQP, x0)

Simulate the closed-loop trajectory of a control affine system under a CBFQP from
initial condition x0.
"""
function (sim::Simulation)(Σ::ControlAffineSystem, k::CBFQP, x0)
	xs = Vector{typeof(x0)}(undef, length(sim))
	xs[1] = x0
	for i in 1:length(sim)-1
		t = sim.ts[i]
		x = xs[i]
		u = k(x)
		xs[i+1] = step(Σ, x, u, t, t + sim.dt)
	end

	return vec2mat(xs)
end

"""
	CBFQP(Σ::ControlAffineSystem, CBF::ControlBarrierFunction)
	CBFQP(Σ::ControlAffineSystem, CBF::ControlBarrierFunction, U)
	CBFQP(Σ::ControlAffineSystem, CBF::ControlBarrierFunction, H, F)
	CBFQP(Σ::ControlAffineSystem, CBF::ControlBarrierFunction, H, F, U)
    CBFQP(Σ::ControlAffineSystem, CBF::ControlBarrierFunction, κ::FeedbackPolicy)
	CBFQP(Σ::ControlAffineSystem, CBF::ControlBarrierFunction, κ::FeedbackPolicy, U)
	CBFQP(Σ::ControlAffineSystem, CBF::ControlBarrierFunction, κ::CLFQP)
	CBFQP(Σ::ControlAffineSystem, CBF::ControlBarrierFunction, κ::CLFQP, U)
	CBFQP(Σ::ControlAffineSystem, CBF::ControlBarrierFunction, CLF::ControlLyapunovFunction, p::Float64=10e2)
	CBFQP(Σ::ControlAffineSystem, CBF::ControlBarrierFunction, CLF::ControlLyapunovFunction, U, p::Float64=10e2)
	CBFQP(Σ::ControlAffineSystem, CBF::ControlBarrierFunction, CLF::ControlLyapunovFunction, H, F, p::Float64=10e2)
	CBFQP(Σ::ControlAffineSystem, CBF::ControlBarrierFunction, CLF::ControlLyapunovFunction, H, F, U, p::Float64=10e2)
	CBFQP(Σ::ControlAffineSystem, CBFs::Vector{ControlBarrierFunction})
	CBFQP(Σ::ControlAffineSystem, CBFs::Vector{ControlBarrierFunction}, U)
	CBFQP(Σ::ControlAffineSystem, CBFs::Vector{ControlBarrierFunction}, H, F)
	CBFQP(Σ::ControlAffineSystem, CBFs::Vector{ControlBarrierFunction}, H, F, U)
	CBFQP(Σ::ControlAffineSystem, CBFs::Vector{ControlBarrierFunction}, κ::FeedbackPolicy)
	CBFQP(Σ::ControlAffineSystem, CBFs::Vector{ControlBarrierFunction}, κ::FeedbackPolicy, U)
	CBFQP(Σ::ControlAffineSystem, CBFs::Vector{ControlBarrierFunction}, CLF::ControlLyapunovFunction, p::Float64=10e2)
	CBFQP(Σ::ControlAffineSystem, CBFs::Vector{ControlBarrierFunction}, CLF::ControlLyapunovFunction, U, p::Float64=10e2)
	CBFQP(Σ::ControlAffineSystem, CBFs::Vector{ControlBarrierFunction}, CLF::ControlLyapunovFunction, H, F, p::Float64=10e2)
	CBFQP(Σ::ControlAffineSystem, CBFs::Vector{ControlBarrierFunction}, CLF::ControlLyapunovFunction, H, F, U, p::Float64=10e2)

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

function CBFQP(Σ::ControlAffineSystem, CBF::ControlBarrierFunction, H, F)
	function control(x)
		model = Model(OSQP.Optimizer)
    	set_silent(model)
    	Σ.m == 1 ? @variable(model, u) : @variable(model, u[1:Σ.m])
		@objective(model, Min, (1/2)u'H(x)*u + F(x)'u)
		@constraint(model, cbf, CBF.∇h(x)*(Σ.f(x) + Σ.g(x)*u) >= -CBF.α(CBF(x)))
		optimize!(model)

		return Σ.m == 1 ? value(u) : value.(u)
	end

	return CBFQP(control)
end

function CBFQP(Σ::ControlAffineSystem, CBF::ControlBarrierFunction, H, F, U)
	function control(x)
		model = Model(OSQP.Optimizer)
    	set_silent(model)
    	Σ.m == 1 ? @variable(model, u) : @variable(model, u[1:Σ.m])
		@objective(model, Min, (1/2)u'H(x)*u + F(x)'u)
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

function CBFQP(Σ::ControlAffineSystem, CBF::ControlBarrierFunction, κ::CLFQP)
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

function CBFQP(Σ::ControlAffineSystem, CBF::ControlBarrierFunction, κ::CLFQP, U)
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

function CBFQP(Σ::ControlAffineSystem, CBF::ControlBarrierFunction, CLF::ControlLyapunovFunction, p::Float64=10e2)
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

function CBFQP(Σ::ControlAffineSystem, CBF::ControlBarrierFunction, CLF::ControlLyapunovFunction, U, p::Float64=10e2)
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

function CBFQP(Σ::ControlAffineSystem, CBF::ControlBarrierFunction, CLF::ControlLyapunovFunction, H, F, p::Float64=10e2)
	function control(x)
		model = Model(OSQP.Optimizer)
    	set_silent(model)
    	Σ.m == 1 ? @variable(model, u) : @variable(model, u[1:Σ.m])
		@variable(model, δ)
		@objective(model, Min, (1/2)u'H(x)*u + F(x)'u + p*δ^2)
		@constraint(model, cbf, CBF.∇h(x)*(Σ.f(x) + Σ.g(x)*u) >= -CBF.α(CBF(x)))
		@constraint(model, clf, CLF.∇V(x)*(Σ.f(x) + Σ.g(x)*u) <= -CLF.α(CLF(x)) + δ)
		optimize!(model)

		return Σ.m == 1 ? value(u) : value.(u)
	end

	return CBFQP(control)
end

function CBFQP(Σ::ControlAffineSystem, CBF::ControlBarrierFunction, CLF::ControlLyapunovFunction, H, F, U, p::Float64=10e2)
	function control(x)
		model = Model(OSQP.Optimizer)
    	set_silent(model)
    	Σ.m == 1 ? @variable(model, u) : @variable(model, u[1:Σ.m])
		@variable(model, δ)
		@objective(model, Min, (1/2)u'H(x)*u + F(x)'u + p*δ^2)
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

function CBFQP(Σ::ControlAffineSystem, CBFs::Vector{ControlBarrierFunction}, H, F)
	function control(x)
		model = Model(OSQP.Optimizer)
    	set_silent(model)
    	Σ.m == 1 ? @variable(model, u) : @variable(model, u[1:Σ.m])
		@objective(model, Min, (1/2)u'H(x)*u + F(x)'u)
		for CBF in CBFs
			@constraint(model, CBF.∇h(x)*(Σ.f(x) + Σ.g(x)*u) >= -CBF.α(CBF(x)))
		end
		optimize!(model)

		return Σ.m == 1 ? value(u) : value.(u)
	end

	return CBFQP(control)
end

function CBFQP(Σ::ControlAffineSystem, CBFs::Vector{ControlBarrierFunction}, H, F, U)
	function control(x)
		model = Model(OSQP.Optimizer)
    	set_silent(model)
    	Σ.m == 1 ? @variable(model, u) : @variable(model, u[1:Σ.m])
		@objective(model, Min, (1/2)u'H(x)*u + F(x)'u)
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

function CBFQP(Σ::ControlAffineSystem, CBFs::Vector{ControlBarrierFunction}, CLF::ControlLyapunovFunction, p::Float64=10e2)
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

function CBFQP(Σ::ControlAffineSystem, CBFs::Vector{ControlBarrierFunction}, CLF::ControlLyapunovFunction, U, p::Float64=10e2)
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

function CBFQP(Σ::ControlAffineSystem, CBFs::Vector{ControlBarrierFunction}, CLF::ControlLyapunovFunction, H, F, p::Float64=10e2)
	function control(x)
		model = Model(OSQP.Optimizer)
    	set_silent(model)
    	Σ.m == 1 ? @variable(model, u) : @variable(model, u[1:Σ.m])
		@variable(model, δ)
		@objective(model, Min, (1/2)u'H(x)*u + F(x)'u + p*δ^2)
		for CBF in CBFs
			@constraint(model, CBF.∇h(x)*(Σ.f(x) + Σ.g(x)*u) >= -CBF.α(CBF(x)))
		end
		@constraint(model, clf, CLF.∇V(x)*(Σ.f(x) + Σ.g(x)*u) <= -CLF.α(CLF(x)) + δ)
		optimize!(model)

		return Σ.m == 1 ? value(u) : value.(u)
	end

	return CBFQP(control)
end

function CBFQP(Σ::ControlAffineSystem, CBFs::Vector{ControlBarrierFunction}, CLF::ControlLyapunovFunction, H, F, U, p::Float64=10e2)
	function control(x)
		model = Model(OSQP.Optimizer)
    	set_silent(model)
    	Σ.m == 1 ? @variable(model, u) : @variable(model, u[1:Σ.m])
		@variable(model, δ)
		@objective(model, Min, (1/2)u'H(x)*u + F(x)'u + p*δ^2)
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
