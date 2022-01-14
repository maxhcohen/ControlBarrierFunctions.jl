"""
	CBFQP(Σ::ControlAffineSystem, HOCBF::HighOrderCBF)
	CBFQP(Σ::ControlAffineSystem, HOCBF::HighOrderCBF, U)
	CBFQP(Σ::ControlAffineSystem, HOCBF::HighOrderCBF, κ::FeedbackPolicy)
	CBFQP(Σ::ControlAffineSystem, HOCBF::HighOrderCBF, κ::FeedbackPolicy, U)
	CBFQP(Σ::ControlAffineSystem, HOCBF::HighOrderCBF, CLF::ControlLyapunovFunction, p::Float64=10e2)
	CBFQP(Σ::ControlAffineSystem, HOCBF::HighOrderCBF, CLF::ControlLyapunovFunction, U, p::Float64=10e2)
	CBFQP(Σ::ControlAffineSystem, HOCBFs::Vector{HighOrderCBF})
	CBFQP(Σ::ControlAffineSystem, HOCBFs::Vector{HighOrderCBF}, U)
	CBFQP(Σ::ControlAffineSystem, HOCBFs::Vector{HighOrderCBF}, κ::FeedbackPolicy)
	CBFQP(Σ::ControlAffineSystem, HOCBFs::Vector{HighOrderCBF}, κ::FeedbackPolicy, U)
	CBFQP(Σ::ControlAffineSystem, HOCBFs::Vector{HighOrderCBF}, CLF::ControlLyapunovFunction, p::Float64=10e2)
	CBFQP(Σ::ControlAffineSystem, HOCBFs::Vector{HighOrderCBF}, CLF::ControlLyapunovFunction, U, p::Float64=10e2)

Construct a High Order Control Barrier Function Quadratic Program.
"""
function CBFQP(Σ::ControlAffineSystem, HOCBF::HighOrderCBF)
	function control(x)
		model = Model(OSQP.Optimizer)
    	set_silent(model)
    	Σ.m == 1 ? @variable(model, u) : @variable(model, u[1:Σ.m])
		@objective(model, Min, (1/2)u'u)
		@constraint(model, cbf, HOCBF.∇ψ(x)*(Σ.f(x) + Σ.g(x)*u) >= -HOCBF.α(HOCBF(x)))
		optimize!(model)
		
		return Σ.m == 1 ? value(u) : value.(u)
	end

	return CBFQP(control)
end

function CBFQP(Σ::ControlAffineSystem, HOCBF::HighOrderCBF, U)
	function control(x)
		model = Model(OSQP.Optimizer)
    	set_silent(model)
    	Σ.m == 1 ? @variable(model, u) : @variable(model, u[1:Σ.ms])
		@objective(model, Min, (1/2)u'u)
		@constraint(model, cbf, HOCBF.∇ψ(x)*(Σ.f(x) + Σ.g(x)*u) >= -HOCBF.α(HOCBF(x)))
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

function CBFQP(Σ::ControlAffineSystem, HOCBF::HighOrderCBF, κ::FeedbackPolicy)
	function control(x)
		model = Model(OSQP.Optimizer)
    	set_silent(model)
    	Σ.m == 1 ? @variable(model, u) : @variable(model, u[1:Σ.m])
		@objective(model, Min, (1/2)u'u - κ(x)'u)
		@constraint(model, cbf, HOCBF.∇ψ(x)*(Σ.f(x) + Σ.g(x)*u) >= -HOCBF.α(HOCBF(x)))
		optimize!(model)
		
		return Σ.m == 1 ? value(u) : value.(u)
	end

	return CBFQP(control)
end

function CBFQP(Σ::ControlAffineSystem, HOCBF::HighOrderCBF, κ::FeedbackPolicy, U)
	function control(x)
		model = Model(OSQP.Optimizer)
    	set_silent(model)
    	Σ.m == 1 ? @variable(model, u) : @variable(model, u[1:Σ.m])
		@objective(model, Min, (1/2)u'u - κ(x)'u)
		@constraint(model, cbf, HOCBF.∇ψ(x)*(Σ.f(x) + Σ.g(x)*u) >= -HOCBF.α(HOCBF(x)))
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

function CBFQP(Σ::ControlAffineSystem, HOCBF::HighOrderCBF, CLF::ControlLyapunovFunction, p::Float64=10e2)
	function control(x)
		model = Model(OSQP.Optimizer)
    	set_silent(model)
    	Σ.m == 1 ? @variable(model, u) : @variable(model, u[1:Σ.m])
		@variable(model, δ)
		@objective(model, Min, (1/2)u'u + p*δ^2)
		@constraint(model, cbf, HOCBF.∇ψ(x)*(Σ.f(x) + Σ.g(x)*u) >= -HOCBF.α(HOCBF(x)))
		@constraint(model, clf, CLF.∇V(x)*(Σ.f(x) + Σ.g(x)*u) <= -CLF.α(CLF(x)) + δ)
		optimize!(model)
		
		return Σ.m == 1 ? value(u) : value.(u)
	end

	return CBFQP(control)
end

function CBFQP(Σ::ControlAffineSystem, HOCBF::HighOrderCBF, CLF::ControlLyapunovFunction, U, p::Float64=10e2)
	function control(x)
		model = Model(OSQP.Optimizer)
    	set_silent(model)
    	Σ.m == 1 ? @variable(model, u) : @variable(model, u[1:Σ.m])
		@variable(model, δ)
		@objective(model, Min, (1/2)u'u + p*δ^2)
		@constraint(model, cbf, HOCBF.∇ψ(x)*(Σ.f(x) + Σ.g(x)*u) >= -HOCBF.α(HOCBF(x)))
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

function CBFQP(Σ::ControlAffineSystem, HOCBFs::Vector{HighOrderCBF})
	function control(x)
		model = Model(OSQP.Optimizer)
    	set_silent(model)
    	Σ.m == 1 ? @variable(model, u) : @variable(model, u[1:Σ.m])
		@objective(model, Min, (1/2)u'u)
		for HOCBF in HOCBFs
			@constraint(model, HOCBF.∇ψ(x)*(Σ.f(x) + Σ.g(x)*u) >= -HOCBF.α(HOCBF(x)))
		end
		optimize!(model)
		
		return Σ.m == 1 ? value(u) : value.(u)
	end

	return CBFQP(control)
end

function CBFQP(Σ::ControlAffineSystem, HOCBFs::Vector{HighOrderCBF}, U)
	function control(x)
		model = Model(OSQP.Optimizer)
    	set_silent(model)
    	Σ.m == 1 ? @variable(model, u) : @variable(model, u[1:Σ.m])
		@objective(model, Min, (1/2)u'u)
		for HOCBF in HOCBFs
			@constraint(model, HOCBF.∇ψ(x)*(Σ.f(x) + Σ.g(x)*u) >= -HOCBF.α(HOCBF(x)))
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

function CBFQP(Σ::ControlAffineSystem, HOCBFs::Vector{HighOrderCBF}, κ::FeedbackPolicy)
	function control(x)
		model = Model(OSQP.Optimizer)
    	set_silent(model)
    	Σ.m == 1 ? @variable(model, u) : @variable(model, u[1:Σ.m])
		@objective(model, Min, (1/2)u'u - κ(x)'u)
		for HOCBF in HOCBFs
			@constraint(model, HOCBF.∇ψ(x)*(Σ.f(x) + Σ.g(x)*u) >= -HOCBF.α(HOCBF(x)))
		end
		optimize!(model)
		
		return Σ.m == 1 ? value(u) : value.(u)
	end

	return CBFQP(control)
end

function CBFQP(Σ::ControlAffineSystem, HOCBFs::Vector{HighOrderCBF}, κ::FeedbackPolicy, U)
	function control(x)
		model = Model(OSQP.Optimizer)
    	set_silent(model)
    	Σ.m == 1 ? @variable(model, u) : @variable(model, u[1:Σ.m])
		@objective(model, Min, (1/2)u'u - κ(x)'u)
		for HOCBF in HOCBFs
			@constraint(model, HOCBF.∇ψ(x)*(Σ.f(x) + Σ.g(x)*u) >= -HOCBF.α(HOCBF(x)))
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

function CBFQP(Σ::ControlAffineSystem, HOCBFs::Vector{HighOrderCBF}, CLF::ControlLyapunovFunction, p::Float64=10e2)
	function control(x)
		model = Model(OSQP.Optimizer)
    	set_silent(model)
    	Σ.m == 1 ? @variable(model, u) : @variable(model, u[1:Σ.m])
		@variable(model, δ)
		@objective(model, Min, (1/2)u'u + p*δ^2)
		for HOCBF in HOCBFs
			@constraint(model, HOCBF.∇ψ(x)*(Σ.f(x) + Σ.g(x)*u) >= -HOCBF.α(HOCBF(x)))
		end
		@constraint(model, clf, CLF.∇V(x)*(Σ.f(x) + Σ.g(x)*u) <= -CLF.α(CLF(x)) + δ)
		optimize!(model)
		
		return Σ.m == 1 ? value(u) : value.(u)
	end

	return CBFQP(control)
end

function CBFQP(Σ::ControlAffineSystem, HOCBFs::Vector{HighOrderCBF}, CLF::ControlLyapunovFunction, U, p::Float64=10e2)
	function control(x)
		model = Model(OSQP.Optimizer)
    	set_silent(model)
    	Σ.m == 1 ? @variable(model, u) : @variable(model, u[1:Σ.m])
		@variable(model, δ)
		@objective(model, Min, (1/2)u'u + p*δ^2)
		for HOCBF in HOCBFs
			@constraint(model, HOCBF.∇ψ(x)*(Σ.f(x) + Σ.g(x)*u) >= -HOCBF.α(HOCBF(x)))
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