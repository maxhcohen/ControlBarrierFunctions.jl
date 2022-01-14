"""
	CBFQP(Σ::ControlAffineSystem, HOCBF::HighOrderCBF)

Construct a CBF-QP from a ControlAffineSystem and ControlBarrierFunction.
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