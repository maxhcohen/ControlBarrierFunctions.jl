"""
	TunableQPSafetyFilter <: SafetyFilter

Controller that solves a control barrier function-based quadratic program (CBF-QP) with tunable class K functions.

Uses OSQP to solve the corresponding QP.

# Fields
- `k::Function` : function that computes safe control actions
"""
struct TunableQPSafetyFilter <: SafetyFilter
	k::Function
end

##### TunableQPSafetyFilter functors #####

"""
	(k::TunableQPSafetyFilter)(x)

Functors for evaluating QP-based safety filter
"""
(k::TunableQPSafetyFilter)(x) = k.k(x)
(k::TunableQPSafetyFilter)(x, t) = k.k(x, t)

##### TunableQPSafetyFilter constructors #####

"""
	TunableQPSafetyFilter(cbfs::Vector{ControlBarrierFunction}, Σ::ControlAffineSystem, kd::Function; tunable=false)

Construct an TunableQPSafetyFilter from a cbf and a desired controller.
"""
function TunableQPSafetyFilter(
	cbfs::Vector{ControlBarrierFunction}, Σ::ControlAffineSystem, kd::Function,
)
	try
		kd(Σ.n == 1 ? rand() : rand(Σ.n), 0.0) # See if desired controller is time-varying
	catch e
		if isa(e, MethodError) # If controller is not time-varying 
			return TunableQPSafetyFilter(x -> solve_tunable_cbf_qp(x, Σ, cbfs, kd))
		else
			return e
		end
	else # If controller is time-varying
		return TunableQPSafetyFilter(
			(x, t) -> solve_time_varying_tunable_cbf_qp(x, t, Σ, cbfs, kd),
		)
	end
end

"""
	TunableQPSafetyFilter(cbf::ControlBarrierFunction, Σ::ControlAffineSystem, kd::Function; kwargs...)

Add ability to pass in single CBF.
"""
TunableQPSafetyFilter(cbf::ControlBarrierFunction, Σ::ControlAffineSystem, kd::Function) =
	TunableQPSafetyFilter([cbf], Σ, kd)

"""
	TunableQPSafetyFilter(hocbf::HighOrderCBF, Σ::ControlAffineSystem, kd::Function)}

Construct a TunableQPSafetyFilter from an HOCBF and desired controller.
"""
function TunableQPSafetyFilter(hocbf::HighOrderCBF, Σ::ControlAffineSystem, kd::Function)
	return TunableQPSafetyFilter(hocbf.cbf, Σ, kd)
end

"""
	TunableQPSafetyFilter(hocbfs::Vector{HighOrderCBF} Σ::ControlAffineSystem, kd::Function)}

Construct a TunableQPSafetyFilter from a multiple HOCBFs and desired controller.
"""
function TunableQPSafetyFilter(hocbfs::Vector{HighOrderCBF}, Σ::ControlAffineSystem, kd::Function)
	# Create vector of cbfs
	cbfs = [hocbf.cbf for hocbf in hocbfs]

	return TunableQPSafetyFilter(cbfs, Σ, kd)
end

##### Funcions for solving QPs #####

"""
	solve_tunable_cbf_qp(x, Σ::ControlAffineSystem, cbfs::Vector{ControlBarrierFunction}, kd::Function)

Solve CBF-QP where coefficients on extended class K functions are decision variables
"""
function solve_tunable_cbf_qp(
	x, Σ::ControlAffineSystem, cbfs::Vector{ControlBarrierFunction}, kd::Function,
)
	model = Model(OSQP.Optimizer)
	set_silent(model)
	u = Σ.m == 1 ? @variable(model, u) : @variable(model, u[1:(Σ.m)])
	@variable(model, αs[1:length(cbfs)])
	@objective(
		model, Min, 0.5 * (u - kd(x))' * (u - kd(x)) + sum([(α - 1.0)^2 for α in αs])
	)
	for (cbf, α) in zip(cbfs, αs)
		@constraint(model, cbf.Lfh(x) + cbf.Lgh(x) * u ≥ -α * cbf.α(cbf(x)))
	end
	optimize!(model)

	return Σ.m == 1 ? value(u) : value.(u)
end

"""
	solve_time_varying_tunable_cbf_qp(x, t, Σ::ControlAffineSystem, cbfs::Vector{ControlBarrierFunction}, kd::Function)

Solve CBF-QP where coefficients on extended class K functions are decision variables and desired controller depends on time
"""
function solve_time_varying_tunable_cbf_qp(
	x, t, Σ::ControlAffineSystem, cbfs::Vector{ControlBarrierFunction}, kd::Function,
)
	model = Model(OSQP.Optimizer)
	set_silent(model)
	u = Σ.m == 1 ? @variable(model, u) : @variable(model, u[1:(Σ.m)])
	@variable(model, αs[1:length(cbfs)])
	@objective(
		model, Min, 0.5 * (u - kd(x, t))' * (u - kd(x, t)) + sum([(α - 1.0)^2 for α in αs])
	)
	for (cbf, α) in zip(cbfs, αs)
		@constraint(model, cbf.Lfh(x) + cbf.Lgh(x) * u ≥ -α * cbf.α(cbf(x)))
	end
	optimize!(model)

	return Σ.m == 1 ? value(u) : value.(u)
end
