"""
	SmoothSafetyFilter <: SafetyFilter

Smooth controller that approximates CBF-QP arbitrarily closely.

# Fields
- `formula::String` : string indicating formula used in smooth safety filter
- `σ::Float64` : smoothing parameter
- `k::Function` : function that computes safe control actions
"""
struct SmoothSafetyFilter <: SafetyFilter
	formula::String
	σ::Float64
	k::Function
end

##### SmoothSafetyFilter functors #####

"""
	(k::SmoothSafetyFilter)(x)

Functors for evaluating smooth safety filter
"""
(k::SmoothSafetyFilter)(x) = k.k(x)
(k::SmoothSafetyFilter)(x, t) = k.k(x, t)

##### SmoothSafetyFilter constructors #####

"""
	SmoothSafetyFilter(cbf::ControlBarrierFunction, Σ::ControlAffineSystem, kd::Function)

Construct an SmoothSafetyFilter from a cbf and a desired controller.

# Keyword Arguments 
- `formula::String` : formula for computing smooth safety filter
- `σ::Float` : smoothing parameter, where larger `σ` is smoother
"""
function SmoothSafetyFilter(
	cbf::ControlBarrierFunction,
	Σ::ControlAffineSystem,
	kd::Function;
	formula = "half sontag",
	σ = 0.1,
)
	try
		kd(Σ.n == 1 ? rand() : rand(Σ.n), 0.0)
	catch e
		if isa(e, MethodError)
			a(x) = cbf.Lfh(x) + cbf.Lgh(x) * kd(x) + cbf.α(cbf(x))
			if formula == "half sontag"
				return SmoothSafetyFilter(
					formula,
					σ,
					x -> kd(x) + λHalfSontag(a(x), norm(cbf.Lgh(x))^2, σ) * cbf.Lgh(x)',
				)
			elseif formula == "sontag"
				return SmoothSafetyFilter(
					formula,
					σ,
					x -> kd(x) + λSontag(a(x), norm(cbf.Lgh(x))^2, σ) * cbf.Lgh(x)',
				)
			elseif formula == "softplus"
				return SmoothSafetyFilter(
					formula,
					σ,
					x -> kd(x) + λSoftplus(a(x), norm(cbf.Lgh(x))^2, σ) * cbf.Lgh(x)',
				)
			else
				@warn "No valid formula provided, defaulting to Half Sontag formula."
				return SmoothSafetyFilter(
					formula,
					σ,
					x -> kd(x) + λHalfSontag(a(x), norm(cbf.Lgh(x))^2, σ) * cbf.Lgh(x)',
				)
			end
		else
			return e
		end
	else
		a(x, t) = cbf.Lfh(x) + cbf.Lgh(x) * kd(x, t) + cbf.α(cbf(x))
		if formula == "half sontag"
			return SmoothSafetyFilter(
				formula,
				σ,
				(x, t) ->
					kd(x, t) + λHalfSontag(a(x, t), norm(cbf.Lgh(x))^2, σ) * cbf.Lgh(x)',
			)
		elseif formula == "sontag"
			return SmoothSafetyFilter(
				formula,
				σ,
				(x, t) -> kd(x, t) + λSontag(a(x, t), norm(cbf.Lgh(x))^2, σ) * cbf.Lgh(x)',
			)
		elseif formula == "softplus"
			return SmoothSafetyFilter(
				formula,
				σ,
				(x, t) ->
					kd(x, t) + λSoftplus(a(x, t), norm(cbf.Lgh(x))^2, σ) * cbf.Lgh(x)',
			)
		else
			@warn "No valid formula provided, defaulting to Half Sontag formula."
			return SmoothSafetyFilter(
				formula,
				σ,
				(x, t) ->
					kd(x, t) + λHalfSontag(a(x, t), norm(cbf.Lgh(x))^2, σ) * cbf.Lgh(x)',
			)
		end
	end
end

"""
	SmoothSafetyFilter(hocbf::HighOrderCBF, Σ::ControlAffineSystem, kd::Function)

Construct a smooth safety filter from an HOCBF.

# Keyword Arguments 
- `formula::String` : formula for computing smooth safety filter
- `σ::Float` : smoothing parameter, where larger `σ` is smoother
"""
function SmoothSafetyFilter(
	hocbf::HighOrderCBF,
	Σ::ControlAffineSystem,
	kd::Function;
	formula = "half sontag",
	σ = 0.1,
)
	return SmoothSafetyFilter(hocbf.cbf, Σ, kd; formula = formula, σ = σ)
end

# Smooth universal formulas, we always default to half Sontag
λSontag(a, b, σ) = b == 0.0 ? 0.0 : (-a + sqrt(a^2 + σ * b^2)) / b
λHalfSontag(a, b, σ) = 0.5 * λSontag(a, b, σ)
λSoftplus(a, b, σ) = b ≤ 0.0 ? 0.0 : σ * log(1.0 + exp(-a / (b * σ)))
