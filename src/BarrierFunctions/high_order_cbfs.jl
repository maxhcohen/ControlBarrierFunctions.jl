"""
    HighOrderCBF

High Order Control barrier function for a control affine system.
# Fields
- `CBF`: CBF used to derive HOCBF
- `ψ`: highest order derivative of CBF
- `Lfψ`: Lie derivative of ψ along f
- `Lgψ`: Lie derivative of ψ along g
"""
struct HighOrderCBF <: BarrierFunction
    ψ
    ∇ψ
	α
end

"""
	HighOrderCBF(CBF::ControlBarrierFunction, Σ::ControlAffineSystem, degree::Int)

Construct a High Order Control Barrier Function of a specified degree.
"""
function HighOrderCBF(CBF::ControlBarrierFunction, Σ::ControlAffineSystem, degree::Int)
	if degree == 2
		return second_order_hocbf(CBF, Σ)
	end
end

"""
	(HOCBF::HighOrderCBF)(x)

Evaluate HOCBF at state x.
"""
function (HOCBF::HighOrderCBF)(x)
	return HOCBF.ψ(x)
end

"""
	second_order_hocbf(CBF::ControlBarrierFunction, Σ::ControlAffineSystem)

Construct a HOCBF of relative degree 2.
"""
function second_order_hocbf(CBF::ControlBarrierFunction, Σ::ControlAffineSystem)
	dfdx(x) = jacobian(Σ.f, x)[1]
	dh2dx2(x) = hessian(CBF.h, x)
	dα(s) = gradient(CBF.α, s)[1]
	ψ1(x) = CBF.∇h(x)*Σ.f(x) + CBF.α(CBF(x))
	∇ψ1(x) = CBF.∇h(x)*dfdx(x) + (dh2dx2(x)*Σ.f(x))' + dα(CBF(x))*CBF.∇h(x)

	return HighOrderCBF(ψ1, ∇ψ1, CBF.α)
end