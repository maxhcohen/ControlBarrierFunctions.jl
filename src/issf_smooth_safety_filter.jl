"""
    ISSfSmoothSafetyFilter <: SafetyFilter

Smooth controller that approximates an input-to-state safe (ISSf) CBF-QP arbitrarily closely.

# Fields
- `formula::String` : string indicating formula used in smooth safety filter
- `σ::Float64` : smoothing parameter
- `k::Function` : function that computes safe control actions
- `ε::Float64` : Issf robustness parameter for matched uncertainties
"""
struct ISSfSmoothSafetyFilter <: SafetyFilter
    formula::String
    σ::Float64
    k::Function
    ε::Float64
end

"""
    (k::ISSfSmoothSafetyFilter)(x)

Functors for evaluating smooth safety filter
"""
(k::ISSfSmoothSafetyFilter)(x) = k.k(x)
(k::ISSfSmoothSafetyFilter)(x, t) = k.k(x, t)

"""
    ISSfSmoothSafetyFilter(cbf::ControlBarrierFunction, Σ::ControlAffineSystem, kd::Function)

Construct an ISSfSmoothSafetyFilter from a cbf and a desired controller.
"""
function ISSfSmoothSafetyFilter(
    cbf::ControlBarrierFunction,
    Σ::ControlAffineSystem,
    kd::Function,
    ε::Float64;
    formula="half sontag",
    σ=0.1,
)
    @assert(ε > 0.0, "ε must be a positive number!")
    @assert(σ > 0.0, "σ must be a positive number!")
    try
        kd(Σ.n == 1 ? rand() : rand(Σ.n), 0.0)
    catch e
        if isa(e, MethodError)
            function a(x)
                return cbf.Lfh(x) + cbf.Lgh(x) * kd(x) + cbf.α(cbf(x)) -
                       (1 / ε) * norm(cbf.Lgh(x))^2
            end
            b(x) = norm(cbf.Lgh(x))^2
            if formula == "half sontag"
                return ISSfSmoothSafetyFilter(
                    formula,
                    σ,
                    x -> kd(x) + λHalfSontag(a(x), norm(cbf.Lgh(x))^2, σ) * cbf.Lgh(x)',
                    ε,
                )
            elseif formula == "sontag"
                return ISSfSmoothSafetyFilter(
                    formula,
                    σ,
                    x -> kd(x) + λSontag(a(x), norm(cbf.Lgh(x))^2, σ) * cbf.Lgh(x)',
                    ε,
                )
            elseif formula == "softplus"
                return ISSfSmoothSafetyFilter(
                    formula,
                    σ,
                    x -> kd(x) + λSoftplus(a(x), norm(cbf.Lgh(x))^2, σ) * cbf.Lgh(x)',
                    ε,
                )
            else
                @warn "No valid formula provided, defaulting to Half Sontag formula."
                return ISSfSmoothSafetyFilter(
                    formula,
                    σ,
                    x -> kd(x) + λHalfSontag(a(x), norm(cbf.Lgh(x))^2, σ) * cbf.Lgh(x)',
                    ε,
                )
            end
        else
            return e
        end
    else
        function a(x, t)
            return cbf.Lfh(x) + cbf.Lgh(x) * kd(x, t) + cbf.α(cbf(x)) -
                   (1 / ε) * norm(cbf.Lgh(x))^2
        end
        if formula == "half sontag"
            return ISSfSmoothSafetyFilter(
                formula,
                σ,
                (x, t) ->
                    kd(x, t) + λHalfSontag(a(x, t), norm(cbf.Lgh(x))^2, σ) * cbf.Lgh(x)',
                ε,
            )
        elseif formula == "sontag"
            return ISSfSmoothSafetyFilter(
                formula,
                σ,
                (x, t) -> kd(x, t) + λSontag(a(x, t), norm(cbf.Lgh(x))^2, σ) * cbf.Lgh(x)',
                ε,
            )
        elseif formula == "softplus"
            return ISSfSmoothSafetyFilter(
                formula,
                σ,
                (x, t) ->
                    kd(x, t) + λSoftplus(a(x, t), norm(cbf.Lgh(x))^2, σ) * cbf.Lgh(x)',
                ε,
            )
        else
            @warn "No valid formula provided, defaulting to Half Sontag formula."
            return ISSfSmoothSafetyFilter(
                formula,
                σ,
                (x, t) ->
                    kd(x, t) + λHalfSontag(a(x, t), norm(cbf.Lgh(x))^2, σ) * cbf.Lgh(x)',
                ε,
            )
        end
    end
end
