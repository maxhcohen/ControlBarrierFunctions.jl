"""
    SafetyFilter

Abstract safety filter.

At a bare minumum, each subtype of SafetyFilter should have a field `k::Function` that computes inputs.
"""
abstract type SafetyFilter end

"""
    simulate(Σ::ControlAffineSystem, k::SafetyFilter, x0, T::Float64)

Simulate under an explicit safety filter 
"""
simulate(Σ::ControlAffineSystem, k::SafetyFilter, x0, T::Float64) = simulate(Σ, k.k, x0, T)